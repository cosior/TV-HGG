%Code for prediction of radial (popularity) or expected degree trajectories.
%Uses simple heuristics for tuning the model parameters for predictions.

clear;
clc;
close all;

warning('off');

%Uncomment the code below based on the network. Adjust parameters as needed.
%If iters > 1 we will be averaging over iters trajectories and
%computing a confidence region based on the value of nstd.

%US Air.
th1=1;
th2=10000;
ksigma=8;
alpha=0.01;
iters=10;
maxlag=10;
nstd=2;

%BTC.
% th1=200;
% th2=1492;
% ksigma=8;
% alpha=0.01;
% iters=100;
% maxlag=10;
% nstd=2;

%PGP.
% th1=1;
% th2=1500;
% ksigma=8;
% alpha=0.01;
% iters=100;
% maxlag=10;
% nstd=2;

%IPv6.
% th1=150;
% th2=450;
% ksigma=8;
% alpha=0.1; 
% iters=100;
% maxlag=10;
% nstd=1;

%arXiv.
% th1=2800;
% th2=6776;
% ksigma=8;
% alpha=0.01;
% iters=10;
% maxlag=20;
% nstd=1;

%Input dir (trajectories location, change accordingly);
% Files=dir('./USAir/RADIAL_TRAJECTORIES/*.txt');
Files=dir('./USAir/KAPPA_TRAJECTORIES/*.txt');

%Output dir (where plots will be stored, change accordingly).
% mydir="./USAir/RADIAL_PREDS/";
mydir="./USAir/KAPPA_PREDS/";
mkdir(mydir);

count=0;

for fn = 1:length(Files)

   %If you want to check a specific node.
   %if Files(fn).name ~= "Node_CLT.txt"
   %       continue;
   %end

   close all;

   clear data X tX;
   filename=strcat(Files(fn).folder, "/", Files(fn).name);
   data=readtable(filename);
   j=1;
   for i=1 :length(data.value)
        if (data.x_snapId(i) >= th1 && data.x_snapId(i) <= th2)
            X(j) = data.value(i);
            tX(j) = data.x_snapId(i);
            j=j+1;
        end
    end

   if (j == 1)
       continue;
   end

   n = length(X);

   %Ignore trajectories with less than 300 points.
   if n < 300
      continue;
   end
   
   %80% of the points are used for training.
   training=floor(0.8*n);

   clear velocities;
   velocities = zeros(n, 1);

   velocities(2:end)=diff(X);

   %Test for velocities' trend-stationarity. 
   if (adftest(velocities)~=1)
       printf("Velocity series likely not stationary for node %s\n", Files(fn).name);
       break;
   end
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%mu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear mu;
   mu=zeros(n,1);      

   %Exponentially weighted moving average. Calculations done up to training.
   for i = 2 : training
       mu(i) = alpha*velocities(i)+(1-alpha)*mu(i-1);
   end

   %mu for predictions.
   %simple heuristic: use the mean mu up to the training 
   %other heuristics could work better depending on the trajectory/network
   mu(training+1:end) = mean(mu(2:training));  

   %more general heuristic: use the mean mu over a previous window of length x
   %(x could be optimized).
   %x = min(training-2, 100);
   %mu(training+1:end) = mean(mu(training-x:training));

   %%%%%%%%%%%%%%%%%%%%%%%%%%hurst index%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear H H_estimates;

   lags = 1: maxlag;    
   q_orders = [1 2 3]; 

   H_estimates = zeros(1, length(q_orders));

   for m = 1 : length(q_orders)

        q = q_orders(m);
        clear M;
        M = zeros(1, length(lags));

        for l = 1:length(lags)
            tau = lags(l);
            clear increments;
            %Compute H from the trajectory up to the training.
            increments = abs(X(tau+1:training) - X(1:training-tau));  
            M(l) = mean(increments.^q);  %q-th order absolute moment
        end

        %Linear regression of log(M) against log(lags)
        p = polyfit(log(lags), log(M), 1);

        %The slope of the line gives q*H, so H = slope / q
        H_estimates(m) = max(p(1) / q, 10^(-6));
   end

   H=mean(H_estimates);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sigma%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear S sigma;
   sigma = zeros(n,1);

   %compute sigmas up to the training.
   for j = 1 : training

         a = max(j-ksigma/2, 1);
         b = min(j+ksigma/2, training);

         clear window_data;

         window_data=velocities(a:b);

         %Using 1st moment.  
         c = sqrt(2)/gamma(0.5);
         S = mean(abs(window_data));
         sigma(j) = S/c;
         sigma(j) = sigma(j)/n^(-H);

    end

    %sigma for predictions.
    %simple heuristic: use the average sigma up to the training.
    sigma(training+1:end) = mean(sigma(1:training));

    %%%%%%%%%%%%%%%%%%%%%%%%%mBm sims%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    clear run;

    t=1;
 
    while t <= iters

           clear Xsim; 
           
           Xsim = mBm(n, H, [ ], mu, X(1), sigma);
           %Expected degrees or radial coordinates cannot be negative.
           run(t,:) = max(Xsim, 0);
         
           t = t+1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear avepos stdpos;
   
    i=1:n;

    if (iters > 1)
        avepos(i) = mean(run(:,i));
        stdpos(i) = sqrt(var(run(:,i)));
    else 
        avepos(i) = run(1,i);
        stdpos(j) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Visible', 'off');

    hold on;
    
    if iters > 1

        clear curve1 curve2;

        %Consider an nstd region around the average. 
        curve1(i) = avepos(i) + nstd*stdpos(i);
        curve2(i) = avepos(i) - nstd*stdpos(i);

        clear i2;
        i2 = [tX(i), fliplr(tX(i))];
        %curve2 cannot be negative.
        inBetween = [curve1, max(fliplr(curve2), 0)];
        fill(i2, inBetween, 'y', 'EdgeColor', 'none', 'DisplayName', 'std region');

    end

    plot(tX, X, 'Color', 'blue', 'DisplayName','real', LineWidth=0.5);

    if iters > 1
        plot(tX, avepos, 'Color', 'red', 'DisplayName','synthetic average', LineWidth=0.5);
    else
        plot(tX, avepos, 'Color', 'red', 'DisplayName','synthetic', LineWidth=0.5);
    end

    xline(tX(training+1), '-', 'DisplayName', 'start of predictions', LineWidth=2);

    ylabel('popularity', 'FontSize', 16);
    %ylabel('expected degree', 'FontSize', 16);
    xlabel('time', 'FontSize', 16)

    legend('Location', 'Best');
    ax = gca; 
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    xlim([tX(1)-0.1 tX(end)+0.1]);
    set(gca, 'FontSize', 16); 
    box on;

    hold off;
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   out=strcat(mydir, Files(fn).name, ".png");
   saveas(gcf,out);

   count=count+1;

   %Print progress.
   fprintf("%d %s\n", count, Files(fn).name);

   %If you want to break earlier.
   % if count==100
   %    break;
   % end

end


