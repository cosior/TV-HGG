%Code for prediction of angular (similarity) trajectories.
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
nstd=1;

%BTC.
% th1=200;
% th2=1492;
% ksigma=8;
% alpha=0.1;
% iters=100;
% maxlag=10;
% nstd=1;

%PGP.
% th1=1;
% th2=1500;
% ksigma=8;
% alpha=0.1;
% iters=100;
% maxlag=10;
% nstd=1;

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
% alpha=0.1;
% iters=10;
% maxlag=20;
% nstd=1;

%Input dir (trajectories location, change accordingly);
Files=dir('./USAir/ANGULAR_TRAJECTORIES/*.txt');

%Output dir (where plots will be stored, change accordingly).
mydir="./USAir/ANGULAR_PREDS/";
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

   %Just in case we passed in an unwrapped trajectory.
   X = mod(X, 2*pi);
   
   clear velocities;
   velocities = zeros(n, 1);

    %Compute angular velocities
    for i = 2 : n
        dtheta = pi-abs(pi-abs(X(i)-X(i-1)));
        ux=cos(X(i-1));
        uy=sin(X(i-1));
        vx=cos(X(i));
        vy=sin(X(i));
        cross_p = ux*vy-uy*vx;
        velocities(i)=sign(cross_p)*dtheta;
    end

    %Test for velocities' trend-stationarity. 
    if (adftest(velocities)~=1)
        printf("Velocity series is not stationary for node %s\n", Files(fn).name);
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
    %simple heuristic: use the mean mu up to the training. 
    %other heuristics could work better depending on the trajectory/network.
    mu(training+1:end) = mean(mu(2:training));  

    %more general heuristic: use the mean mu over a previous window of length x 
    %(x could be optimized). Can work better, e.g., in BTC. Using median
    %instead of mean can also help, especially in cases where there are jumps, 
    %like in PGP and arXiv. Helps also in IPv6.
    % x = min(training-2, 200);
    % mu(training+1:end) = median(mu(training-x:training));
    %%mu(training+1:end) = mean(mu(training-x:training));

    %%%%%%%%%%%%%%%%%%%%%%%%%Hurst index%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear H Xun lags H_estimates;

    %Unwrapped trajectory.
    Xun=cumsum(velocities);

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
            %Compute H from the (unwrapped) trajectory up to the training.
            increments = abs(Xun(tau+1:training) - Xun(1:training-tau));  
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
    clear runun;

    t=1;

    while t <= iters

           clear Xsim; 

           Xsim = mBm(n, H, [ ], mu, X(1), sigma);
           runun(t,:) = Xsim;
           t = t+1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear aveposun variance stdposun;
   
    i=1:n;

    if (iters > 1)
        aveposun(i) = mean(runun(:,i));
        stdposun(i) = sqrt(var(runun(:,i)));
    else 
        aveposun(i) = runun(1,i);
        stdposun(i) = 0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Visible', 'off');

    hold on;
    
    if iters > 1

        clear curve1 curve2;

        %Consider an nstd region around the average. 
        curve1(i) = aveposun(i) + nstd*stdposun(i);
        curve2(i) = aveposun(i) - nstd*stdposun(i);

        clear i2;
        i2 = [tX(i), fliplr(tX(i))];
        inBetween = [curve1, fliplr(curve2)];
        fill(i2, inBetween, 'y', 'EdgeColor', 'none', 'DisplayName', 'std region');

    end
    
    %Unwrapped real trajectory.
    clear Xun;
    Xun=X(1)+cumsum(velocities);

    plot(tX, Xun, 'Color', 'blue', 'DisplayName','real', LineWidth=0.5);

    if iters > 1
        plot(tX, aveposun, 'Color', 'red', 'DisplayName','synthetic average', LineWidth=0.5);
    else
        plot(tX, aveposun, 'Color', 'red', 'DisplayName','synthetic', LineWidth=0.5);
    end

    xline(tX(training+1), '-', 'DisplayName', 'start of predictions', LineWidth=2);

    ylabel('similarity', 'FontSize', 16);
    xlabel('time', 'FontSize', 16)
 
    legend('Location', 'Best');
    ax = gca; 
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    xlim([tX(1)-0.1 tX(end)+0.1]);
    set(gca, 'FontSize', 16); 
    box on;

    hold off;
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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



