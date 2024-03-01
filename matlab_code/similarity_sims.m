%Code for creating simulated counterparts of angular (similarity) trajectories.

%Requires the subplot_tight function from here: 
%https://www.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot

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
chunks=100;
alpha=0.01;
iters=10;
maxlag=10;
nstd=1;

%BTC.
% th1=200;
% th2=1492;
% ksigma=8;
% chunks=10;
% alpha=0.01;
% iters=1;
% maxlag=10;
% nstd=1;

%PGP.
% th1=1;
% th2=1500;
% ksigma=8;
% chunks=10;
% alpha=0.05;
% iters=1;
% maxlag=10;
% nstd=1;

%IPv6.
% th1=150;
% th2=450;
% ksigma=8;
% chunks=10;
% alpha=0.1;
% iters=1;
% maxlag=10;
% nstd=1;

%arXiv.
% th1=2800;
% th2=6776;
% ksigma=8;
% chunks=10;
% alpha=0.01;
% iters=1;
% maxlag=20;
% nstd=1;
% flag=1;

%Input dir (trajectories' location, change accordingly).
Files=dir('./USAir/ANGULAR_TRAJECTORIES/*.txt');

%Output dir (where plots will be stored, change accordingly).
mydir1="./USAir/ANGULAR_SIMS/";
mkdir(mydir1);
%Output dir for storing hurst exponents (change accordingly).
mydir2="./USAir/";
outh=strcat(mydir2, "sim_hursts.txt");
fileIDh = fopen(outh, 'w');

count=0;

for fn = 1:length(Files)

   %%If you want to check a specific node.
   % if Files(fn).name ~= "Node_CLT.txt"
   %       continue;
   % end

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
   
   %Just in case we passed in an unwrapped trajectory.
   X = mod(X, 2*pi);
   
   clear velocities;
   velocities = zeros(n, 1);

   %Angular velocities. 
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

    %Exponentially weighted moving average.
    for i = 2:n
        mu(i) = alpha*velocities(i)+(1-alpha)*mu(i - 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%Hurst exponent%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            increments = abs(Xun(tau+1:end) - Xun(1:end-tau));  
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
    S = zeros(n,1);
    sigma = zeros(n,1);

    for j = 1 : n

         a = max(j-ksigma/2, 1);
         b = min(j+ksigma/2, n);

         clear window_data;

         window_data=velocities(a:b);
    
         %Using 1st moment.  
         c = sqrt(2)/gamma(0.5);
         S(j) = mean(abs(window_data));
         sigma(j) = S(j)/c;
         sigma(j) = sigma(j)/n^(-H);
         
    end

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
    clear aveposun stdposun;
   
    i=1:n;

    if (iters > 1)
        aveposun(i) = mean(runun(:,i));
        stdposun(i) = sqrt(var(runun(:,i)));
    else 
        aveposun(i) = runun(1,i);
        stdposun(i) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%trajectory plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iters > 1

        clear curve1 curve2;

        %Consider an nstd region around the average. 
        curve1(i) = aveposun(i) + nstd*stdposun(i);
        curve2(i) = aveposun(i) - nstd*stdposun(i);

        figure('Visible', 'off');
        vmargin=0.22;
        hmargin=0.05;
        subplot_tight(1, 6, 1, [vmargin hmargin]);
        hold on;

        clear i2;
        i2 = [tX(i), fliplr(tX(i))];
        inBetween = [curve1, fliplr(curve2)];
        fill(i2, inBetween, 'y', 'EdgeColor', 'none', 'DisplayName', 'std region');

    else

        figure('Visible', 'off');
        vmargin=0.22;
        hmargin=0.05;
        subplot_tight(1, 6, 1, [vmargin hmargin]);
        hold on;

    end

    %Unwrapped real trajectory.
    clear Xun;
    Xun=X(1)+cumsum(velocities);

    plot(tX, Xun, 'Color', 'blue', 'DisplayName','real', LineWidth=0.5);

    if iters > 1
        plot(tX, aveposun, 'Color', 'red', 'DisplayName','synth. ave.', LineWidth=0.5);
    else
        plot(tX, aveposun, 'Color', 'red', 'DisplayName','synthetic', LineWidth=0.5);
    end

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

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%velocities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   subplot_tight(1, 6, 2, [vmargin hmargin]);
   hold on;
  
   clear Vsynth tmp
    
   Vsynth=zeros(n, 1);
   %Plot a single run here.
   tmp = runun(1,:);
   Vsynth(2:end)=diff(tmp);

   plot(tX, Vsynth, 'Color', 'red', 'DisplayName','synthetic', LineWidth=0.5);
   plot(tX, velocities, 'Color', 'blue', 'DisplayName','real', LineWidth=0.5);

   ylabel('velocity', 'FontSize', 16);
   xlabel('time', 'FontSize', 16);

   ax = gca; 
   ax.XAxis.MinorTick = 'on';
   ax.YAxis.MinorTick = 'on';
   xlim([tX(1)-0.1 tX(end)+0.1]);
   set(gca, 'FontSize', 16); 
   box on;

   hold off;

   %%%%%%%%%%%%%%%%%%%%%%%%%%PDF of velocities%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear kde evaluated_kdes average_kde;

   subplot_tight(1, 6, 5, [vmargin hmargin]);
   hold on;

   mygrid = linspace(-pi, pi, 1000);

   kde = ksdensity(velocities, mygrid);

   %For arXiv the bandwidth better set to 0.01.
   if flag==1
       clear kde;
       kde = ksdensity(velocities, mygrid, 'Bandwidth', 0.01);
   end

   plot(mygrid, kde,  'Color', 'blue', 'DisplayName', 'real', LineWidth=2);

   %std of real velocities.
   std_v_real = sqrt(var(velocities));
    
   std_v=0;
   std_t=0;

   for t = 1:iters
    
         clear Vsynth tmp kde;
    
         Vsynth=zeros(n, 1);
         tmp = runun(t,:);
         Vsynth(2:end)=diff(tmp);

         std_v = std_v+sqrt(var(Vsynth));   

         clear kde;
         kde=ksdensity(Vsynth, mygrid);

         %For arXiv.
         if flag==1
             clear kde;
             kde = ksdensity(Vsynth, mygrid, 'Bandwidth', 0.01);
         end

         evaluated_kdes(:, t) = kde;
    end

    %Average std of synthetic velocities.
    std_v = std_v/iters;

    %Compute the average density.
    clear average_kde;
    average_kde = mean(evaluated_kdes, 2);
    plot(mygrid, average_kde, 'Color', 'red', 'DisplayName', 'synth.', 'LineStyle', '-.', LineWidth=2 );

    %Juxtapose a normal pdf with the same mean and std.
    clear mu sigma;
    [mu,sigma] = normfit(velocities);

    pdf = normpdf(mygrid, mu, sigma);
    plot(mygrid, pdf, 'Color', 'black', 'DisplayName', 'Gaus.', 'LineStyle', '--', LineWidth=2);

    xlim([min(velocities) max(velocities)]);
   
    ylabel('pdf', 'FontSize', 16);
    xlabel('velocity', 'FontSize', 16);
    legend('Location', 'Best');

    ax = gca; 
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    set(gca, 'FontSize', 16); 
    box on;

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%variance-time plot%%%%%%%%%%%%%%%%%%%%%%%%%
    clear myvar myvar2;

    subplot_tight(1, 6, 4, [vmargin hmargin]);

    m = 1;
 	while (m <= ceil(n/chunks))
        myvar(m) = var(sum(reshape(velocities(1:end-mod(length(velocities), m)), m, []))/m);
        myvar2(m)=0;
   	    m = m + 1;
    end
   
    for t = 1:iters
    
         clear Vsynth tmp

         Vsynth=zeros(n, 1);
         tmp = runun(t,:);
         Vsynth(2:end)=diff(tmp);
        
         m = 1;
 	     while (m <= ceil(n/chunks))
               myvar2(m) = myvar2(m)+var(sum(reshape(Vsynth(1:end-mod(length(Vsynth), m)), m, []))/m);
   	           m = m + 1;
         end
    end

    myvar2 = myvar2/iters;
     
    m = 1:ceil(n/chunks);

    clear binCenters averages;
    [binCenters, averages] = computeLogBinnedAverages(m, myvar(m), 15);

    clear validIndices xValid yValid;
    validIndices = ~isnan(binCenters) & ~isnan(averages);
    xValid = binCenters(validIndices);
    yValid = averages(validIndices);

    loglog(xValid, yValid, 'k--', 'DisplayName', 'real', 'MarkerFaceColor', 'blue', LineWidth=1, MarkerSize=20, Marker='square');
    hold on;
    
    clear binCenters averages;
    [binCenters, averages] = computeLogBinnedAverages(m, myvar2(m), 15);

    clear validIndices xValid yValid;
    validIndices = ~isnan(binCenters) & ~isnan(averages);
    xValid = binCenters(validIndices);
    yValid = averages(validIndices);

    loglog(xValid, yValid, 'k--', 'DisplayName', 'synthetic', 'MarkerFaceColor', 'red', LineWidth=1, MarkerSize=15, Marker='o');

    loglog(m, var(velocities)*m.^(2*(mean(H)-1)), '-', 'Color', 'k', 'DisplayName', '\sigma^2 m^{2H-2}', LineWidth=2);
   
    ylabel('variance', 'FontSize', 16);
    xlabel('block size (m)', 'FontSize', 16);
    legend('Location', 'Best');

    xlim([2 max(m)+1]);
    ax = gca; 
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    set(gca, 'FontSize', 16); 
    box on;

    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%Autocorrelation plot%%%%%%%%%%%%%%%%%%%%%%%%%
    clear acfr lagsr acfs lagss acfsave

    [acfr,lagsr, bounds] = autocorr(velocities);

    subplot_tight(1, 6, 3, [vmargin hmargin]);
    hold on;

    %Consider lags 1-5 (lag 0 is just 1)
    plot(lagsr(2:6), acfr(2:6), 'k--', 'DisplayName', 'real', 'MarkerFaceColor', 'blue', LineWidth=1, MarkerSize=20, Marker='square' );
    
    acfsave = zeros(max(lagsr)+1, 1);

    for t = 1:iters
    
         clear Vsynth;

         Vsynth=zeros(n, 1);
         tmp = runun(t,:);
         Vsynth(2:end)=diff(tmp);
    
         clear acfs lagss;
         [acfs,lagss] = autocorr(Vsynth);
         acfsave = acfsave+acfs;
    end

    acfsave = acfsave/iters;
    plot(lagsr(2:6), acfsave(2:6), 'k--', 'DisplayName', 'synthetic', 'MarkerFaceColor', 'red', LineWidth=1, MarkerSize=15, Marker='o' );
  
    %Fractional Gaussian noise.
    autocorr_FGN = 0.5*((lagsr-1).^(2*mean(H)) + (lagsr+1).^(2*mean(H)) - 2*lagsr.^(2*mean(H)));
    plot(lagsr(2:6), autocorr_FGN(2:6),'DisplayName', 'fGn', 'Color', 'black', 'LineStyle', '-', LineWidth=2) ;  

    ylabel('autocorrelation','FontSize', 16);
    xlabel('lag', 'FontSize', 16);
    legend('Location', 'Best');

    xlim([1 5]);
    ax = gca; 
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    set(gca, 'FontSize', 16); 
    box on;

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%Polar histogram%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot_tight(1, 6, 6, [0.035 0.035]);
    
    clear Hist;
    Hist=polarhistogram(X, 'BinWidth', 0.1, 'Normalization', 'probability', 'FaceColor', 'blue');   

    hold on;
    
    %If iters > 1 we plot the distribution of the average trajectory.
    Hist=polarhistogram(aveposun, 'BinWidth', 0.1, 'Normalization', 'probability', 'FaceColor', 'red');   

    title('polar histogram','FontWeight','normal');
    set(gca, 'FontSize', 16); 

    hold off;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %Save the plots.
   desiredWidth = 1800;
   desiredHeight = 350;
   set(gcf, 'Position', [100, 100, desiredWidth, desiredHeight]);
   out=strcat(mydir1, Files(fn).name, ".png");
   saveas(gcf,out);

   %Output some stats.
   fprintf("%d %s H=%.4f sigma_real=%.2f sigma_sim=%.2f mu_real=%f mu_sim=%f\n", fn, Files(fn).name, mean(H), std_v_real, std_v, mean(velocities), mean(diff(aveposun)));
   
   %Save the Hurst exponent.
   fprintf(fileIDh,"%s %f\n", Files(fn).name, H);

   count=count+1;

   %%If you want to break earlier.
   % if count==100
   %    break;
   % end

end

fclose(fileIDh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For log-binning the variance-time plot.
function [binCenters, averages] = computeLogBinnedAverages(xData, yData, numBins)
   
    % Compute the bin edges based on xData
    minX = min(xData);
    maxX = max(xData);
    binEdges = logspace(log10(minX), log10(maxX), numBins + 1);

    % Initialize arrays to store bin-wise sum and count
    binSum = zeros(1, numBins);
    binCount = zeros(1, numBins);

    % Iterate over the xData and yData
    for i = 1:numel(xData)

         % Find the bin index for the current xData point
         binIndex = find(xData(i) >= binEdges, 1, 'last');

          % Check if binIndex is within valid range
         if ~isempty(binIndex) && binIndex <= numBins

            % Accumulate the sum and count for the corresponding bin
            binSum(binIndex) = binSum(binIndex) + yData(i);
            binCount(binIndex) = binCount(binIndex) + 1;
        
         end
       
    end

    % Compute bin centers
    binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));

    % Compute averages per bin
    averages = binSum ./ binCount;

end