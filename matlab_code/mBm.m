function [mbm, ts, hs] = mBm(n, H, interval, mu, init, mysigma, fig)

%  Riemann-Liouville Multrifractional Brownian motion used in the paper
%  Fundamental dynamics of popularity-similarity trajectories in real networks, 
%  E. S. Papaefthymiou, C. Iordanou, F. Papadopoulos (2023), 
%  arXiv:2309.01675 

% This is a modification of the code from: 
% giannit (2023). Fractional and Multifractional Brownian motion generator 
% (https://github.com/Rabelaiss/mBm/releases/tag/1.2).

% Adapted from:
%    S. V. Muniandy and S. C. Lim (2001)
%    Modeling of locally self-similar processes using multifractional Brownian
%    motion of Riemann-Liouville type.
%    Physical Review E 63(4 Pt 2):046104
%    DOI: 10.1103/PhysRevE.63.046104


if nargin < 6
    error("At least 6 inputs required: mBm(n, H, interval, mu, init, mystd) with 'n' = length of the path, 'H' = Hurst function, 'interval' = vector with two components, 'mu' = average of increments (slope), 'init' = starting position of mBm, 'mystd'= standard deviation of increments")
elseif nargin == 6
    fig = 0;
end

if ~isscalar(n)
    error("'n' must be a number")
elseif n < 2
    error("'n' must be bigger than 1")
elseif n - floor(n) > 0
    n = round(n);
    warning("'n' must be an integer, it was rounded to the nearest integer (%d)", n)
end

if ~isa(H,'function_handle')    %If H is not a function.
    if (isscalar(H) || length(H) > 1)
        if ( min(H) <= 0 ) || (max(H) >= 1)
             error("'H' must be 0 < H < 1")
        else
            H = @(t) 0*t + H; 
        end
    else
        error("'H' must be a function or a number")
    end
end

if ~isa(mysigma,'function_handle')  %If mysigma is not a function.
    if ( min(mysigma) < 0 )
         error("'sigma' must be >= 0")
    else
        mysigma = @(t) 0*t + mysigma;
    end
end

if ~isnumeric(interval)
    error("'interval' must be a numeric vector")
elseif isempty(interval)
    interval = [0 1];
elseif (max(size(interval)) ~= 2) || (min(size(interval)) ~= 1)
    error("'interval' must be a vector with two components")
elseif diff(interval) <= 0
    error("'interval' must be a vector with two increasing components")
end

if ~isa(fig,'logical')
    if isscalar(fig)
        if (fig ~= 0) && (fig ~= 1)
            error("'fig' must be 0 or 1 (false or true)")
        end
    else
        error("'fig' must be 0, 1, false or true")
    end
end

ts = linspace(interval(1), interval(2), n)'; % time steps
hs = H(ts);                                  % Hurst steps;
sigma = mysigma(ts);                         %sigma steps;

if min(hs) <= 0
    warning("The Hurst function goes below 0 (min = %.2f) while it should be in the interval (0,1)", min(hs))
elseif max(hs) >= 1
    warning("The Hurst function goes above 1 (max = %.2f) while it should be in the interval (0,1)", max(hs))
end

t0 = 0;                                     % starting time for the mBm 
tn = 1;                                     % final time for the mBm 
dt = (tn - t0) / (n - 1);                   % step size
ts = linspace(t0, tn, n)';                  % time steps for the mBm

mbm = zeros(n,1);

%%%%%%%%%%Setup variance of gaussian noises%%%%%%%%%%%%%%
for i = 1 : length(hs)-1
    myvar(i) = (sigma(i)^2)*hs(i)*pi/(gamma(1-2*hs(i))*cos(hs(i)*pi));
    mystd(i) = sqrt(myvar(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Set the std of the noises.
rng('shuffle');                            % Set a random seed based on current time
xi  = mystd'.*randn(n-1,1);                % gaussian white noise = vector of random numbers from the standard normal distribution (i.e. mean=0 and std=1).
w   = @(t, H) sqrt((t .^ (2 * H) - (t - dt) .^ (2 * H)) ./ (2 * H * dt)) ./ gamma(H + 1/2); % eq 19 of the paper https://sci-hub.se/10.1103/PhysRevE.63.046104

%%%%%%%%%%Setup initial position of mBm%%%%%%%%%%
mbm(1)=init;

for k = 2:n                                          % skip k=1.
    weights  = w(ts(2:k), hs(k));
    mbm(k) = sqrt(dt)*sum(xi(1:k-1).*flip(weights)); % eq 17 of the paper 

    %%%%%%%%%%Give trend to the mBm and elevation (initial position)
    mbm(k)=mbm(k)+sum(mu(1:k-1))+init;
end

%%%%%%%%%%Output results in time 1 to n.
ts=1:n;
if fig
    figure
    %plot(ts,mbm,ts,hs)
    plot(ts,mbm,ts,sigma)
    ylim([min(mbm)-0.1 max(1,max(mbm)+0.1)])
end
