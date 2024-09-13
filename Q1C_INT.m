clear all 
close all 
clc
format short

% MU

a = -3; 
b = 3; 
N = 10^4;
xT = linspace(a, b, N);

%INT method
utr = 50; % truncation level of the integral

% Define characteristic function
cf = @(u, r) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, r, 1, 0.05);

ugrid = (0:N-1)*pi/(b-a);
CharFn = cf(ugrid, 0); % Initialize with r = 0

for j = 1:length(xT)
    utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)) - 10^(-15), utr); % Initial guess with r = 0
    utr = utrunc(j);
    % Numerical integration
    INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)), 0, utrunc(j))/pi, 0);
end

figure;
title('PDF estimated via INT Method for different values of \mu');
xlabel('Log Returns');

% Loop over different values of r for sensitivity analysis
r_values = [-1, -0.5, 0, 0.1, 0.5, 1];
num_r = length(r_values);

hold on;
for k = 1:num_r
    CharFn = cf(ugrid, r_values(k));
    INTpdf = zeros(size(xT));
    for j = 1:length(xT)
        utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, r_values(k))) - 10^(-15), utr);
        utr = utrunc(j);
        INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, r_values(k))), 0, utrunc(j))/pi, 0);
    end
    plot(xT, INTpdf, 'DisplayName', ['\mu = ' num2str(r_values(k))]);
end
hold off;

legend('show');

% THETA

clear all
clc

a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

%INT method
utr = 50; % truncation level of the integral

% Define characteristic function
cf = @(u, theta) cfHeston(u, 5, theta, 0.5, -0.7, 1, 0.1, 1, 0.05);

ugrid = (0:N-1)*pi/(b-a);
CharFn = cf(ugrid, 0); % Initialize with theta = 0

for j = 1:length(xT)
    utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)) - 10^(-15), utr); % Initial guess with theta = 0
    utr = utrunc(j);
    % Numerical integration
    INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)), 0, utrunc(j))/pi, 0);
end

figure;
title('PDF estimated via INT Method for different values of \theta');
xlabel('Log Returns');

% Loop over different values of r for sensitivity analysis
theta_values = [0.01, 0.03, 0.05, 0.08];
num_theta = length(theta_values);

hold on;
for k = 1:num_theta
    CharFn = cf(ugrid, theta_values(k));
    INTpdf = zeros(size(xT));
    for j = 1:length(xT)
        utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, theta_values(k))) - 10^(-15), utr);
        utr = utrunc(j);
        INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, theta_values(k))), 0, utrunc(j))/pi, 0);
    end
    plot(xT, INTpdf, 'DisplayName', ['\theta = ' num2str(theta_values(k))]);
end
hold off;

legend('show');

% ETA

clear all
clc

a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

%INT method
utr = 50; % truncation level of the integral

% Define characteristic function
cf = @(u, eta) cfHeston(u, 5, 0.05, eta, -0.7, 1, 0.1, 1, 0.05);

ugrid = (0:N-1)*pi/(b-a);
CharFn = cf(ugrid, 0.1); % Initialize with eta = 0.1

for j = 1:length(xT)
    utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0.1)) - 10^(-15), utr); % Initial guess with eta = 0.1
    utr = utrunc(j);
    % Numerical integration
    INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0.1)), 0, utrunc(j))/pi, 0);
end

figure;
title('PDF estimated via INT Method for different values of \eta');
xlabel('Log Returns');

% Loop over different values of r for sensitivity analysis
eta_values = [0.01, 0.05, 0.5, 1];
num_eta = length(eta_values);

hold on;
for k = 1:num_eta
    CharFn = cf(ugrid, eta_values(k));
    INTpdf = zeros(size(xT));
    for j = 1:length(xT)
        utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, eta_values(k))) - 10^(-15), utr);
        utr = utrunc(j);
        INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, eta_values(k))), 0, utrunc(j))/pi, 0);
    end
    plot(xT, INTpdf, 'DisplayName', ['\eta = ' num2str(eta_values(k))]);
end
hold off;

legend('show');

% RHO

clear all
clc

a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

%INT method
utr = 50; % truncation level of the integral

% Define characteristic function
cf = @(u, rho) cfHeston(u, 5, 0.05, 0.5, rho, 1, 0.1, 1, 0.05);

ugrid = (0:N-1)*pi/(b-a);
CharFn = cf(ugrid, 0); % Initialize with rho = 0

for j = 1:length(xT)
    utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)) - 10^(-15), utr); % Initial guess with rho = 0
    utr = utrunc(j);
    % Numerical integration
    INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)), 0, utrunc(j))/pi, 0);
end

figure;
title('PDF estimated via INT Method for different values of \rho');
xlabel('Log Returns');

% Loop over different values of r for sensitivity analysis
rho_values = [-0.9, -0.7, 0, 0.9];
num_rho = length(rho_values);

hold on;
for k = 1:num_rho
    CharFn = cf(ugrid, rho_values(k));
    INTpdf = zeros(size(xT));
    for j = 1:length(xT)
        utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, rho_values(k))) - 10^(-15), utr);
        utr = utrunc(j);
        INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, rho_values(k))), 0, utrunc(j))/pi, 0);
    end
    plot(xT, INTpdf, 'DisplayName', ['\rho = ' num2str(rho_values(k))]);
end
hold off;

legend('show');

% KAPPA

clear all
clc

a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

%INT method
utr = 50; % truncation level of the integral

% Define characteristic function
cf = @(u, kappa) cfHeston(u, kappa, 0.05, 0.5, -0.7, 1, 0.1, 1, 0.05);

ugrid = (0:N-1)*pi/(b-a);
CharFn = cf(ugrid, 0); % Initialize with kappa = 0

for j = 1:length(xT)
    utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)) - 10^(-15), utr); % Initial guess with kappa = 0
    utr = utrunc(j);
    % Numerical integration
    INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, 0)), 0, utrunc(j))/pi, 0);
end

figure;
title('PDF estimated via INT Method for different values of \kappa');
xlabel('Log Returns');

% Loop over different values of r for sensitivity analysis
kappa_values = [0.01, 0.7, 5, 6];
num_kappa = length(kappa_values);

hold on;
for k = 1:num_kappa
    CharFn = cf(ugrid, kappa_values(k));
    INTpdf = zeros(size(xT));
    for j = 1:length(xT)
        utrunc(j) = fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u, kappa_values(k))) - 10^(-15), utr);
        utr = utrunc(j);
        INTpdf(j) = max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u, kappa_values(k))), 0, utrunc(j))/pi, 0);
    end
    plot(xT, INTpdf, 'DisplayName', ['\kappa = ' num2str(kappa_values(k))]);
end
hold off;

legend('show');

function cf= cfHeston(u, kappa, theta, sigma, rho, tau, r, S0, v0)
% Heston parameters:
% kappa = variance mean reversion speed parameter
% theta = variance long−run level parameter
% rho = correlation between two Brownian motions
% sigma = volatility of variance
% v0 = initial variance
% S0 = initial stock price


% Log of the stock price.
x = log(S0);

% Parameter transformation
a = kappa*theta;
sg2=sigma^2;
d = sqrt((rho*sigma*1i*u - kappa).^2 - sg2*(1i*u - u.^2));
g = (kappa - rho*sigma*1i*u - d)./ (kappa - rho*sigma*1i*u + d);

% "Little Heston Trap" formulation
D1 = (kappa - rho*sigma*1i*u - d)/sg2;
D2 = ((1-exp(-d*tau))./(1-g.*exp(-d*tau)));
D = D1.*D2;
G = (1-g.*exp(-d*tau))./(1-g);
C = 1i*u*r*tau + a/sg2*((kappa - rho*sigma*1i*u- d)*tau - 2*log(G));
% The characteristic function.

cf = exp(C + D*v0 + 1i*u*x);
end