clear all 
close all 
clc
format short

% MU

% Define parameters
a = -3; 
b = 3; 
N = 10^4;
xT = linspace(a, b, N);

% Define characteristic function
cf = @(u, r) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, r, 1, 0.05);

% Create grid for numerical integration
ugrid = (0:N-1)*pi/(b-a);

% Loop over different values of r
r_values = [-1, -0.5, 0, 0.1,  0.5, 1];
num_r = length(r_values);

figure;
hold on;
for k = 1:num_r
    CharFn = cf(ugrid, r_values(k));
    COSpdf = zeros(size(xT));
    for j = 1:length(xT)
        V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
        COSpdf(j) = max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
    end
    plot(xT, COSpdf, 'DisplayName', ['\mu = ' num2str(r_values(k))]);
end
hold off;

xlabel('x');
ylabel('COS PDF');
legend('show');
title('COS PDF for different values of \mu');

clear all 
clc

% THETA

% Define parameters
a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

% Define characteristic function
cf = @(u, theta) cfHeston(u, 5, theta, 0.5, -0.7, 1, 0.1, 1, 0.05);

% Create grid for numerical integration
ugrid = (0:N-1)*pi/(b-a);

% Loop over different values of 
theta_values = [0.01, 0.03, 0.05, 0.08];
num_theta = length(theta_values);

figure;
hold on;
for k = 1:num_theta
    CharFn = cf(ugrid, theta_values(k));
    COSpdf = zeros(size(xT));
    for j = 1:length(xT)
        V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
        COSpdf(j) = max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
    end
    plot(xT, COSpdf, 'DisplayName', ['\theta = ' num2str(theta_values(k))]);
end
hold off;

xlabel('x');
ylabel('COS PDF');
legend('show');
title('COS PDF for different values of \theta');

% ETA

clear all 
clc

% Define parameters
a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

% Define characteristic function
cf = @(u, sigma) cfHeston(u, 5, 0.05, sigma, -0.7, 1, 0.1, 1, 0.05);

% Create grid for numerical integration
ugrid = (0:N-1)*pi/(b-a);

% Loop over different values of 
sigma_values = [0.01, 0.5, 1];
num_sigma = length(sigma_values);

figure;
hold on;
for k = 1:num_sigma
    CharFn = cf(ugrid, sigma_values(k));
    COSpdf = zeros(size(xT));
    for j = 1:length(xT)
        V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
        COSpdf(j) = max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
    end
    plot(xT, COSpdf, 'DisplayName', ['\sigma = ' num2str(sigma_values(k))]);
end
hold off;

xlabel('x');
ylabel('COS PDF');
legend('show');
title('COS PDF for different values of \sigma');

% RHO

clear all 
clc

% Define parameters
a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

% Define characteristic function
cf = @(u, rho) cfHeston(u, 5, 0.05, 0.5, rho, 1, 0.1, 1, 0.05);

% Create grid for numerical integration
ugrid = (0:N-1)*pi/(b-a);

% Loop over different values of 
rho_values = [-0.9, -0.7, 0, 0.9];
num_rho = length(rho_values);

figure;
hold on;
for k = 1:num_rho
    CharFn = cf(ugrid, rho_values(k));
    COSpdf = zeros(size(xT));
    for j = 1:length(xT)
        V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
        COSpdf(j) = max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
    end
    plot(xT, COSpdf, 'DisplayName', ['\rho = ' num2str(rho_values(k))]);
end
hold off;

xlabel('x');
ylabel('COS PDF');
legend('show');
title('COS PDF for different values of \rho');

% KAPPA

clear all 
clc

% Define parameters
a = -1; 
b = 1; 
N = 10^4;
xT = linspace(a, b, N);

% Define characteristic function
cf = @(u, kappa) cfHeston(u, kappa, 0.05, 0.5, -0.7, 1, 0.1, 1, 0.05);

% Create grid for numerical integration
ugrid = (0:N-1)*pi/(b-a);

% Loop over different values of 
kappa_values = [0.01, 0.7, 5, 6];
num_kappa = length(kappa_values);

figure;
hold on;
for k = 1:num_kappa
    CharFn = cf(ugrid, kappa_values(k));
    COSpdf = zeros(size(xT));
    for j = 1:length(xT)
        V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
        COSpdf(j) = max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
    end
    plot(xT, COSpdf, 'DisplayName', ['\kappa = ' num2str(kappa_values(k))]);
end
hold off;

xlabel('x');
ylabel('COS PDF');
legend('show');
title('COS PDF for different values of \kappa');

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