clear all 
close all 
clc
format short

% COS method
a=-1; b=1; N=10^4;
xT=(a:1/(10^4):b);

cf=@(u) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, 0.1, 1, 0.05);
ugrid=(0:N-1)*pi/(b-a);
CharFn=cf(ugrid);
for j=1:length(xT);
V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
COSpdf(j)=max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
end

plot(xT,COSpdf,'r')
title('PDF estimated via COS Method')
xlabel('Log Returns')

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



