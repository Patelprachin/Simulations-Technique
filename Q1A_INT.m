clear all 
close all 
clc
format short

a=-1; b=1; N=10^4;
xT=(a:1/(10^4):b);

%INT method
utr=50; %truncation level of the integral
cf=@(u) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, 0.1, 1, 0.05);
ugrid=(0:N-1)*pi/(b-a);
CharFn=cf(ugrid);

for j=1:length(xT)
utrunc(j)=fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u)) -10^(-15),utr);
utr=utrunc(j);
%numerical integration
INTpdf(j)=max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u)),0, utrunc(j))/pi,0);
end

plot(xT,INTpdf,'-g')
title('PDF estimated via INT Method')
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




