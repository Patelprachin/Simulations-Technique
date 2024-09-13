clear all 
close all 
clc

% COS method
a=-1; b=1; N=10^4;
xT=(a:1/(10^4):b);

cf = @(u) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, 0.01, 1, 0.05);
ugrid=(0:N-1)*pi/(b-a);
CharFn=cf(ugrid);
for j=1:length(xT);
V = (2/(b-a))*cos((xT(j)-a)*(0:N-1)*pi/(b-a));
COSpdf(j)=max(real(sum(CharFn.*V.*exp(1i*(0:N-1)*pi*(-a)/(b-a)))-0.5*CharFn(1)*1*V(1)),0);
end

%INT method
utr=50; %truncation level of the integral
cf = @(u) cfHeston(u, 5, 0.05, 0.5, -0.7, 1, 0.01, 1, 0.05);
ugrid=(0:N-1)*pi/(b-a);
CharFn=cf(ugrid);

for j=1:length(xT)
utrunc(j)=fsolve(@(u) real(exp(-1i*u*xT(j)).*cf(u)) -10^(-15),utr);
utr=utrunc(j);
%numerical integration``````````````
qgkdf(j)=max(quadgk(@(u) real(exp(-1i*u*xT(j)).*cf(u)),0, utrunc(j))/pi,0);
end

%Moments using INT method
momINT=[trapz(xT,qgkdf) trapz(xT,xT.*qgkdf) trapz(xT,(xT.^2).*qgkdf) trapz(xT,(xT.^3).*qgkdf) trapz(xT,(xT.^4).*qgkdf)];

%Moments using COS method
momCOS=[trapz(xT,COSpdf) trapz(xT,xT.*COSpdf) trapz(xT,(xT.^2).*COSpdf) trapz(xT,(xT.^3).*COSpdf) trapz(xT,(xT.^4).*COSpdf)];

r=0.01;
tau=0.5;
df=exp(-r*tau);
K = (1);
european = [trapz(xT,max(exp(xT)-K,0).*qgkdf) trapz(xT,max(K-exp(xT),0).*qgkdf)
trapz(xT,max(exp(xT)-K,0).*COSpdf) trapz(xT,max(K-exp(xT),0).*COSpdf)]*exp(-r*tau)
%Compute BS implied vol
ivol = blsimpv(1,K,r,0.5,european(:,1))

% Parameters
S0 = 1;

% Define the function to calculate implied volatility using COS method
impliedVolatilityCOS = @(K, tau) blsimpv(S0, K, r, tau, european(:,1));

% Define strike prices and time to maturity grid
K_values = 0.5:0.1:1.5; % Strike prices
tau_values = 0:0.1:1.0; % Time to maturity

% Initialize matrices to store implied volatility values
impliedCallVolatilityMatrix = zeros(length(K_values), length(tau_values));
impliedPutVolatilityMatrix = zeros(length(K_values), length(tau_values));

% Loop through strike prices and time to maturity
for i = 1:length(K_values)
    for j = 1:length(tau_values)
        K = K_values(i);
        tau = tau_values(j);
        
        % Calculate implied volatility using COS method
        european = [trapz(xT,max(exp(xT)-K,0).*qgkdf) trapz(xT,max(K-exp(xT),0).*qgkdf)
                    trapz(xT,max(exp(xT)-K,0).*COSpdf) trapz(xT,max(K-exp(xT),0).*COSpdf)]*exp(-r*tau);
        impliedVolatility = impliedVolatilityCOS(K, tau);
        
        % Store implied volatility in the matrices
        impliedCallVolatilityMatrix(i, j) = impliedVolatility(1); % Call option implied volatility
        impliedPutVolatilityMatrix(i, j) = impliedVolatility(2); % Put option implied volatility
    end
end

% Plot the implied volatility surface for call options
figure;
surf(K_values,tau_values , impliedCallVolatilityMatrix);    
xlabel('Strike Price');                             
ylabel('Time to Maturity');                         
zlabel('Implied Volatility (Call Options)');
title('Implied Volatility Surface (Call Options)');

% Plot the implied volatility surface for put options
figure;
surf(tau_values, K_values, impliedPutVolatilityMatrix);
xlabel('Time to Maturity');
ylabel('Strike Price');
zlabel('Implied Volatility (Put Options)');
title('Implied Volatility Surface (Put Options)');

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




