clear all
close all
clc

rng(786)

% Parameters
theta = -0.01;
kappa = 0.045;
sigma = 0.2;
S0 = 100;
mu = 0.11;
r = 0.01;
K = 80:5:120;
t = [1 10];
p=0.05 ; 

% Calculating Xt value such that Stock Price is 80 at end of 1 year
syms x
Xt_365 = double(vpasolve(80 == S0*exp(mu + x),x));

nsteps = 365; % Using actual day count convention for stress testing purposes
nsimul = 1000000; % specify the number of simulations
expiry = 1/4;
dt = expiry/nsteps; % time step
timestep = (0:dt:expiry)';

% the increments of the Gamma clock
GY1T = gamrnd(dt/kappa,kappa,[nsteps,nsimul]);

% the increments of the Gaussian part
eY1 = randn(nsteps,nsimul);

% VG increments with target sum adjustment
target_sum = Xt_365;
dX = theta*GY1T + sigma*sqrt(GY1T).*eY1;
dX = dX - mean(dX, 1) + target_sum/nsteps;

X = [zeros(1,nsimul); cumsum(dX)];

% Sorting VG bridge to get quantile values
Sorted_X = [sort(X(2,:),'ascend')
            sort(X(11,:),'ascend')]; 

% Simulating 1 day and 10 day stock prices using real probability measure 
St1 = S0*exp(mu*t(1,1)/365 + Sorted_X(1,:));    % 1 day simulated stock prices via VG bridge
St10 = S0*exp(mu*t(1,2)/365 + Sorted_X(2,:));   % 10 day simulated stock prices via VG bridge 

% Checking Stock Price convergence to target price (80)
St=zeros(nsteps,nsimul);

for h = 1:nsimul
    for i = 1:nsteps
        St(i,h) = S0*exp(mu*i/365 + X(i+1,h));
    end
end

% Plotting simulated Stock Prices  
% to see the plot change nsimul to 10,000 or else code takes too long to run
% h =figure('Color',[1 1 1])
% plot((1/nsteps:1/nsteps:1), St)
% title('Simulated Paths of the Stock Prices via Variance-Gamma Bridge')
% xlabel('Time (years)')

Sh_1 =St1(1,p*nsimul);      % 1 day 5% quantile of simulated stock prices
Sh_10 = St10(1,p*nsimul);    % 10 day 5% quantile of simulated stock prices

%Initializing Variables
C0=zeros(1,length(K));
C1=zeros(1,length(K));
C10=zeros(1,length(K));

% Calculating call option prices over a range of strike prices
for h = 1:length(K);
     K1 = K(h);

    for i =1;
        C0(i,h) = bs_call(S0(i),K1,r*92/365,92/365,sigma);  % Actual day count convention as per black scholes formula
        C1(i,h)= bs_call(Sh_1(i),K1,r*91/365,91/365,sigma);
        C10(i,h)= bs_call(Sh_10(i),K1,r*81/365,81/365,sigma);

    end
end

% Calculating VaR as per formulation : C(h) - C(0)
VaR_1= C1-C0;
VaR_10= C10-C0;

% Create tables
var_table = table(K', VaR_1', VaR_10', 'VariableNames', {'Strike Price', '1 day 95% VaR', '10 day 95% VaR'});

% Display tables
disp('VaR Table:');
disp(var_table);

function c = bs_call(s,k,r,tau,sigma)
d2=(log(s./(k.*exp(-r.*tau))))./(sigma.*sqrt(tau))-0.5*sigma.*sqrt(tau);
d1=d2+sigma.*sqrt(tau);
c=s.*normcdf(d1)-k.*exp(-r.*tau).*normcdf(d2);
end
