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

nsteps = 63; % trading days convention
nsimul = 1000000; % specify the number of simulations
expiry = 1/4;
dt = expiry/nsteps; % time step
timestep = (0:dt:expiry)';

% the increments of the Gamma clock
GY1T = gamrnd(dt/kappa,kappa,[nsteps,nsimul]);

% the increments of the Gaussian part
eY1 = randn(nsteps,nsimul);

% VG increments with target sum adjustment
dX = theta*GY1T + sigma*sqrt(GY1T).*eY1;

% Putting all the parts together
X = [zeros(1,nsimul); cumsum(dX)];

% Sorting VG bridge to get quantile values
Sorted_X = [sort(X(2,:),'ascend')
            sort(X(11,:),'ascend')]; 

% Simulating 1 day and 10 day stock prices using real probability measure
St1 = S0*exp(mu*t(1,1)/365 + Sorted_X(1,:));    % 1 day simulated stock prices via VG bridge
St10 = S0*exp(mu*t(1,2)/365 +Sorted_X(2,:));   % 10 day simulated stock prices via VG bridge

% Calculating 1 day & 10 day log returns
rt1 = log(St1/S0);
rt10 = log(St10/S0);

% Calculating 5% quantiles for 1 & 10 day stock prices and returns generated via VG bridge
Sh_1 =St1(1,p*nsimul);      
Sh_10 = St10(1,p*nsimul);   

rh_1 = rt1(1,p*nsimul);
rh_10 = rt10(1,p*nsimul);

disp('5% quantiles for 1 & 10 day holding period log returns generated via variance gamma bridge respectively:')
disp(rh_1)
disp(rh_10)

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
% title('Simulated Paths of the Stock Prices via VG Bridge')
% xlabel('Time (years)')


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
