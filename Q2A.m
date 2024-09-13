clear all
close all 
clc
format short

rng(123)

nsimul=1000000; expiry=1/4; nsteps=63; % Trading days convention
dt=expiry/nsteps;
timestep=[0:dt:expiry]';
Wt=zeros(nsteps+1,nsimul);

%Simulate the Brownian motion at T:
eY = randn(1,nsimul);
Wt(nsteps+1,:)= sqrt(expiry).*eY;

%Simulate the Brownian motion W(t):
for j=2:nsteps
deltat1=(nsteps+1-j)/(nsteps+1-j+1);
eYt = randn(1,nsimul);
Wt(j,:)=deltat1*Wt(j-1,:)+(1-deltat1)*Wt(nsteps+1,:)+sqrt(deltat1*dt)*eYt;
end

%Parameters
S0=100;
mu=0.1;
sigma=0.2;
r=0.01;
K = 80:5:120;
t = [1 10];
p=0.05 ;   % quantile for VaR calculations

% Calculate 1 day and 10 day quantiles using the closed-form expression
Z_alpha = [norminv(p,0,1/nsteps)
           norminv(p,0,10/nsteps)];  

% Sorting Brownian bridge to get quantile values
Sorted_Wt = [sort(Wt(2,:),'ascend')
            sort(Wt(11,:),'ascend')]; 

St=zeros(nsteps,nsimul);

for h = 1:nsimul
    for i = 1:nsteps
        St(i,h) = S0*exp(mu*i/365 + sigma*Wt(i+1,h));
    end
end

% Plotting simulated Stock Prices  
% to see the plot change nsimul to 10,000 or else code takes too long to run
% h =figure('Color',[1 1 1])
% plot((1/nsteps:1/nsteps:1), St)
% title('Simulated Paths of the Stock Prices via Brownian Bridge')
% xlabel('Time (years)')

% Simulating 1 day and 10 day stock prices using risk neutral valuation
St1 = S0*exp(mu*t(1,1)/365 + sigma*Sorted_Wt(1,:));    % 1 day simulated stock prices via brownian bridge
St10 = S0*exp(mu*t(1,2)/365 + sigma*Sorted_Wt(2,:));   % 10 day simulated stock prices via brownian bridge 

% Calculating 1 day & 10 day log returns
rt1 = log(St1/S0);
rt10 = log(St10/S0);

% Simulating 1 day & 10 day stock prices and returns via Gaussian quantiles
St1_gauss = S0*exp(mu*t(1,1)/365 + sigma*Z_alpha(1,1)); 
St10_gauss = S0*exp(mu*t(1,2)/365 + sigma*Z_alpha(2,1)); 

rt1_gauss = log(St1_gauss/S0);
rt10_gauss = log(St10_gauss/S0);

% Calculating 5% quantiles for 1 & 10 day stock prices and returns generated via brownian bridge
Sh_1 =St1(1,p*nsimul);      
Sh_10 = St10(1,p*nsimul);   

rh_1 = rt1(1,p*nsimul);
rh_10 = rt10(1,p*nsimul);
disp('5% quantiles for 1 & 10 day holding period log returns generated via brownian bridge respectively:')
disp(rh_1)
disp(rh_10)


% Testing accuracy of holding period returns against gaussian returns
error_gauss = [(rt1_gauss-rh_1)
                (rt10_gauss-rh_10)]; 

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
error_table = table( error_gauss(1, :)', error_gauss(2, :)', 'VariableNames', {'1 Day', '10 Day'});

% Display tables
disp('VaR Table:');
disp(var_table);

disp('Error between Brownian Bridge & Gaussian');
disp(error_table);

function c = bs_call(s,k,r,tau,sigma)
d2=(log(s./(k.*exp(-r.*tau))))./(sigma.*sqrt(tau))-0.5*sigma.*sqrt(tau);
d1=d2+sigma.*sqrt(tau);
c=s.*normcdf(d1)-k.*exp(-r.*tau).*normcdf(d2);
end
