% Simulated Method of Moments Estimation

%Order of moments
%1. Variance of exchange rate return
%2. First-order autocovarinace of exchange rate returns
%3. Fourth-order autocovarinace of exchange rate returns
%3. Variance of deviation of fundamentals
%4. First-order autocovariance of deviation of fundamentals

%Order of parameters
%1. Degree of heterogeneity
%2. Standard deviation of expectation shock
%3. AR coefficient on expectation shock
%4. learning gain parameter (agent 1)
%5. learning gain parameter (agent 2)

%Functions used by this code:
% weightmatrix
% covnw
% cf
% numderiv
% g_st
% simulatemoments
% =========================================================================

clc;
clear all;
warning off

%% Set specifications for estimation
% United Kingdom = 'UK'; Germany = 'GER'; Japan = 'JPN'; Switzerland = 'SWS' 
ctry = 'SWS';
[parameters,T,burn,EM] = country(ctry);
s = 2; %Multiple of T periods to simulate (should be greater than one)
%tol = .001; %tolerance for numerical derivatives in SE calc

%Set specifications for simulation after estimation
totalsimulations = 2000;
timeperiods = burn + T;

%Calculation of sT (should be integer greater than one)
sT = s*T; %Number of time periods in each simulation

%Country Specific Parameters
alpha = parameters(1); %Fundamentals equation constant
delta = parameters(2); %Fundamentals equation parameter on time
rho_f = parameters(3); %Fundamentals equation parameter on lagged fundamentals term
sigma_epsilon = parameters(4); %Standard deviation of white noise term in Fundamentals equation
deepparameters = [alpha;delta;rho_f;sigma_epsilon];

%Calculate weighting matrix W
W = eye(5);

%% Calculate criterion function for theta_initial
mu_lb = .1;
mu_ub = .3;

eshock_std_lb = .1;
eshock_std_ub = .25;

eshock_ro_lb = .85;
eshock_ro_ub = .95;

gain1_lb = .15;
gain1_ub = .25;

gain2_lb = .001;
gain2_ub = .25;

theta_initial = [0.2,0.2,0.9,0.2,.01];
cfvalue_theta_initial = cf(theta_initial,sT,EM,W,deepparameters,totalsimulations);

%% Estimate parameters
%Constraints/Options on minimization
A = [];
b = [];
Aeq = [];
beq = [];
lb = [mu_lb,eshock_std_lb,eshock_ro_lb,gain1_lb,gain2_lb];
ub = [mu_ub,eshock_std_ub,eshock_ro_ub,gain1_ub,gain2_ub];
options = optimset('Algorithm','interior-point');
%Assign function handle
f = @(input)cf(input,sT,EM,W,deepparameters,totalsimulations);
%Perform minimization
[theta_hat,cfvalue_theta_hat] = fmincon(f,theta_initial,A,b,Aeq,beq,lb,ub,[],options);


%% Calculate Standard errors
% B = numderiv(theta_hat,tol,sT,EM,deepparameters);
% cov1 = ((B'*W*B)^-1);
% theta_hat_cov = ((1+(1/s))/T)*cov1;
% theta_hat_se = [sqrt(theta_hat_cov(1,1)),sqrt(theta_hat_cov(2,2)),sqrt(theta_hat_cov(3,3)),sqrt(theta_hat_cov(4,4))];


%% Simulate model with estimates
%statistics = simulatemoments(theta_hat,timeperiods,totalsimulations,burn,deepparameters);

