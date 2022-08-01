%% Monetary Model of Exchange Rates

% Fundamentals
% f_t = mu + delta*t + rho_f*f_tm1 + e_t
% e_t iid N(0,sigma^2_e)

% Misspecification
% Agent 2 includes the lag of the exchange rate in its PLM

% Variables
% 1. f_t: fundamentals (logs)
% 2. s_t: exchange rate (logs)
% 3. xi_t = f_t - s_t: deviation from fundamentals (logs)
% 4. exchange rate (levels)
% 5. fundamentals (levels)
% 6. Deviation from fundamentals (levels)

% Funtions used by this code:
% solve_ree
% lom_ree
% lom_learning_ls
% lom_learning_cg
% stats
%   returns
%   difference
% simstats

% ====================================================================================

clear all
clc;
warning off

%% Simulation Specifications
% United Kingdom = 'UK'; Germany = 'GER'; Japan = 'JPN'; Switzerland =
% 'SWS'
ctry = 'SWS';
% Law of motion type; 1 = RE; 2 = CG Learning
lom_type = 1;
N = 2000; % Number of simulations
ahoc = .9; %Value for initial condition of learning parameters
seedstate = 14; %Set state value for pseudo random number generator
stream = RandStream('mt19937ar','Seed',seedstate);
%RandStream.setDefaultStream(stream);
RandStream.setGlobalStream(stream);
s1 = zeros(N,9); %Storage for exchange rate returns statistics
s2 = zeros(N,9); %Storage for deviations from fundamentals statistics
s3 = zeros(N,16); %Storage for in sample predictability statistics
s4 = zeros(N,4); %Storage for out of sample forecasting statistics
sim_mom = zeros(N,6); %Storage for simulated moment statistics
data_predict_errors = zeros(N,2);

%% Learning/Expectation Shock Parameters
mu = 1; 
%mu = 0.1616; %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = 0.1886; %Standard deviation of the expectation shock
eshock_ro = 0.8502;
gain1 = 0.1504; %Gain parameter for CG learning - agent 1
gain2 = 0.0989; %Gain parameter for CG learning - agent 2
gain = [gain1;gain2];

eshock_1_on = 1; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated

%% Deep/Country specific parameters
lambda = 8; %Interest rate semi-elasticity of Money Demand
phi = 1; %Income elasticity of money demand
a1 = lambda/(1+lambda);
b1 = 1/(1+lambda);
[parameters,periods,burn,mom] = country(ctry);
% periods = 11000;
% burn = 1000;
alpha = parameters(1); %Fundamentals equation constant
delta = parameters(2); %Fundamentals equation parameter on time
rho_f = parameters(3); %Fundamentals equation parameter on lagged fundamentals term
sigma_epsilon = parameters(4); %Standard deviation of white noise term in Fundamentals equation
f_0 = .1; %Initial value for fundamentals
s_0 = .1; %Initial value for exchange rate
T = periods + burn;

%% Expectation Shock Parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock_std;  %Standard deviation of expectation shock for agent 2
eshock_1_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 1
eshock_2_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_1_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_2_ro];


%% Find REE solutions
[ree] = solve_ree(lambda,rho_f,delta,alpha);
 
%% Simulation
%Initialize variables for simulation
K = ahoc.*ree; %Starting values for parameter estimates;
 
for n = 1:N
     
    % Initialize variables
    wn = sigma_epsilon*randn(T+2,1); %T White noise values for fundamentals
    eshock_wn = eshock_std*randn(T+2,2); %Expectation shock white noise values

    % GET LAW OF MOTION FOR VARIABLES
    if lom_type == 1
    [lom] = lom_ree(burn,T,ree,a1,b1,f_0,alpha,delta,rho_f,wn,eshock,eshock_wn);
    else
    [lom,pf,count,Rtest1,Rtest2,Kappa1,Kappa2] = lom_learn_cg(burn,T,f_0,s_0,a1,b1,alpha,delta,rho_f,wn,K,mu,eshock,eshock_wn,gain); 
    end

    %GET STATISTICS
    s1(n,:) = stats1(lom); %Calculate exchange rate returns statistics
    s2(n,:) = stats2(lom); %Calculate deviations from fundamentals statistics
    s3(n,:) = stats3(lom); %Calculate in sample predictability statistics
    s4(n,:) = stats4(lom); %Calculate out of sample forecasting statistics
    sim_mom(n,:) = sim_moments(lom,mom); %Calculate simulated moment statistics
    
    if lom_type == 2
    [data_predict_errors(n,:)] = statistics_predictions(lom); %Calculate prediction errors
    end
end


%Calculate simulation statistics
statistics_1 = simstats1(s1); %Exchange rate returns
statistics_2 = simstats2(s2); %Deviations from fundamentals
statistics_3 = simstats3(s3,n); %In sample predictability 
statistics_4 = simstats4(s4); %Out of sample forecasting
statistics_moments = simstats_moments(sim_mom); %Simulated moments

if lom_type == 2
statistics_5 = simstats5(data_predict_errors); %Prediction errors
end


