function [g_sT] = g_st(input,sT,EM,deepparameters)
%Simulates model and outputs criterion function

% Misspecification
% Agent 2 omits the fundamentals in its PLM

% Variables   1. fundamentals 
%             2. exchange rate (logs)         
%             3. deviation from fundamentals   
%             4. exchange rate (levels) 
%Input:
% input: parameters used in msm estimation
%   input(1) = mu;  Proportion of agents in economy who have a correctly specified PLM (agent 1)
%   input(2) = eshock_std;  Standard deviation of expectation shock
%   input(3) = eshock_ro;  AR coefficient on expectation shock
%   input(4) = gain; parameter for CG learning
% sT: Number of time periods for each simulation
% EM: 1x4 matrix of empirical moments 
%   EM(1) = exchange rate return standard deviation
%   EM(2) = exchange rate first-order autocorrelation
%   EM(3) = deviation from fundamentals standard deviation
%   EM(4) = deviation from fundamentals first-order autocorrelation
% Weighting Matrix (W)
% deepparameters: deep parameters of model
%   parameters(1) = alpha; %Fundamentals equation constant
%   parameters(2) = delta; %Fundamentals equation parameter on time
%   parameters(3) = rho_f; %Fundamentals equation parameter on lagged fundamentals term
%   parameters(4) = sigma_epsilon; %Standard deviation of white noise term in Fundamentals equation

%Output:
% g_sT: 4x1 matrix of deviations from empirical moments
%   g_sT(1) = deviation from exchange rate return standard deviation
%   g_sT(2) = deiviation from exchange rate return first-order AC
%   g_sT(3) = deviation from deviation from fundamentals standard deviation
%   g_sT(4) = deviation from deviation from fundamentals first-order AC
% -------------------------------------------------------------------------

%% DECLARATIONS
mu = input(1); %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = input(2); %Standard deviation of expectation shock
eshock_ro = input(3); %AR coefficient on expectation shock
gain = input(4); %Gain parameter for CG learning

totalsimulations = 1; %Number of simulations to run
eshock_1_on = 1; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated
burn = sT/2; %Number of initial observations to remove

seed_state = 98; %Set state value for pseudo random number generator
randn('state',seed_state);
ahoc = .9; %Value for initial condition of learning parameters
data = zeros(totalsimulations,20);

%% DEEP PARAMETERS
lambda = 8; %Interest rate semi-elasticity of Money Demand

%Country Specific Parameters
alpha = deepparameters(1); %Fundamentals equation constant
delta = deepparameters(2); %Fundamentals equation parameter on time
rho_f = deepparameters(3); %Fundamentals equation parameter on lagged fundamentals term
sigma_epsilon = deepparameters(4); %Standard deviation of white noise term in Fundamentals equation
f_0 = 0; %Initial value for fundamentals

%% Expectation Shock Parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock_std;  %Standard deviation of expectation shock for agent 2
eshock_1_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 1
eshock_2_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_1_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_2_ro];

%% FIND REE SOLUTIONS
[ree] = SOLVE(lambda,rho_f,delta,alpha);

%% Simulation

%INITIALIZE VARIABLES FOR SIMULATION
K = ahoc.*ree; %Starting values for parameter estimates;

for count_simulations = 1:totalsimulations
    
    % Initialize variables
    wn = sigma_epsilon*randn(sT,1); %T White noise values for fundamentals
    eshock_wn = eshock_std*randn(sT,2); %Expectation shock white noise values

    % GET LAW OF MOTION FOR VARIABLES 
    [lom] = LEARN_CG(burn,sT,f_0,alpha,delta,rho_f,lambda,wn,K,mu,eshock,eshock_wn,gain); 

    %GET STATISTICS
    [s] = STATS(lom);
    data(count_simulations,:) = s;
    
end

%% Calculate Simulation Statistics an Criterion Function

statistics = SIMSTATS(data);

%Calculate criterion function
g_sT = [statistics.Returns_std-EM(1);statistics.Returns_rho1-EM(2);statistics.xi_std-EM(3);statistics.xi_rho1-EM(4)];

end