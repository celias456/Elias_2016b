function [CF] = cf(input,sT,EM,W,deepparameters,totalsimulations)
%Simulates model and outputs criterion function

%Input:
% input: parameters used in msm estimation
%   input(1) = mu;  Proportion of agents in economy who have a correctly specified PLM (agent 1)
%   input(2) = eshock_std;  Standard deviation of expectation shock
%   input(3) = eshock_ro;  AR coefficient on expectation shock
%   input(4) = gain; parameter for CG learning (agent 1)
%   input(5) = gain; parameter for CG learning (agent 2)
% sT: Number of time periods for each simulation
% EM: 1x4 matrix of empirical moments 
%   EM(1) = exchange rate return variance
%   EM(2) = exchange rate first-order autocovariance
%   EM(3) = exchange rate fourth-order autocovariance
%   EM(4) = deviation from fundamentals variance
%   EM(5) = deviation from fundamentals first-order autocovariance
% Weighting Matrix (W) (identity)
% deepparameters: deep parameters of model
%   parameters(1) = alpha; %Fundamentals equation constant
%   parameters(2) = delta; %Fundamentals equation parameter on time
%   parameters(3) = rho_f; %Fundamentals equation parameter on lagged fundamentals term
%   parameters(4) = sigma_epsilon; %Standard deviation of white noise term in Fundamentals equation

%Output:
% CF: criterion function value     
% -------------------------------------------------------------------------

%% DECLARATIONS
mu = input(1); %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = input(2); %Standard deviation of expectation shock
eshock_ro = input(3); %AR coefficient on expectation shock
gain1 = input(4); %Gain parameter for CG learning (agent 1)
gain2 = input(5); %Gain parameter for CG learning (agent 2)
gain = [gain1;gain2];

eshock_1_on = 1; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated
burn = sT/2; %Number of initial observations to remove

seedstate = 14; %Set state value for pseudo random number generator
stream = RandStream('mt19937ar','Seed',seedstate);
%RandStream.setDefaultStream(stream); %Not used in newer version of Matlab
RandStream.setGlobalStream(stream);
ahoc = .9; %Value for initial condition of learning parameters
sim_mom = zeros(totalsimulations,5); %Storage for simulated moment statistics

%% Deep/Country specific parameters
lambda = 8; %Interest rate semi-elasticity of Money Demand
a1 = lambda/(1+lambda);
b1 = 1/(1+lambda);

%Country Specific Parameters
alpha = deepparameters(1); %Fundamentals equation constant
delta = deepparameters(2); %Fundamentals equation parameter on time
rho_f = deepparameters(3); %Fundamentals equation parameter on lagged fundamentals term
sigma_epsilon = deepparameters(4); %Standard deviation of white noise term in Fundamentals equation
f_0 = 0.1; %Initial value for fundamentals
s_0 = 0.1; %Initial value for exchange rate

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

for count_simulations = 1:totalsimulations
    
    % Initialize variables
    wn = sigma_epsilon*randn(sT+2,1); %T White noise values for fundamentals
    eshock_wn = eshock_std*randn(sT+2,2); %Expectation shock white noise values

    % GET LAW OF MOTION FOR VARIABLES 
    [lom] = lom_learn_cg(burn,sT,f_0,s_0,a1,b1,alpha,delta,rho_f,wn,K,mu,eshock,eshock_wn,gain); 

    %GET Statistics
    sim_mom(count_simulations,:) = sim_moments(lom); %Calculate simulated moment statistics
    
end

%% Calculate Simulation Statistics and Criterion Function

statistics_moments = simstats_moments(sim_mom); %Simulated moments

%Calculate criterion function
g_sT = [statistics_moments.Returns_var-EM(1);statistics_moments.Returns_auto_1-EM(2);statistics_moments.Returns_auto_4-EM(3);statistics_moments.xi_var-EM(4);statistics_moments.xi_auto_1-EM(5)];
CF = (g_sT'*W*g_sT);

end