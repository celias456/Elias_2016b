function [statistics] = simulatemoments(input,T,totalsimulations,burn,deepparameters)
%Simulates model and outputs simulated momeents

% Fundamentals
% f_t = mu + delta*t + rho_f*f_tm1 + e_t
% e_t iid N(0,sigma^2_e)

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
% T: Number of time periods for each simulation
% totalsimulations: number of replications
% burn: number of initial observiations to remove in each replication
% deepparameters: deep parameters of model
%   parameters(1) = alpha; %Fundamentals equation constant
%   parameters(2) = delta; %Fundamentals equation parameter on time
%   parameters(3) = rho_f; %Fundamentals equation parameter on lagged fundamentals term
%   parameters(4) = sigma_epsilon; %Standard deviation of white noise term in Fundamentals equation

%Output:
% 1.Exchange rate return mean
% 2.Exchange rate return standard deviation
% 3.Exchange rate autocorrelation 1
% 4.Exchange rate autocorrelation 4
% 5.Exchange rate autocorrelation 8
% 6 Exchange rate autocorrelation 16
% 7 Deviation from fundamentals mean
% 8 Deviation from fundamentals standard deviation
% 9 Deviation from fundamentals autocorrelation 1
% 10 Deviation from fundamentals autocorrelation 4
% 11 Deviation from fundamentals autocorrelation 8
% 12 Deviation from fundamentals autocorrelation 16
% 13 Exchange rate return variance
% 14 Exchange rate return autocovariance 1
% 15 Exchange rate return autocovariance 4
% 16 Deviation from fundamentals variance
% 17 Deviation from fundamentals autocovariance 1
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

seedstate = 14; %Set state value for pseudo random number generator
stream = RandStream('mt19937ar','Seed',seedstate);
RandStream.setDefaultStream(stream);
ahoc = .9; %Value for initial condition of learning parameters
data = zeros(totalsimulations,17);

%% DEEP PARAMETERS
lambda = 8; %Interest rate semi-elasticity of Money Demand
a1 = lambda/(1+lambda);
b1 = 1/(1+lambda);

%Country Specific Parameters
alpha = deepparameters(1); %Fundamentals equation constant
delta = deepparameters(2); %Fundamentals equation parameter on time
rho_f = deepparameters(3); %Fundamentals equation parameter on lagged fundamentals term
sigma_epsilon = deepparameters(4); %Standard deviation of white noise term in Fundamentals equation
f_0 = 0; %Initial value for fundamentals
s_0 = 0; %Initial value for exchange rate

%% Expectation Shock Parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock_std;  %Standard deviation of expectation shock for agent 2
eshock_1_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 1
eshock_2_ro = eshock_ro; %First-order autocorrelation term for expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_1_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_2_ro];

%% FIND REE SOLUTIONS
[ree] = solve_ree(lambda,rho_f,delta,alpha);

%% Simulation

%INITIALIZE VARIABLES FOR SIMULATION
K = ahoc.*ree; %Starting values for parameter estimates;

for count_simulations = 1:totalsimulations
    
    % Initialize variables
    wn = sigma_epsilon*randn(T+2,1); %T White noise values for fundamentals
    eshock_wn = eshock_std*randn(T+2,2); %Expectation shock white noise values

    % GET LAW OF MOTION FOR VARIABLES 
    [lom] = lom_learn_cg(burn,T,f_0,s_0,a1,b1,alpha,delta,rho_f,wn,K,mu,eshock,eshock_wn,gain); 

    %GET STATISTICS
    [s] = stats(lom);
    data(count_simulations,:) = s;
    
end

%% Calculate Simulation Statistics an Criterion Function

statistics = simstats(data);


