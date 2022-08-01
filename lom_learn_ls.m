function [lom,Rtest1,Rtest2] = lom_learn_ls(burn,T,f_0,alpha,delta,rho_f,lambda,wn,K,mu,eshock,eshock_wn) 
%% Simulates Law of motion for all variables under RLS Adaptive Learning
%
%Input: 
% burn: Number of initial observations to remove
% T: Number of time periods
% f_0: starting value of fundamentals
% alpha: Fundamentals equation constant
% delta: Fundamentals equation parameter on time
% rho_f: Fundamentals equation parameter on lagged fundamentals term
% wn: Tx1 matrix of white noise terms for fundamentals equation 
% K: 3x1 matrix of initial values for parameter estimates
% mu: Proportion of agents with correctly specified PLM
% eshock: 2x4 matrix of expectation shock specifications
% eshock_wn: white noise terms for expectation shock
%
%Output: 
% lom (T-burn)x4 matrix of variables
%     allvariables(:,1) = fundamentals
%     allvariables(:,2) = exchange rate (logs)
%     allvariables(:,3) = deviation from fundamentals

% -------------------------------------------------------------------------


%% Create storage matrices, declare variables, and initialize counters

lom = zeros(T,3); %Stores law of motion of all variables
Kappa1 = zeros(3,T); %Store values for the parameter updates (agent 1) 
Kappa2 = zeros(2,T); %Store values for the parameter updates (agent 2) 
R1 = zeros(3,3,T); %Store values for R matrix (agent 1)
R2 = zeros(2,2,T); %Store values for R matrix (agent 2)
expshocks = zeros(T,2); %Stores law of motion of expectation shocks

minmu = 1 - mu;
a = lambda/(1+lambda); %Parameter on average expectations
b = 1/(1+lambda); %Parameter on fundamentals

Rtest1 = 0; %Counts number of times the R1 matrix is not invertible
Rtest2 = 0; %Counts number of times the R2 matrix is not invertible

%% Initialize learning parameters

R1(:,:,1) = .01*ones(3,3); %Initial value for R matrix (agent 1)
R2(:,:,1) = .01*ones(2,2); %Initial value for R matrix (agent 2)

Kappa1(:,1) = K(1:3); %Initial values for parameters which are updated each period (agent 1)
Kappa2(:,1) = K(1:2); %Initial values for parameters which are updated each period (agent 2)

Kappa_C1_0 = Kappa1(1,1); %Initial parameter value for constant (agent 1)
Kappa_T_0 = Kappa1(2,1); %Initial parameter value on time (agent 1)
Kappa_f_0 = Kappa1(3,1); %Initial parameter value on current fundamentals (agent 1)

Kappa_C2_0 = Kappa2(1,1); %Initial parameter value for constant (agent 2)
Kappa_0 = Kappa2(2,1); %Initial parameter value on time shock (agent 2)

expshocks(1,1) = eshock(1,1)*eshock_wn(1,1); %Store initial value of expectation shock for agent 1
expshocks(1,2) = eshock(2,1)*eshock_wn(1,2); %Store initial value of expectation shock for agent 2

%% Initial value of fundamentals, expectations of future exchange rate,
% value of exchange rate, and deviation from fundamentals

f_initial = alpha + delta*1 + rho_f*f_0 + wn(1); %Initial value of fundamentals

%Agents form initial expectations
ES1_0 = Kappa_C1_0 + Kappa_T_0 + Kappa_f_0*(alpha+delta) + (Kappa_T_0 + Kappa_f_0*delta)*1 + Kappa_f_0*rho_f*f_initial + expshocks(1,1); %Agent 1
ES2_0 = Kappa_C2_0 + Kappa_0  + Kappa_0*1 + expshocks(1,2);%Agent 2
ES_0 = mu*ES1_0 + (minmu)*ES2_0; %Average expectations in economy

%exchange rate is realized from reduced form equation
s_initial = a*ES_0 + b*f_initial; 

xi_initial = f_initial - s_initial; %deviation from fundamentals


% Store initial values of variables
lom(1,1) = f_initial; %Store initial value of fundamentals
lom(1,2) = s_initial; %Store initial value of exchange rate (logs)
lom(1,3) = xi_initial; %Store initial value of deviation of fundamentals



%% Begin RLS learning recursion

for t = 2:T
       
   %% Set values for recursion
   tm1 = t-1;
   f_tm1 = lom(tm1,1); %Fundamentals in previous period
   s_tm1 = lom(tm1,2); %Exchange rate one period previously
   
   gain = 1/t; %Gain 
   
   Kappa1_tm1 = Kappa1(:,tm1); %Values of agent 1 coefficient estimates one period previously
   Kappa_C1_tm1 = Kappa1_tm1(1); %Agent 1 constant one period previously
   Kappa_T_tm1 = Kappa1_tm1(2); %Agent 1 coefficient on time one period previously
   Kappa_f_tm1 = Kappa1_tm1(3); %Agent 1 coefficient on fundamentals one period previously
   
   Kappa2_tm1 = Kappa2(:,tm1); %Values of agent 2 coefficient estimates one period previously
   Kappa_C2_tm1 = Kappa2_tm1(1); %Agent 2 constant one period previously
   Kappa_tm1 = Kappa2_tm1(2); %Agent 2 coefficient on time one period previously
   
   R1_tm1 = R1(:,:,tm1); %Value of agent 1 R matrix one period previously
   R2_tm1 = R2(:,:,tm1); %Value of agent 2 R matrix one period previously

   shat1_tm1 = Kappa_C1_tm1 - Kappa_T_tm1 + Kappa_T_tm1*t + Kappa_f_tm1*f_tm1; %Agent 1 most recent predicted value for equity price
   prederror1 = s_tm1 - shat1_tm1; %Agent 1 most recent prediction error
   
   shat2_tm1 = Kappa_C2_tm1 - Kappa_tm1 + Kappa_tm1*t; %Agent 2 most recent predicted value for equity price
   prederror2 = s_tm1 - shat2_tm1; %Agent 2 most recent prediction error

   %% Agents update parameters by adding s_tm1 to the information set and
   %run a LS regression of s_tm1 on a constant, t, and f_tm1.

   %Update the R matrix - Agent 1
   R1_t = R1_tm1 + gain * ( ([1;t;f_tm1] * [1,t,f_tm1]) - R1_tm1); 
   R1inv = (R1_t)^-1;
   %Update the coefficient matrix - Agent 1
   Kappa1_t = Kappa1_tm1 + gain * R1inv * ([1;t;f_tm1]*(prederror1)); 

   %Update the R matrix - Agent 2
   R2_t = R2_tm1 + gain * ( ([1;t] * [1,t]) - R2_tm1); 
   R2inv = (R2_t)^-1;
   %Update the coefficient matrix - Agent 1
   Kappa2_t = Kappa2_tm1 + gain * R2inv * ([1;t]*(prederror2)); 
   
   Kappa_C1_t = Kappa1_t(1,1);
   Kappa_T_t = Kappa1_t(2,1);
   Kappa_f_t = Kappa1_t(3,1);
   
   Kappa_C2_t = Kappa2_t(1,1);
   Kappa_t = Kappa2_t(2,1);

   %This is a check to make sure that R1 and R2 are invertible
   if R1inv(1,1) == inf || R1inv(1,2) == inf || R1inv(1,3) == inf || R1inv(2,1) == inf || R1inv(2,2) == inf || R1inv(2,3) == inf || R1inv(3,1) == inf || R1inv(3,2) == inf || R1inv(3,3) == inf
      Kappa_C1_t = Kappa_C1_tm1; 
      Kappa_T_t = Kappa_T_tm1;
      Kappa_f_t = Kappa_f_tm1;
      Rtest1 = Rtest1 + 1;
   end
   if R2inv(1,1) == inf || R2inv(1,2) == inf || R2inv(2,1) == inf || R2inv(2,2) == inf
      Kappa_C2_t = Kappa_C2_tm1; 
      Kappa_t = Kappa_tm1;
      Rtest2 = Rtest2 + 1;
   end

   
   %Update Kappa1, R1, Kappa2, and R2 with values of new parameters
   Kappa1(1,t) = Kappa_C1_t;
   Kappa1(2,t) = Kappa_T_t;
   Kappa1(3,t) = Kappa_f_t;
   R1(:,:,t) = R1_t;
   
   Kappa2(1,t) = Kappa_C2_t;
   Kappa2(2,t) = Kappa_t;
   R2(:,:,t) = R2_t;
   
   %% The value of the current fundamentals (f_t) is realized
   f_t = alpha + delta*t + f_tm1*rho_f + wn(t);
   
   %Expectation shocks are realized
   expshocks(t,1) = eshock(1,1)*(expshocks(tm1,1)*eshock(1,4) + eshock_wn(t,1)); %Current value of expectation shock for agent 1
   expshocks(t,2) = eshock(2,1)*(expshocks(tm1,2)*eshock(2,4) + eshock_wn(t,2)); %Current value of expectation shock for agent 2
   
   %Agents form expectations of the future price of exchange rate 
   EP1_tp1 = Kappa_C1_t + Kappa_T_t + Kappa_f_t*(alpha+delta) + (Kappa_T_t+Kappa_f_t*delta)*t + Kappa_f_t*rho_f*f_t + expshocks(t,1); %Agent 1
   EP2_tp1 = Kappa_C2_t + Kappa_t + Kappa_t*t + expshocks(t,2); %Agent 2
   EP_tp1 = mu*EP1_tp1 + (minmu)*EP2_tp1; %Average expectations in economy
   
   %s_T is realized from the reduced form of the model (or from ALM)
   s_t = a*EP_tp1 + b*f_t;
   
   
   %Obtain law of motion for all variables
   lom(t,1) = f_t; %Price
   lom(t,2) = s_t; %Exchange rate (logs)
   lom(t,3) = f_t - s_t; %Deviation from fundamentals
   
   
   
end
   
%% Drop the first value in all the series
lom = lom(burn+1:end,:);

end