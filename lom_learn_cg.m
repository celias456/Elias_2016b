function [lom,pf,count,Rtest1,Rtest2,Kappa1,Kappa2] = lom_learn_cg(burn,T,f_0,s_0,a1,b1,alpha,delta,rho_f,wn,K,mu,eshock,eshock_wn,gain) 
%% Simulates Law of motion for all variables under CG Adaptive Learning
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
% gain: gain in learning recursion
%
%Output: 
% lom (T-burn)x4 matrix of variables
%     allvariables(:,1) = fundamentals
%     allvariables(:,2) = exchange rate (logs)
%     allvariables(:,3) = deviation from fundamentals
%     allvariables(:,4) = exchange rate (levels)
%     allvariables(:,5) = fundamentals (levels)
%     allvariables(:,6) = deviation from fundamentals (levels)
%     allvariables(:,7) = previous period exchange rate
%     allvariables(:,8) = Agent 1 previous period prediction of er
%     allvariables(:,9) = Agent 2 previous period prediction of er
% -------------------------------------------------------------------------


%% Create storage matrices, declare variables, and initialize counters

lom = zeros(T+2,9); %Stores law of motion of all variables
Kappa1 = zeros(3,T+2); %Store values for the parameter updates (agent 1) 
Kappa2 = zeros(2,T+2); %Store values for the parameter updates (agent 2) 
R1 = zeros(3,3,T+2); %Store values for R matrix (agent 1)
R2 = zeros(2,2,T+2); %Store values for R matrix (agent 2)
expshocks = zeros(T+2,2); %Stores law of motion of expectation shocks

minmu = 1 - mu;

Rtest1 = 0; %Counts number of times the R1 matrix is not invertible
Rtest2 = 0; %Counts number of times the R2 matrix is not invertible
pf = 0; %Counts number of times projection facility is used

%% Initialize learning parameters

R1(:,:,2) = .01*ones(3,3); %Initial value for R matrix (agent 1)
R2(:,:,2) = .01*ones(2,2); %Initial value for R matrix (agent 2)

Kappa1(:,2) = K(1:3); %Initial values for parameters which are updated each period (agent 1)
Kappa2(:,2) = [K(1,1);.02]; %Initial values for parameters which are updated each period (agent 2)

Kappa_C1_0 = Kappa1(1,2); %Initial parameter value for constant (agent 1)
Kappa_T_0 = Kappa1(2,2); %Initial parameter value on time (agent 1)
Kappa_f_0 = Kappa1(3,2); %Initial parameter value on current fundamentals (agent 1)

Kappa_C2_0 = Kappa2(1,2); %Initial parameter value for constant (agent 2)
Kappa_0 = Kappa2(2,2); %Initial parameter value on lagged exchange rate (agent 2)

expshocks(2,1) = eshock(1,1)*eshock_wn(2,1); %Store initial value of expectation shock for agent 1
expshocks(2,2) = eshock(2,1)*eshock_wn(2,2); %Store initial value of expectation shock for agent 2

%% Initial value of fundamentals, expectations of future exchange rate,
% value of exchange rate, and deviation from fundamentals

f_initial = alpha + delta*1 + rho_f*f_0 + wn(2); %Initial value of fundamentals

%Agents form initial expectations
ES1_0 = Kappa_C1_0 + Kappa_T_0 + Kappa_f_0*(alpha+delta) + (Kappa_T_0 + Kappa_f_0*delta)*1 + Kappa_f_0*rho_f*f_initial + expshocks(2,1); %Agent 1
ES2_0 = Kappa_C2_0*(1+Kappa_0)  + (Kappa_0^2)*s_0 + expshocks(2,2);%Agent 2
ES_0 = mu*ES1_0 + (minmu)*ES2_0; %Average expectations in economy

%exchange rate is realized from reduced form equation
s_initial = a1*ES_0 + b1*f_initial; 

xi_initial = f_initial - s_initial; %deviation from fundamentals


% Store initial values of variables
lom(2,1) = f_initial; %Store initial value of fundamentals
lom(1,2) = s_0; %Store first period value of exchange rate (logs)
lom(2,2) = s_initial; %Store second period value of exchange rate (logs)
lom(2,3) = xi_initial; %Store initial value of deviation of fundamentals



%% Begin RLS learning recursion

for t = 3:(T+2)
       
   %% Set values for recursion
   tm2 = t-2;
   tm1 = t-1;
   f_tm1 = lom(tm1,1); %Fundamentals in previous period
   s_tm1 = lom(tm1,2); %Exchange rate one period previously
   s_tm2 = lom(tm2,2); %Exchange rate two periods previously
   
   Kappa1_tm1 = Kappa1(:,tm1); %Values of agent 1 coefficient estimates one period previously
   Kappa_C1_tm1 = Kappa1_tm1(1); %Agent 1 constant one period previously
   Kappa_T_tm1 = Kappa1_tm1(2); %Agent 1 coefficient on time one period previously
   Kappa_f_tm1 = Kappa1_tm1(3); %Agent 1 coefficient on fundamentals one period previously
   
   Kappa2_tm1 = Kappa2(:,tm1); %Values of agent 2 coefficient estimates one period previously
   Kappa_C2_tm1 = Kappa2_tm1(1); %Agent 2 constant one period previously
   Kappa_tm1 = Kappa2_tm1(2); %Agent 2 coefficient on lagged exchange rate one period previously
   
   R1_tm1 = R1(:,:,tm1); %Value of agent 1 R matrix one period previously
   R2_tm1 = R2(:,:,tm1); %Value of agent 2 R matrix one period previously

   shat1_tm1 = Kappa_C1_tm1 - Kappa_T_tm1 + Kappa_T_tm1*t + Kappa_f_tm1*f_tm1; %Agent 1 most recent predicted value for exchange rate
   prederror1 = s_tm1 - shat1_tm1; %Agent 1 most recent prediction error
   
   shat2_tm1 = Kappa_C2_tm1 + Kappa_tm1*s_tm2; %Agent 2 most recent predicted value for exchange rate
   prederror2 = s_tm1 - shat2_tm1; %Agent 2 most recent prediction error

   %% Agents update parameters by adding s_tm1 to the information set and
   %run a LS regression of s_tm1 on a constant, t, f_tm1, and s_tm2.

   %Update the R matrix - Agent 1
   R1_t = R1_tm1 + gain(1) * ( ([1;tm1;f_tm1] * [1,tm1,f_tm1]) - R1_tm1); 
   R1inv = (R1_t)^-1;
   %Update the coefficient matrix - Agent 1
   Kappa1_t = Kappa1_tm1 + gain(1) * R1inv * ([1;tm1;f_tm1]*(prederror1)); 

   %Update the R matrix - Agent 2
   R2_t = R2_tm1 + gain(2) * ( ([1;s_tm2] * [1,s_tm2]) - R2_tm1); 
   R2inv = (R2_t)^-1;
   %Update the coefficient matrix - Agent 1
   Kappa2_t = Kappa2_tm1 + gain(2) * R2inv * ([1;s_tm2]*(prederror2)); 
   
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
   
   %Projection Facility for agent 2
   if abs(Kappa_t) >= 1
      Kappa_C2_t = Kappa_C2_tm1; 
      Kappa_t = Kappa_tm1;
      pf = pf + 1;
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
   EP2_tp1 = Kappa_C2_t*(1+Kappa_t) + (Kappa_t^2)*s_tm1 + expshocks(t,2); %Agent 2
   EP_tp1 = mu*EP1_tp1 + (minmu)*EP2_tp1; %Average expectations in economy
   
   %s_T is realized from the reduced form of the model (or from ALM)
   s_t = a1*EP_tp1 + b1*f_t;
   
   %Get levels of variables
   S_t = exp(s_t); %exchange rate
   F_t = exp(f_t); %fundamentals
   DFF = F_t - S_t; %Deviation from fundamentals
   
   %Obtain law of motion for all variables
   lom(t,1) = f_t; %Fundamentals (logs)
   lom(t,2) = s_t; %Exchange rate (logs)
   lom(t,3) = f_t - s_t; %Deviation from fundamentals
   lom(t,4) = S_t; %Exchange rate (levels)
   lom(t,5) = F_t; %Fundamentals (levels)
   lom(t,6) = DFF; %Deviation from fundamentals (levels)
   
   %Obtain law of motion for agent forecast errors (levels)
   lom(t,7) = exp(s_tm1); %Previous period exchange rate
   lom(t,8) = exp(shat1_tm1); %Agent 1 previous period prediction of er
   lom(t,9) = exp(shat2_tm1); %Agent 2 previous period prediction of er
   
   
end
   
%% Drop the first two values in all the series
lom = lom(3:T+2,:);
lom = lom(burn+1:end,:);

Kappa1 = Kappa1(:,3:T+2);
Kappa1 = Kappa1(:,burn+1:end);
Kappa1 = Kappa1';

Kappa2 = Kappa2(:,3:T+2);
Kappa2 = Kappa2(:,burn+1:end);
Kappa2 = Kappa2';

if pf > 0
   count = 1;
else
   count = 0;
end

end