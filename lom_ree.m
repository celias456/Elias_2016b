function [lom] = lom_ree(burn,T,ree,a1,b1,f_0,alpha,delta,rho_f,wn,eshock,eshock_wn)
%Simulates model under REE
%
%Input: 
% burn: Number of initial observations to remove
% T: Number of time periods
% ree: 3x1 matrix of REE solutions
%     ree(1) = constant
%     ree(2) = parameter on time
%     ree(3) = coefficient on fundamentals
% f_0: initial value for fundamentals
% alpha: Fundamentals equation constant
% delta: Fundamentals equation parameter on time
% rho_f: Fundamentals equation parameter on lagged fundamentals term
% wn: Tx1 matrix of white noise terms for fundamentals equation 
%
%Output: 
% lom: (T-burn)x4 matrix of variables
%     lom(:,1) = fundamentals
%     lom(:,2) = exchange rate (logs)
%     lom(:,3) = deviation from fundamentals

% -------------------------------------------------------------------------

%Initialize Storage
f = zeros(T+2,1); %Initialize storage for fundamentals
s = zeros(T+2,1); %Initialize storage for exchange rate 
xi = zeros(T+2,1); %Initialize storage for deviation from fundamentals;
S = zeros(T+2,1); %Initialize storage for exchange rate (levels);
F = zeros(T+2,1); %Initialize storage for fundamentals (levels);
DFF = zeros(T+2,1); %Initialize storage for deviation from fundamentals (levels);

expshocks = zeros(T+2,2); %Initialize storage of expectation shocks

a = ree(1);
b = ree(2);
c = ree(3);

%Set first period values
f_1 = alpha + delta*1 + rho_f*f_0 + wn(1); %First period value of fundamentals

expshocks(2,1) = eshock(1,1)*(eshock_wn(1,1)*eshock(1,4) + eshock_wn(2,1)); %Store initial value of expectation shock

Exp_1 = a + b + c*(alpha+delta) + (b+(c*delta))*1 + c*rho_f*f_1 + expshocks(2,1); %Agent expectations 

s_1 = a1*Exp_1 + b1*f_1; %Initial value of exchange rate (logs)

xi_1 = f_1 - s_1; %Initial value for deviation from fundamentals


%Store values
f(2,1) = f_1; %Fundamentals
s(2,1) = s_1; %Exchange Rate (logs)
xi(2,1) = xi_1; %Deviation From Fundamentals
S(2,1) = exp(s_1); %Exchange Rate (levels)
F(2,1) = exp(f_1); %Fundamentals (levels)
DFF(2,1) = F(2,1) - S(2,1); %Deviation from fundamentals (levels)


for t = 3:(T+2)
    tm1 = t-1;
    f_tm1 = f(tm1,1); %Previous period value of fundamentals
    f_t = alpha + delta*t + rho_f*f_tm1 + wn(t); %Current period value of fundamentals is realized
    
    expshocks(t,1) = eshock(1,1)*(expshocks(tm1,1)*eshock(1,4) + eshock_wn(t,1)); %Current value of expectation shock
    
    Exp_t = a + b + c*(alpha+delta) + (b+(c*delta))*t + c*rho_f*f_t + expshocks(t,1); %Agent expectations
       
    s_t = a1*Exp_t + b1*f_t; %Current period exchange rate is realized (logs)
    
    xi_t = f_t - s_t; %Current period deviation from fundamentals
    
    
    %Store law of motion for all variables
    f(t,1) = f_t; %Store current period value of fundamentals
    s(t,1) = s_t; %Store current period value of exchange rate in logs
    xi(t,1) = xi_t; %Store current period value the deviation from fundamentals
    S(t,1) = exp(s_t); %Store current period value of exchange rate in levels
    F(t,1) = exp(f_t); %Store current period value of fundamentals in levels
    DFF(t,1) = F(t,1) - S(t,1); %Store current period value of deviation from fundamentals in levels
    
end

%Remove initial (burn) observations
data = [f,s,xi,S,F,DFF];
data = data(3:T+2,:);
lom = data(burn+1:end,:);


end