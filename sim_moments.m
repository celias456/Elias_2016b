function [data] = sim_moments(allvariables,EM)
%Calculates moments of simulated data

%Input: 
%allvariables: 
%allvariables(:,1) = fundamentals
%allvariables(:,2) = exchange rate (logs)
%allvariables(:,3) = deviation from fundamentals
%allvariables(:,4) = exchange rate (levels)

%Output:
% data: row vector of statistics
% 1. Exchange rate return variance
% 2. Exchange rate return autocovariance 1
% 3. Exchange rate return autocovariance 4
% 4. Deviation from fundamentals variance
% 5. Deviation from fundamentals autocovariance 1

%% Exchange rate returns
s = allvariables(:,2); %Get log exchange rates
s_returns = returns_log(s); %Calculate exchange rate returns (percentage)

%Calculate exchange rate return variance and autocovariance
s_returns_var = var(s_returns,1);
s_auto1_matrix = cov(s_returns(2:end),s_returns(1:end-1),1);
s_auto1 = s_auto1_matrix(2,1);
s_auto4_matrix = cov(s_returns(5:end),s_returns(1:end-4));
s_auto4 = s_auto4_matrix(2,1);

%% Deviation from fundamentals
xi = allvariables(:,3); %Get deviation from fundamentals

%Calculate deviation from fundamentals variance and autocovariance
xi_var = var(xi,1);
xi_auto1_matrix = cov(xi(2:end),xi(1:end-1));
xi_auto1 = xi_auto1_matrix(2,1);

%Calculate criterion function
g_sT = [s_returns_var-EM(1);s_auto1-EM(2);s_auto4-EM(3);xi_var-EM(4);xi_auto1-EM(5)];
W = eye(5);
CF = (g_sT'*W*g_sT);

data = [s_returns_var,s_auto1,s_auto4,xi_var,xi_auto1,CF];

end

