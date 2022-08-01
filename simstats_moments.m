function [stats] = simstats_moments(data)
%Calculations simulated moment statistics from a number of simulations

%Input:
% data: row vector of data
% %Statistics are listed as:
% 1. Exchange rate return variance
% 2. Exchange rate return autocovariance 1
% 3. Exchange rate return autocovariance 4
% 4. Deviation from fundamentals variance
% 5. Deviation from fundamentals autocovariance 1
% 6. Criterion function

s = median(data,1);


stats = struct('Returns_var',s(1),'Returns_auto_1',s(2),'Returns_auto_4',s(3),'xi_var',s(4),'xi_auto_1',s(5),'CF',s(6));


end

