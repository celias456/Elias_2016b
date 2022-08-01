function [stats] = simstats1(data)
%Calculates exchange rate returns statistics from a number of simulations


%Input:
% row vector of statistics
% 1.Exchange rate return standard deviation
% 2.Exchange rate vrstat 1
% 3.Exchange rate vrstat 8
% 4.Exchange rate vrstat 16
% 5 Exchange rate vrstat 32
% 6.Exchange rate autocorrelation coefficient 1
% 7.Exchange rate autocorrelation coefficient 8
% 8.Exchange rate autocorrelation coefficient 16
% 9. Exchange rate autocorrelation coefficient 32

%Calculate median of statistics
s = median(data,1);

stats = struct('Returns_std',s(1),'Returns_vr_1',s(2),'Returns_vr_8',s(3),'Returns_vr_16',s(4),'Returns_vr_32',s(5),'Returns_acor_1',s(6),'Returns_acor_4',s(7),'Returns_acor_8',s(8),'Returns_acor_16',s(9));
end

