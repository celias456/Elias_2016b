function [data] = stats1(allvariables)
%Calculates Exchange Rate Statistics

%Input: 
%allvariables: 
%allvariables(:,1) = fundamentals
%allvariables(:,2) = exchange rate (logs)
%allvariables(:,3) = deviation from fundamentals


%Output:
% row vector of statistics
% 1.Exchange rate return standard deviation
% 2.Exchange rate vrstat 1
% 3.Exchange rate vrstat 8
% 4.Exchange rate vrstat 16
% 5.Exchange rate vrstat 32
% 6.Exchange rate autocorrelation coefficient 1
% 7.Exchange rate autocorrelation coefficient 8
% 8.Exchange rate autocorrelation coefficient 16
% 9. Exchange rate autocorrelation coefficient 32

s = allvariables(:,2); %Get logs exchange rates
[s_returns] = returns_log(s); %Calculate exchange rate returns (percentage)


%Calculate statistics for exchange rate returns
s_returns_std = std(s_returns); %returns standard deviation
test = isinf(s_returns_std);

%Variance Ratio Statistics
s_vrstat_1 = vrstat(s,1);
s_vrstat_8 = vrstat(s,8);
s_vrstat_16 = vrstat(s,16);
s_vrstat_32 = vrstat(s,32);

%Returns Autocorrelation Coefficients
s_corr_1 = corr(s_returns(2:end),s_returns(1:end-1));
s_corr_4 = corr(s_returns(5:end),s_returns(1:end-4));
s_corr_8 = corr(s_returns(9:end),s_returns(1:end-8));
s_corr_16 = corr(s_returns(17:end),s_returns(1:end-16));

%Test to see if exchange rate return standard deviation is infinity
if test == 1
   s_returns_std = 100000000000000000000000000000000000000000;
   s_corr_1 = 0; 
   s_corr_4 = 0;
   s_corr_8 = 0;
   s_corr_16 = 0;
end

%% Return the statistics data
% %Statistics are listed as:
% 1.Exchange rate return standard deviation
% 2.Exchange rate vrstat 1
% 3.Exchange rate vrstat 8
% 4.Exchange rate vrstat 16
% 5. Exchange rate vrstat 32
% 6.Exchange rate autocorrelation coefficient 1
% 7.Exchange rate autocorrelation coefficient 8
% 8.Exchange rate autocorrelation coefficient 16
% 9. Exchange rate autocorrelation coefficient 32

data = [s_returns_std,s_vrstat_1,s_vrstat_8,s_vrstat_16,s_vrstat_32,s_corr_1,s_corr_4,s_corr_8,s_corr_16];

end

