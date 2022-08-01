function [r] = returns_log(series)
%Calculate percentage returns of a series in natural logs
%Input
% series: time series

%Output
% r: vector of returns (1 element less than series)
% -------------------------------------------------------------------------

%Get first-difference series (log growth rate)
diff = difference(series,1);

%Find percentages of log growth rate series
r = (exp(diff)-1)*100;


end

