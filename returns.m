function [r,count] = returns(series)
%Calculate percentage returns of a series in logs
%Input
% series: time series

%Output
% r: vector of returns (1 element less than series)
% count: number of times denominator is equal to zero
% -------------------------------------------------------------------------

count = 0;
T = length(series);
r = zeros(T,1); %Store exchange rate returns

for t = 2:T
    tm1 = t-1;
    r(t) = (series(t)/series(tm1)-1)*100;
    
    if series(tm1) == 0
       series(tm1) = 1e-323;
       count = count + 1;
    end
end

%Drop the first entry (it is equal to zero)
r = r(2:end);



