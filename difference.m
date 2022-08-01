function [series] = difference(data,k)
%Computes the k-period difference of a series
%Input:
% data: time series
% k: lag

%Output:
% series: differenced series
% -------------------------------------------------------------------------


T = length(data); %Length of original series
series = zeros(T,1);
start = k+1;

for i = start:T
    series(i,1) = data(i) - data(i-k);
end

series = series(k+1:end,1);
end

