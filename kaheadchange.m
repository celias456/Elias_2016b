function [newdata] = kaheadchange(data,k)
%Gets the k-period ahead change of a series
%Input:
    %k: number of periods ahead
    %data: data series (column vector)

%%

%Obtain length of original series T
T = length(data);

%Obtain k-period ahead values from original data
k_period_ahead = data(1+k:end,1);

% Drop last k elemnts in original data
current_period = data(1:(T-k),1);

%Calculate k-period ahead change
newdata = k_period_ahead - current_period;

end

