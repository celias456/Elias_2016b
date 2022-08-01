function [vr] = vrstat(series,k)
%Calculates Variance Ratio Statistic
%k: lag length
%series: time series

%Find first-differenced series
fd = difference(series,1);

%Find variance of first-differenced series
var_fd = var(fd,1);

%Find k-differenced series
kd = difference(series,k);

%Find variance of k-differenced series
var_kd = var(kd,1);

%Calculate VR stat
vr = (1/k)*(var_kd/var_fd);

end

