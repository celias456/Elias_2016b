function [Beta,t,R,tstat] = nwestreg(s,xi,k)
%Conduct Newey-West Regression

%Input:
    %s: exchange rate (column vector)
    %xi: deviation from fundamentals (column vector)
    %k: number of periods ahead 

%Output:
    %Beta = coefficient estimate
    %t = significance of t-stat (1=yes, 2=no)
    %R = Adjusted R-squared

%% Prepare Data

%Get length of original data series
T = length(s);

%Get length of new series
tmk = T-k;

%Get k-period ahead change in exchange rate
Y = kaheadchange(s,k);

%Obtain corresponding deviation from fundamentals data
X = xi(1:tmk,1);

%Add constant term to deviation from fundamentals data
X1 = [ones(tmk,1),X];    

    
%% Conduct Newey West regression
%Find lag-length for Newey-West Regression
l = floor(4*(tmk/100)^(2/9));

%Conduct Newey West regression
reg = nwest(Y,X1,l);


%% Get data from regression

%Get beta coefficient estimate
Beta = reg.beta(2);

%Get significance of t statistic
if reg.tstat(2) >= 1.96
    t = 1;
else
    t = 0;
end

%Get adjusted R-squared
R = reg.rbar;

if R < 0
   R = 0;
end

%Get t statistic
tstat = reg.tstat(2);

end

