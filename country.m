function [parameters,T,burn,mom] = country(string)
%Assigns parameters based on country

%Input
% string: country code United Kingdom = 'UK'; Germany = 'GER'; Japan = JAP

%Output
% parameters: 
%   parameters(1) = alpha; %Fundamentals equation constant
%   parameters(2) = delta; %Fundamentals equation parameter on time
%   parameters(3) = rho_f; %Fundamentals equation parameter on lagged fundamentals term
%   parameters(4) = sigma_epsilon; %Standard deviation of white noise term in Fundamentals equation
% T: number of time periods observed in each country's data set
% NW_data: Newey West data matrix
% mom: moments used in MSM estimation
%   mom(1) = exchange rate return standard deviation
%   mom(2) = exchange rate first-order autocorrelation
%   mom(3) = deviation of fundamentals standard deviation
%   mom(4) = deviation of fundamentals first-order autocorrelation
% -------------------------------------------------------------------------

if strcmpi(string,'UK') == 1
   alpha = 0.0981422; %Fundamentals equation constant
   delta = -0.0001038; %Fundamentals equation parameter on time
   rho_f = 0.9822844; %Fundamentals equation parameter on lagged fundamentals term
   sigma_epsilon = 0.02498; %Standard deviation of white noise term in Fundamentals equation
   %load NWUK.csv
   %NW_data = NWUK;
   T = 141; %Number of time periods in empirical data
   burn = T;
   mom = [22.5971;4.2855;1.1324;0.0464;0.0449]; %Empirical moments in data

elseif strcmpi(string,'GER') == 1
   alpha = 0.0194671; %Fundamentals equation constant
   delta = -0.0001544; %Fundamentals equation parameter on time
   rho_f = 0.9893052; %Fundamentals equation parameter on lagged fundamentals term
   sigma_epsilon = 0.02792; %Standard deviation of white noise term in Fundamentals equation   
   %load NWGER.csv
   %NW_data = NWGER;
   T = 112; %Number of time periods in empirical data
   burn = T;
   mom = [31.0221;6.1623;3.9782;0.1677;0.1630]; %Empirical moments in data

elseif strcmpi(string,'JPN') == 1
   alpha = -0.1718731; %Fundamentals equation constant
   delta = 0.0000818; %Fundamentals equation parameter on time
   rho_f = 0.9647813; %Fundamentals equation parameter on lagged fundamentals term
   sigma_epsilon = 0.02511; %Standard deviation of white noise term in Fundamentals equation  
   %load NWJAP.csv
   %NW_data = NWJAP;
   T = 176; %Number of time periods in empirical data
   burn = T;
   mom = [31.5086;4.0653;2.1094;0.2173;0.2127]; %Empirical moments in data

else
   alpha = 0.0528478; %Fundamentals equation constant
   delta = 0.0001275; %Fundamentals equation parameter on time
   rho_f = 0.9708941; %Fundamentals equation parameter on lagged fundamentals term
   sigma_epsilon = 0.0291; %Standard deviation of white noise term in Fundamentals equation  
   %load NWSWS.csv
   %NW_data = NWSWS;
   T = 173; %Number of time periods in empirical data
   burn = T;
   mom = [41.9369;0.8799;2.5127;0.0703;0.0666]; %Empirical moments in data 
end

parameters = [alpha;delta;rho_f;sigma_epsilon];

end

