function [data] = stats2(allvariables)
%Calculates Deviation from Fundamentals Statistics

%Input: 
%allvariables: 
%allvariables(:,1) = fundamentals
%allvariables(:,2) = exchange rate (logs)
%allvariables(:,3) = deviation from fundamentals


%Output:
% data: row vector of statistics
% 1. Deviation from fundamentals standard deviation
% 2.Deviation from fundamentals vrstat 1
% 3.Deviation from fundamentals vrstat 8
% 4.Deviation from fundamentals vrstat 16
% 5.Deviation from fundamentals vrstat 32
% 6. Deviation from fundamentals autocorrelation 1
% 7. Deviation from fundamentals autocorrelation 4
% 8. Deviation from fundamentals autocorrelation 8
% 9. Deviation from fundamentals autocorrelation 16

xi = allvariables(:,3); %Get deviations from fundamentals

%Calculate statistics for deviations from fundamentals
xi_std = std(xi);

%Variance Ratio Statistics
xi_vrstat_1 = vrstat(xi,1);
xi_vrstat_8 = vrstat(xi,8);
xi_vrstat_16 = vrstat(xi,16);
xi_vrstat_32 = vrstat(xi,32);

%Deviation from Fundamentls autocorrelation coefficients
xi_corr_1 = corr(xi(2:end),xi(1:end-1));
xi_corr_4 = corr(xi(5:end),xi(1:end-4));
xi_corr_8 = corr(xi(9:end),xi(1:end-8));
xi_corr_16 = corr(xi(17:end),xi(1:end-16));

data = [xi_std,xi_vrstat_1,xi_vrstat_8,xi_vrstat_16,xi_vrstat_32,xi_corr_1,xi_corr_4,xi_corr_8,xi_corr_16];


end

