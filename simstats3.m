function [stats] = simstats3(data,n)
%Calculates in sample productivity statistics from a number of simulations

%Input:
% 1. In sample predictability regression coefficient (k=1)
% 2. In sample predictability significance (k=1)
% 3. In sample predictability R-squared (k=1)
% 4. In sample predictability regression coefficient (k=4)
% 5. In sample predictability significance (k=4)
% 6. In sample predictability R-squared (k=4)
% 7. In sample predictability regression coefficient (k=8)
% 8. In sample predictability significance (k=8)
% 9. In sample predictability R-squared (k=8)
% 10. In sample predictability regression coefficient (k=16)
% 11. In sample predictability significance (k=16)
% 12. In sample predictability R-squared (k=16)
% 13. In sample predictability t-stat (k=1)
% 14. In sample predictability t-stat (k=4)
% 15. In sample predictability t-stat (k=8)
% 16. In sample predictability t-stat (k=16)


%% In Sample Predictability
%Calculate median regression coefficients
b1 = median(data(:,1));
b4 = median(data(:,4));
b8 = median(data(:,7));
b16 = median(data(:,10));

%Calculate percentage of t-stats that are significant
sig1 = (sum(data(:,2))/n);
sig4 = (sum(data(:,5))/n);
sig8 = (sum(data(:,8))/n);
sig16 = (sum(data(:,11))/n);

%Calculate median R-squared
r1 = median(data(:,3));
r4 = median(data(:,6));
r8 = median(data(:,9));
r16 = median(data(:,12));

%Calculate median t-statistics
tstat_1 = median(data(:,13));
tstat_4 = median(data(:,14));
tstat_8 = median(data(:,15));
tstat_16 = median(data(:,16));

stats = struct('beta_1',b1,'t_1',sig1,'R_1',r1,'beta_4',b4,'t_4',sig4,'R_4',r4,'beta_8',b8,'t_8',sig8,'R_8',r8,'beta_16',b16,'t_16',sig16,'R_16',r16, 'tstat_1', tstat_1,'tstat_4', tstat_4,'tstat_8', tstat_8,'tstat_16', tstat_16);
end

