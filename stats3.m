function [data] = stats3(allvariables)
%Calculates in-sample predictability statistics

%Input: 
%allvariables: 
%allvariables(:,1) = fundamentals
%allvariables(:,2) = exchange rate (logs)
%allvariables(:,3) = deviation from fundamentals


%Output:
% data: matrix of statistics
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

s = allvariables(:,2); %Get log exchange rates
xi = allvariables(:,3); %Get deviations from fundamentals

%Conduct Regressions 
[Beta_1,t_1,R_1,tstat_1] = nwestreg(s,xi,1);
[Beta_4,t_4,R_4,tstat_4] = nwestreg(s,xi,4);
[Beta_8,t_8,R_8,tstat_8] = nwestreg(s,xi,8);
[Beta_16,t_16,R_16,tstat_16] = nwestreg(s,xi,16);

data = [Beta_1,t_1,R_1,Beta_4,t_4,R_4,Beta_8,t_8,R_8,Beta_16,t_16,R_16,tstat_1,tstat_4,tstat_8,tstat_16];

end

