function [stats] = simstats4(data)
%Calculates out of sample forecasting statistics for a number of
%simulations

%Input:
% data: row vector of data
% 1. Out of sample predictability U-stat (1 period ahead)
% 2. Out of sample predictability U-stat (4 period ahead)
% 3. Out of sample predictability U-stat (8 period ahead)
% 4. Out of sample predictability U-stat (16 period ahead)

s = median(data,1);

stats = struct('U_1',s(1),'U_4',s(2),'U_8',s(3),'U_16',s(4));

end

