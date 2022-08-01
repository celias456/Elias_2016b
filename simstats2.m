function [stats] = simstats2(data)
%Calculates deviations from fundamentals statistics from a number of simulations

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


%Calculate median of statistics
s = median(data,1);

stats = struct('xi_std',s(1),'xi_vr_1',s(2),'xi_vr_8',s(3),'xi_vr_16',s(4),'xi_vr_32',s(5),'xi_acor_1',s(6),'xi_acor_4',s(7),'xi_acor_8',s(8),'xi_acor_16',s(9));



end

