function [stats_prediction_errors] = simstats5(data_predict_errors)
%Calculate statistics on prediction errors

%Calculate prediction error medians;
prediction_errors = median(data_predict_errors,1);

stats_prediction_errors = struct('Agent_1', prediction_errors(1), 'Agent_2', prediction_errors(2));


end

