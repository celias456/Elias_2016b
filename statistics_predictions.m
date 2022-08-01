function [error] = statistics_predictions(allvariables)
%Calculates prediction statistics

price = allvariables(:,7); %previous period exchange rate
predict1 = allvariables(:,8); %Agent 1 previous period prediction
predict2 = allvariables(:,9); %Agent 2 previous period prediction
n = length(price);

error1 = (price - predict1).^2;
error1_sum = sum(error1);
mse_1 = error1_sum/n;

error2 = (price - predict2).^2;
error2_sum = sum(error2);
mse_2 = error2_sum/n;

error = [mse_1, mse_2];

end

