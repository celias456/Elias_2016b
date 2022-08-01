%Graphs the Agent learning coefficients
%Order of variables in coefficients.csv
%1 (A). Period
%2 (B). Speculator beliefs learning constant
%3 (C). Speculator beliefs learning coefficient on time
%4 (D). Speculator beliefs learning coefficient on fundamentals
%5 (E). Non-Speculator beliefs learning constant
%6 (F). Non-Speculator beliefs learning lagged exchange rate
%7 (G). REE constant
%8 (H). REE coefficient on time
%9 (I). REE coefficient on fundamentals
%10 (K). REE coefficient on lagged exchange rate (all zeros)

load coefficients.csv
coefficient_data = coefficients;

period = coefficient_data(:,1);

speculator_constant_learning = coefficient_data(:,2);
constant_ree = coefficient_data(:,7);

speculator_time_learning = coefficient_data(:,3);
time_ree = coefficient_data(:,8);

speculator_fundamentals_learning = coefficient_data(:,4);
fundamentals_ree = coefficient_data(:,9);

nonspeculator_constant_learning = coefficient_data(:,5);
nonspeculator_exchangerate_learning = coefficient_data(:,6);
exchangerate_ree = coefficient_data(:,10);


figure
subplot(3,1,1)
plot(period,speculator_constant_learning,'k',period,constant_ree,'--k')
xlabel('Period')
ylabel('Constant')
legend('AL','RE');
%title('Constant')

subplot(3,1,2)
plot(period,speculator_time_learning,'k',period,time_ree,'--k')
xlabel('Period')
ylabel('Time')
legend('AL','RE');
axis([0,150,-0.06,0.08])
%title('Time')

subplot(3,1,3)
plot(period,speculator_fundamentals_learning,'k',period,fundamentals_ree,'--k')
xlabel('Period')
ylabel('Fundamentals')
legend('AL','RE');
axis([0,150,-6,6])
%title('Fundamentals')

% mtit('Forecasting Model Paramater Values: Speculator');


figure
subplot(2,1,1)
plot(period,nonspeculator_constant_learning,'k',period,constant_ree,'--k')
xlabel('Period')
ylabel('Constant')
legend('AL','RE');
axis([0,150,-0.1,1])
%title('Constant')

subplot(2,1,2)
plot(period,nonspeculator_exchangerate_learning,'k',period,exchangerate_ree,'--k')
xlabel('Period')
ylabel('Lagged Exchange Rate')
legend('AL','RE');
%title('Lagged Exchange Rate')
axis([0,150,-0.1,1.4])

% mtit('Forecasting Model Paramater Values: Non-Speculator');
