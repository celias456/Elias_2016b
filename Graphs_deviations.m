%Graphs the trajectory of the deviation from fundamentals
%Order of variables in deviations.csv
%1 (A). Period
%2 (B). Deviation from fundamentals for RE1 Model
%3 (C). Deviation from fundamentals for SN Model

load deviations.csv
deviations_data = deviations;

period = deviations_data(:,1);
re1 = deviations_data(:,2);
sn = deviations_data(:,3);

figure
plot(period,sn,'k',period,re1,'--k')
xlabel('Period')
ylabel('Deviation from Fundamentals');
legend('SN','RE1');


