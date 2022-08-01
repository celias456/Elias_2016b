%Graph Code
%Order of variables in exchangerate.csv:
%1. period
%2. exchange rate level (RE1)
%3. exchange rate level (SN)
%4. deviation from fundamental log (RE1)
%5. deviation from fundamental log (SN)
%6. % Deviation of actual exchange rate (bunch of zeros) (SN)
%7. % Deviation of Speculator forecast (SN)
%8. % Deviation of Non-speculator forecast (SN)
%9. Speculator Beliefs constant
%10. Speculator Beliefs time
%11. Speculator Beliefs fundamentals
%12. Non-Speculator Beliefs constant
%13. Non-Speculator Beliefs lagged exchange rate

%% Exchange Rate Graph
load exchangerate.csv
exchangerate_data = exchangerate;

x1 = exchangerate_data(:,1);
y1 = exchangerate_data(:,2);
y2 = exchangerate_data(:,3);

plot(x1,y1,'k',x1,y2,'--k')
xlabel('Period')
ylabel('Exchange Rate')

legend('RE1','SN')

%% Deviation From Fundamentals Graph

load exchangerate.csv
exchangerate_data = exchangerate;

x1 = exchangerate_data(:,1);
y1 = exchangerate_data(:,4);
y2 = exchangerate_data(:,5);

plot(x1,y1,'k',x1,y2,'--k')
xlabel('Period')
ylabel('Deviation From Fundamentals (log)')

legend('RE1','AP')

%% Forecast Graph

load exchangerate.csv
exchangerate_data = exchangerate;

x1 = exchangerate_data(:,1);
y1 = exchangerate_data(:,6);
y2 = exchangerate_data(:,7);
y3 = exchangerate_data(:,8);

plot(x1,y1,'k',x1,y2,'--k',x1,y3,':k')
xlabel('Period')
ylabel('Percent Deviation')

legend('Exchange Rate','Speculator','Non-Speculator')

%% Beliefs

load exchangerate.csv
exchangerate_data = exchangerate;

x1 = exchangerate_data(:,1);
y1 = exchangerate_data(:,9);
y2 = exchangerate_data(:,10);
y3 = exchangerate_data(:,11);
y4 = exchangerate_data(:,12);
y5 = exchangerate_data(:,13);

figure
subplot(2,1,1)
plot(x1,y1,'k',x1,y2,'--k',x1,y3,':k')
xlabel('Period')
ylabel('Parameter Value')
% legend('Constant','Price','Prod. Shock')
h1 = legend('Constant','$t$','$f_{t}$');
set(h1,'Interpreter','latex')
title('Speculator')
% axis([0,300,-1,3])

subplot(2,1,2)
plot(x1,y4,'k',x1,y5,'--k')
xlabel('Period')
ylabel('Parameter Value')
% legend('Constant','Price')
h1 = legend('Constant','$s_{t-1}$');
set(h1,'Interpreter','latex')
title('Non-Speculator')
axis([0,150,0,1.6])
