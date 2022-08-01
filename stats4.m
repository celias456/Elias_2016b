function [data] = stats4(allvariables)
%Calculates out of sample forecasting statistics
%Input: 
% allvariables: 
%     allvariables(:,1) = fundamentals
%     allvariables(:,2) = exchange rate (logs)
%     allvariables(:,3) = deviation from fundamentals


%Output:
% data: row vector of statistics
% U-stat (1 period ahead)
% U-stat (4 period ahead)
% U-stat (8 period ahead)
% U-stat (16 period ahead)

%% Get variables
s = allvariables(:,2); %Get Log exchange rate
xi = allvariables(:,3); %Get deviation from fundamentals

T = length(s);

%Calculate estimate of drift parameter (delta-hat)
s_available = s(1:52); %log exchange rate
delta_hat = (s_available(52) - s_available(1))/51; % estimate of drift parameter

%% Calculate parameters of monetary model

%1-period ahead 
delta_s_1 = kaheadchange(s,1); %1 period ahead change of log exchange rate
x_1 = [ones(51,1),xi(1:51)]; %constant and deviation from fundamentals 
y_1 = delta_s_1(1:51); %"available" exchange rate depreciation data
b_1 = regress(y_1,x_1); %conduct regression

%4-period ahead 
delta_s_4 = kaheadchange(s,4); %4 period ahead change of log exchange rate
x_4 = [ones(48,1),xi(1:48)]; %constant and deviation from fundamentals 
y_4 = delta_s_4(1:48); %"available" exchange rate depreciation data
b_4 = regress(y_4,x_4); %conduct regression

%8-period ahead 
delta_s_8 = kaheadchange(s,8); %8 period ahead change of log exchange rate
x_8 = [ones(44,1),xi(1:44)]; %constant and deviation from fundamentals 
y_8 = delta_s_8(1:44); %"available" exchange rate depreciation data
b_8 = regress(y_8,x_8); %conduct regression

%16-period ahead 
delta_s_16 = kaheadchange(s,16); %16 period ahead change of log exchange rate
x_16 = [ones(36,1),xi(1:36)]; %constant and deviation from fundamentals 
y_16 = delta_s_16(1:36); %"available" exchange rate depreciation data
b_16 = regress(y_16,x_16); %conduct regression

%% Forecast errors, RMSE, and U-Statistic

%1-period ahead
%Get length of "future" data series
T_future_1 = T-52;
%Get actual depreciation
actual_dep_1 = delta_s_1(52:end); %Actual depreciation
%rw
dep_rw_1 = delta_hat.*ones(T_future_1,1); %RW forecast of deprecation
error_rw_1 = (actual_dep_1 - dep_rw_1).^2; %squared forecast error
rmse_rw_1 = sqrt(mean(error_rw_1)); %RMSE
%monetary model
alpha_1 = b_1(1).*ones(T_future_1,1); %make vector of constant parameters
xi_1 = xi(52:(end-1)); %Get deviation from fundamentals data used for monetary model forecast
dep_mm_1 = alpha_1 + b_1(2).*xi_1; %monetary model forecast of deprecation
error_mm_1 = (actual_dep_1 - dep_mm_1).^2; %squared forecast error
rmse_mm_1 = sqrt(mean(error_mm_1)); %RMSE
%Calculate U statistic
U_1 = rmse_mm_1/rmse_rw_1;

%4-period ahead
%Get length of "future" data series
T_future_4 = T-55;
%Get actual depreciation
actual_dep_4 = delta_s_4(52:end); %Actual depreciation
%rw
dep_rw_4 = 4*delta_hat.*ones(T_future_4,1); %RW forecast of deprecation
error_rw_4 = (actual_dep_4 - dep_rw_4).^2; %squared forecast error
rmse_rw_4 = sqrt(mean(error_rw_4)); %RMSE
%monetary model
alpha_4 = b_4(1).*ones(T_future_4,1); %make vector of constant parameters
xi_4 = xi(52:(end-4)); %Get deviation from fundamentals data used for monetary model forecast
dep_mm_4 = alpha_4 + b_4(2).*xi_4; %monetary model forecast of deprecation
error_mm_4 = (actual_dep_4 - dep_mm_4).^2; %squared forecast error
rmse_mm_4 = sqrt(mean(error_mm_4)); %RMSE
%Calculate U statistic
U_4 = rmse_mm_4/rmse_rw_4;


%8-period ahead
%Get length of "future" data series
T_future_8 = T-59;
%Get actual depreciation
actual_dep_8 = delta_s_8(52:end); %Actual depreciation
%rw
dep_rw_8 = 8*delta_hat.*ones(T_future_8,1); %RW forecast of deprecation
error_rw_8 = (actual_dep_8 - dep_rw_8).^2; %squared forecast error
rmse_rw_8 = sqrt(mean(error_rw_8)); %RMSE
%monetary model
alpha_8 = b_8(1).*ones(T_future_8,1); %make vector of constant parameters
xi_8 = xi(52:(end-8)); %Get deviation from fundamentals data used for monetary model forecast
dep_mm_8 = alpha_8 + b_8(2).*xi_8; %monetary model forecast of deprecation
error_mm_8 = (actual_dep_8 - dep_mm_8).^2; %squared forecast error
rmse_mm_8 = sqrt(mean(error_mm_8)); %RMSE
%Calculate U statistic
U_8 = rmse_mm_8/rmse_rw_8;


%16-period ahead
%Get length of "future" data series
T_future_16 = T-67;
%Get actual depreciation
actual_dep_16 = delta_s_16(52:end); %Actual depreciation
%rw
dep_rw_16 = 16*delta_hat.*ones(T_future_16,1); %RW forecast of deprecation
error_rw_16 = (actual_dep_16 - dep_rw_16).^2; %squared forecast error
rmse_rw_16 = sqrt(mean(error_rw_16)); %RMSE
%monetary model
alpha_16 = b_16(1).*ones(T_future_16,1); %make vector of constant parameters
xi_16 = xi(52:(end-16)); %Get deviation from fundamentals data used for monetary model forecast
dep_mm_16 = alpha_16 + b_16(2).*xi_16; %monetary model forecast of deprecation
error_mm_16 = (actual_dep_16 - dep_mm_16).^2; %squared forecast error
rmse_mm_16 = sqrt(mean(error_mm_16)); %RMSE
%Calculate U statistic
U_16 = rmse_mm_16/rmse_rw_16;

data = [U_1,U_4,U_8,U_16];
end

