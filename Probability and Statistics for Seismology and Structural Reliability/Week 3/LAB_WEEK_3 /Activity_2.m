clear
clc
close all

%% Activity_2
% Data
number_of_samples = 100;
mean_data = 2.346; %mm
std_data = 0.047;  %mm

% 90% confidence interval of the mean?

lowerbound = norminv(0.05);
upperbound = norminv(0.95);

lowerbounb_mean = mean_data + lowerbound*std_data/sqrt(number_of_samples);
upperbound_mean = mean_data + upperbound*std_data/sqrt(number_of_samples);

disp(['The lower bound of the mean is = ',num2str(lowerbounb_mean)])
disp(['The upper bound of the mean is = ',num2str(upperbound_mean)])


