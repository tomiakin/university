%% SECTION 1.1
close all;
clear;
clc;
% You've done plot now see if they are joint/individual
% Discuss everthing you can take from plot, and take most useful ones

% Upload data from excel file
concrete_sample = readtable("0_Data_Section_1_and_2.xlsx", ...
    'Sheet','Data_Section_1', 'Range', 'A2:C64');
concrete_data = table2array(concrete_sample);

% Sorting of the data in ascending order
data_id = sort(concrete_data(:,1));
s1_data = sort(concrete_data(:,2));
s2_data = sort(concrete_data(:,3));

% Central measures (Faber 3.2.1)
disp('Central Measures')
% mean
s1_mean = mean(s1_data);
s2_mean = mean(s2_data);
disp([' The sample mean of s 1 is: ', num2str(s1_mean)])
disp([' The sample mean of s 2 is: ', num2str(s2_mean)])
% median
s1_median = median(s1_data);
s2_median = median(s2_data);
disp([' The sample median of s 1 is: ', num2str(s1_median)])
disp([' The sample median of s 2 is: ', num2str(s2_median)])
% mode (see Week 2.1 Part 1 p11 and Faber p30)

% Dispersion measures (Faber 3.2.4)
disp('Dispersion Measures')
% variance
s1_var = var(s1_data);
s2_var = var(s2_data);
disp([' The sample variance of s 1 is: ', num2str(s1_var)])
disp([' The sample variance of s 2 is: ', num2str(s2_var)])
% standard dev can also be found through (sqrt(s1_var));
s1_std = std(s1_data);
s2_std = std(s2_data);
disp([' The standard deviation of s 1 is: ', num2str(s1_std)])
disp([' The standard deviation of s 2 is: ', num2str(s2_std)])
% coefficient of Variation
s1_CoV = s1_std/s1_mean;
s2_CoV = s1_std/s2_mean;
disp([' The sample CoV of s 1 is: ', num2str(s1_CoV)])
disp([' The sample CoV of s 2 is: ', num2str(s2_CoV)])

% Other Measures (Faber 3.2.5)
disp('Other Measures')
% skewness
s1_skewness = skewness(s1_data);
s2_skewness = skewness(s2_data);
disp([' The skewness of s 1 is: ', num2str(s1_skewness)])
disp([' The skewness of s 2 is: ', num2str(s2_skewness)])
% kurtosis
s1_kurtosis = kurtosis(s1_data);
s2_kurtosis = kurtosis(s2_data);
disp([' The kurtosis of s 1 is: ', num2str(s1_kurtosis)])
disp([' The kurtosis of s 2 is: ', num2str(s2_kurtosis)])