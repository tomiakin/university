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

%% SECTION 1.2.1 PROBABILITY PLOT

% General

% Quantification of the number of data N
N = length(s1_data);

% Creation of the SEED vector: 1,2,...,N
SEED = 1:1:N;

% Calculation of the percentiles
Percentiles = SEED./(1+N);

% Normal
s1_norm_inv = norminv(Percentiles);
s1_norm_linear = fittype('c + m*x', 'coefficients', {'c', 'm'});
s1_norm_fit = fit(s1_norm_inv', s1_data, s1_norm_linear);
n_x_s1_fit = linspace(min(s1_norm_inv), max(s1_norm_inv), 100);
n_y_s1_fit = feval(s1_norm_fit, n_x_s1_fit);

s2_norm_inv = norminv(Percentiles);
s2_norm_linear = fittype('c + m*x', 'coefficients', {'c', 'm'});
s2_norm_fit = fit(s2_norm_inv', s2_data, s2_norm_linear);
n_x_s2_fit = linspace(min(s2_norm_inv), max(s2_norm_inv), 100);
n_y_s2_fit = feval(s2_norm_fit, n_x_s2_fit);

% Lognormal
s1_log_inv = logninv(Percentiles);
s1_log_linear = fittype('c + m*x', 'coefficients', {'c', 'm'});
s1_log_fit = fit(s1_norm_inv', log(s1_data), s1_log_linear);
l_x_s1_fit = linspace(min(s1_norm_inv), max(s1_norm_inv), 100);
l_y_s1_fit = feval(s1_log_fit, l_x_s1_fit);

s2_log_inv = logninv(Percentiles);
s2_log_linear = fittype('c + m*x', 'coefficients', {'c', 'm'});
s2_log_fit = fit(s2_norm_inv', log(s2_data), s2_log_linear);
l_x_s2_fit = linspace(min(s2_norm_inv), max(s2_norm_inv), 100);
l_y_s2_fit = feval(s2_log_fit, l_x_s2_fit);

% Plot
% figure;
% tiledlayout(1,4)
a = 20;
z = 13;
l = 3;

% Storey 1 - Normal Model
% nexttile
figure(1)
plot(s1_norm_inv, s1_data, 'ks','markersize',z,'linewidth',2,'markerfacecolor','g');
hold on;
plot(n_x_s1_fit, n_y_s1_fit, 'r--', 'LineWidth', l);
coefficients_s1 = coeffvalues(s1_norm_fit);
equation_str = sprintf('$x = %.2f %+0.2f \\phi^{-1} F_x(x)$',coefficients_s1(1), coefficients_s1(2));
text_pos_s1 = [0.05, 0.9];  % Adjust this position based on your preference
text(text_pos_s1(1), text_pos_s1(2), equation_str, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', a);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = a;
ylabel('x - Sorted Data', 'Interpreter', 'latex');
xlabel('$$\phi^{-1} F_{x}(x)$$', 'Interpreter', 'latex');
title('Storey 1 - Normal Model', 'Interpreter', 'latex');
axis square
grid on

% Storey 2 - Normal Model
% nexttile
figure(2)
plot(s2_norm_inv, s2_data, 'ks','markersize',z,'linewidth',2,'markerfacecolor','g');
hold on;
plot(n_x_s2_fit, n_y_s2_fit, 'r--', 'LineWidth', l);
coefficients_s2 = coeffvalues(s2_norm_fit);
equation_str = sprintf('$x = %.2f %+0.2f \\phi^{-1} F_x(x)$',coefficients_s2(1), coefficients_s2(2));
text_pos_s2 = [0.05, 0.9];  % Adjust this position based on your preference
text(text_pos_s2(1), text_pos_s2(2), equation_str, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', a);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = a;
ylabel('x - Sorted Data', 'Interpreter', 'latex');
xlabel('$$\phi^{-1} F_{x}(x)$$', 'Interpreter', 'latex');
title('Storey 2 - Normal Model', 'Interpreter', 'latex');
axis square
grid on

% Storey 1 - Lognormal Model
% nexttile
figure(3)
plot(s1_norm_inv, log(s1_data), 'ko','markersize',z,'linewidth',2,'markerfacecolor','y')
hold on;
plot(l_x_s1_fit, l_y_s1_fit, 'r--', 'LineWidth', l);
coefficients_ls1 = coeffvalues(s1_log_fit);
equation_str = sprintf('$$\\ln(x) = %0.2f %+0.2f \\phi^{-1} F_x(x)$$', coefficients_ls1(1), coefficients_ls1(2));
text_pos_ls1 = [0.05, 0.9];  % Adjust this position based on your preference
text(text_pos_ls1(1), text_pos_ls1(2), equation_str, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', a);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = a;
ylabel('$\ln(x)$', 'Interpreter', 'latex');
xlabel('$$\phi^{-1} F_{x}(x)$$', 'Interpreter', 'latex');
title('Storey 1 - Lognormal Model', 'Interpreter', 'latex');
axis square
grid on

% Storey 2 - Lognormal Model
% nexttile
figure(4)
plot(s2_norm_inv, log(s2_data), 'ko','markersize',z,'linewidth',2,'markerfacecolor','y')
hold on;
plot(l_x_s2_fit, l_y_s2_fit, 'r--', 'LineWidth', l);
coefficients_ls2 = coeffvalues(s2_log_fit);
equation_str = sprintf('$$\\ln(x) = %0.2f %+0.2f \\phi^{-1} F_x(x)$$', coefficients_ls2(1), coefficients_ls2(2));
text_pos_ls2 = [0.05, 0.9];  % Adjust this position based on your preference
text(text_pos_ls2(1), text_pos_ls2(2), equation_str, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', a);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = a;
ylabel('$\ln(x)$', 'Interpreter', 'latex');
xlabel('$$\phi^{-1} F_{x}(x)$$', 'Interpreter', 'latex');
title('Storey 2 - Lognormal Model', 'Interpreter', 'latex');
axis square
grid on




