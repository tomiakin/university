clear
clc
close all

% THE DATA SET IS CONCRETE CUBE COMPRESSIVE STRENGTHS (MPa)

%% SECTION 1.1
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

% See Faber 3.2.6 Sample moments 

% Measures of correlation (Faber 3.2.7)
disp('Measures of correlation')
sample_covariance = cov(s1_data, s2_data);
disp([' The sample covariance of s 1 and 2 or cov(1,2) is: ', num2str(sample_covariance(1,2))])
sample_coeff = corrcoef(s1_data, s2_data);
disp([' The sample coefficient of correlation of s 1 and 2 is: ', num2str(sample_coeff(1,2))])

% Graphical representations (Faber 3.3)

% 1D scatter plot (Faber 3.3.1) 
% PERHAPS MARK OUT CENTRAL MEASURES AND THEN USE THE
% GRAPH TO DISCUSS CENTRAL MEASURES ETC see Faber 3.3.1 for when to use
% (data set should be small so may not be applicable here!)
figure(1)
scatter(s1_data, 2 * ones(size(s1_data)), 50, 'k^', 'DisplayName', 's 1 Data');
hold on;
scatter(s2_data, -2 * ones(size(s2_data)), 50, 'ko', 'DisplayName', 's 2 Data');
% Set the y-axis limits
ylim([-5, 5]);
set(gca, 'YTick', [], 'YColor', 'none');
legend('s 1', 's 2')
title('1D Scatter Plot') % remove from final plot

% Histograms (Faber 3.3.2)
% USE when observations are unknown but freq is known so may not
% be suitable here Faber 3.3.2! so you get k as no iof intervals why not
% then use integer bins like 10-20 and 30- like you will do fo rchi square
% also note this is sturges rule reference it
k = int32 (1 + ( 3.3 * log10(length(concrete_data)) ) ); % k is no of intervals
figure(2)
tiledlayout(2,1)
nexttile
histogram(s1_data, k); title('s 1'); xlabel('MPa'); ylabel('Number of observations');
nexttile
histogram(s2_data, k); title('s 2'); xlabel('MPa'); ylabel('Number of observations');



% Quantile Plots (Faber 3.3.3)
% Discuss interpolation?
% and extraction of median, see p35 on how to describe stuff its key!
% see paragraph under eq 3.9 says why this is favourable over hist
% see p35 for a nice way to format figure, maybe set the size of each one
% individually so its scaled nicely, see how they wrote caption
figure(3)
quantile_values = data_id/(length(data_id)+1); %eq 3.9 Faber
tiledlayout(2,1)
nexttile
scatter(quantile_values, s1_data, 20, 'k', 'filled', 'Marker', 'd'); 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
title('s 1'); xlabel('quantile index, \nu'); ylabel('Compressive strength (MPa)'); ylim([0 50])
nexttile
scatter(quantile_values, s2_data, 20, 'k', 'filled', 'Marker', 'd');  
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
title('s 2'); xlabel('quantile index, \nu'); ylabel('Compressive strength (MPa)'); ylim([0 50])

% Tukey box plot Faber 3.3.4
combined_data = [s1_data, s2_data];
figure(4)
boxplot(combined_data)
title('Boxplots of s 1 and s 2 (MPa)') % remove from final put in caption


% Q-Q plots Faber 3.3.5
% https://uk.mathworks.com/help/stats/qqplot.html Faber 3.3.5
% how to interpret see p39, look at where the points fall (above/below
% line)
figure(5)
qqplot(s1_data, s2_data); xlabel('Compressive strength, s 1') 
ylabel('Compressive strength, s 2') 
title('Scatter Plot with Line of Best Fit')

% Tukey Mean difference Faber 3.3.5
tukey_difference = (s1_data - s2_data);
tukey_mean = (s1_data + s2_data)/2;
figure(6)
scatter(tukey_mean, tukey_difference, 20, 'k', 'filled', 'Marker', 'd'); 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
xlabel('mean')
ylabel('difference')


%% SECTION 1.2

% Lognormal Probability plot values (Week 4 activity 6 for parameters) ADD EQ OF
log_s1 = log(s1_data);
zeta_s1 = log( 1 + (s1_CoV^2) ) ^ 0.5;
lambda_s1 = log(s1_mean) - ( (zeta_s1^2)/2 );
lcdf_s1 = logncdf(s1_data, zeta_s1, lambda_s1);
loginv_s1 = (1/zeta_s1 * log_s1) - (lambda_s1/zeta_s1);
p_s1 = polyfit(log_s1, loginv_s1, 1);
lxfit_s1 = linspace(min(log_s1),max(log_s1),100);
lyfit_s1 = polyval(p_s1, lxfit_s1);

log_s2 = log(s2_data);
zeta_s2 = log( 1 + (s2_CoV^2) ) ^ 0.5;
lambda_s2 = log(s2_mean) - ( (zeta_s2^2)/2 );
lcdf_s2 = logncdf(s2_data, zeta_s2, lambda_s2);
loginv_s2 = (1/zeta_s2 * log_s2) - (lambda_s2/zeta_s2);
p_s2 = polyfit(log_s2, loginv_s2, 1);
lxfit_s2 = linspace(min(log_s2),max(log_s2),100);
lyfit_s2 = polyval(p_s2, lxfit_s2);

% Normal Probability plot values
zscore_s1 = (s1_data - s1_mean) / s1_std;
zscore_s2 = (s2_data - s2_mean) / s2_std;
n_s1 = polyfit(s1_data, zscore_s1, 1);
nx_s1 = linspace(min(s1_data),max(s1_data),100);
ny_s1 = polyval(n_s1, nx_s1);
n_s2 = polyfit(s2_data, zscore_s2, 1);
nx_s2 = linspace(min(s2_data),max(s2_data),100);
ny_s2 = polyval(n_s2, nx_s2);
%METHOD 2 CAN ALSO DO | CAN I DO THE SAME FOR LOGNORMAL
% x = normcdf(s1_data, s1_mean, s1_std);
% y = norminv(x);
% figure(10);
% plot(s1_data, y);

% Normal and Lognormal Plots , the stronger the gradient the better the
% fit?
figure()
tiledlayout(2,2)
% normal
nexttile
scatter(s1_data, zscore_s1)
hold on; grid on; plot(nx_s1, ny_s1);
xlabel('Storey 1 Data'); ylabel('Inverse of normal distribution')
nexttile
scatter(s2_data, zscore_s2)
hold on; grid on; plot(nx_s2, ny_s2);
xlabel('Storey 2 Data'); ylabel('Inverse of normal distribution')
% lognormal
nexttile
scatter(log_s1, loginv_s1)
hold on; grid on; plot(lxfit_s1, lyfit_s1)
xlabel('Log(Storey 1 Data)'); ylabel('Inverse of lognormal distribution')
nexttile
scatter(log_s2, loginv_s2)
hold on; grid on; plot(lxfit_s2, lyfit_s2)
xlabel('Log(Storey 2 Data)'); ylabel('Inverse of lognormal distribution')


% Chi Square

% Number of bins
num_bins = k;

% Discretize the data into 7 bins
[~, edges_s1] = discretize(s1_data, num_bins);
[~, edges_s2] = discretize(s2_data, num_bins);

% Calculate observed frequencies
obs_freq_s1 = histcounts(s1_data, edges_s1);
obs_freq_s2 = histcounts(s2_data, edges_s2);

% Calculate mean and standard deviation
mu_s1 = mean(s1_data);
sigma_s1 = std(s1_data);

mu_s2 = mean(s2_data);
sigma_s2 = std(s2_data);

% Calculate the cumulative probabilities for each bin
cum_prob_s1 = normcdf(edges_s1, mu_s1, sigma_s1);
cum_prob_s2 = normcdf(edges_s2, mu_s2, sigma_s2);

% Calculate the probabilities for each bin
prob_s1 = diff(cum_prob_s1);
prob_s2 = diff(cum_prob_s2);

% Calculate the expected frequencies
exp_freq_s1 = prob_s1 * length(s1_data);
exp_freq_s2 = prob_s2 * length(s2_data);

% Perform the Chi-Square test
chi_sq_stat_s1 = sum((obs_freq_s1 - exp_freq_s1).^2 ./ exp_freq_s1);
chi_sq_stat_s2 = sum((obs_freq_s2 - exp_freq_s2).^2 ./ exp_freq_s2);

% Degrees of freedom
df = num_bins - 1;









