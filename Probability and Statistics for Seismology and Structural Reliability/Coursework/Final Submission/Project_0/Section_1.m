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
normplot(s1_data);

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
% tiledlayout(1,2)
% nexttile
scatter(quantile_values, s1_data, 80, 'k', 'filled', 'Marker', 'd'); 
ax = gca;
ax.FontSize = 20;
set(gca, 'TickLabelInterpreter', 'latex');
title('Storey 1','Interpreter', 'latex')
box on;
x.XGrid = 'off'; ax.YGrid = 'on';
xlabel('quantile index','Interpreter', 'latex'); ylabel('Comp. strength (MPa)','Interpreter', 'latex'); ylim([0 50])
% nexttile
figure(4)
scatter(quantile_values, s2_data, 80, 'k', 'filled', 'Marker', 'd');
ax = gca;
set(gca, 'TickLabelInterpreter', 'latex');
ax.FontSize = 20;
box on;
x.XGrid = 'off'; ax.YGrid = 'on';
title('Storey 2','Interpreter', 'latex')
xlabel('quantile index','Interpreter', 'latex'); ylabel('Comp. strength (MPa)','Interpreter', 'latex'); ylim([0 50])

% % Tukey box plot Faber 3.3.4
% combined_data = [s1_data, s2_data];
% figure(4)
% boxplot(combined_data)
% title('Boxplots of s 1 and s 2 (MPa)') % remove from final put in caption
% 
% 
% % Q-Q plots Faber 3.3.5
% % https://uk.mathworks.com/help/stats/qqplot.html Faber 3.3.5
% % how to interpret see p39, look at where the points fall (above/below
% % line)
% figure(5)
% qqplot(s1_data, s2_data); xlabel('Compressive strength, s 1') 
% ylabel('Compressive strength, s 2') 
% title('Scatter Plot with Line of Best Fit')
% 
% % Tukey Mean difference Faber 3.3.5
% tukey_difference = (s1_data - s2_data);
% tukey_mean = (s1_data + s2_data)/2;
% figure(6)
% scatter(tukey_mean, tukey_difference, 20, 'k', 'filled', 'Marker', 'd'); 
% ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
% xlabel('mean')
% ylabel('difference')
% 
% 
% %% SECTION 1.2.1 PROBABILITY PLOT
% 
% % Lognormal Probability plot values (Week 4 activity 6 for parameters) ADD EQ OF
% log_s1 = log(s1_data);
% zeta_s1 = log( 1 + (s1_CoV^2) ) ^ 0.5;
% lambda_s1 = log(s1_mean) - ( (zeta_s1^2)/2 );
% lcdf_s1 = logncdf(s1_data, zeta_s1, lambda_s1);
% loginv_s1 = (1/zeta_s1 * log_s1) - (lambda_s1/zeta_s1);
% p_s1 = polyfit(log_s1, loginv_s1, 1);
% lxfit_s1 = linspace(min(log_s1),max(log_s1),100);
% lyfit_s1 = polyval(p_s1, lxfit_s1);
% 
% log_s2 = log(s2_data);
% zeta_s2 = log( 1 + (s2_CoV^2) ) ^ 0.5;
% lambda_s2 = log(s2_mean) - ( (zeta_s2^2)/2 );
% lcdf_s2 = logncdf(s2_data, zeta_s2, lambda_s2);
% loginv_s2 = (1/zeta_s2 * log_s2) - (lambda_s2/zeta_s2);
% p_s2 = polyfit(log_s2, loginv_s2, 1);
% lxfit_s2 = linspace(min(log_s2),max(log_s2),100);
% lyfit_s2 = polyval(p_s2, lxfit_s2);
% 
% % Normal Probability plot values
% zscore_s1 = (s1_data - s1_mean) / s1_std;
% zscore_s2 = (s2_data - s2_mean) / s2_std;
% n_s1 = polyfit(s1_data, zscore_s1, 1);
% nx_s1 = linspace(min(s1_data),max(s1_data),100);
% ny_s1 = polyval(n_s1, nx_s1);
% n_s2 = polyfit(s2_data, zscore_s2, 1);
% nx_s2 = linspace(min(s2_data),max(s2_data),100);
% ny_s2 = polyval(n_s2, nx_s2);
% %METHOD 2 CAN ALSO DO | CAN I DO THE SAME FOR LOGNORMAL
% % x = normcdf(s1_data, s1_mean, s1_std);
% % y = norminv(x);
% % figure(10);
% % plot(s1_data, y);
% 
% % Normal and Lognormal Plots , the stronger the gradient the better the
% % fit?
% figure(7)
% tiledlayout(1,4)
% % normal
% nexttile
% scatter(s1_data, zscore_s1)
% hold on; grid on; plot(nx_s1, ny_s1);
% xlabel('Storey 1 Data'); ylabel('Inverse N Dist')
% nexttile
% scatter(s2_data, zscore_s2)
% hold on; grid on; plot(nx_s2, ny_s2);
% xlabel('Storey 2 Data'); ylabel('Inverse N Dist')
% % lognormal
% nexttile
% scatter(log_s1, loginv_s1)
% hold on; grid on; plot(lxfit_s1, lyfit_s1)
% xlabel('Log(Storey 1 Data)'); ylabel('Inverse of LogN Dist')
% nexttile
% scatter(log_s2, loginv_s2)
% hold on; grid on; plot(lxfit_s2, lyfit_s2)
% xlabel('Log(Storey 2 Data)'); ylabel('Inverse of LogN Dist')
% 
% 
% %% SECTION 1.2.2 CHI SQUARE
% a = 25;
% zz = 1;
% % STATISTICS, these variables are the same for normal/lognormal
% 
% % Discretize the data into specified edges/bins found using Sturge rule and
% % adjusting when the freq < 5
% edges_s1 = [10, 15.6, 21.2, 26.8, 32.4, 38, Inf];
% edges_s2 = [12, 15.4, 18.8, 22.2, Inf];
% [~, edges_s1] = discretize(s1_data, edges_s1);
% [~, edges_s2] = discretize(s2_data, edges_s2);
% 
% % Calculate observed frequencies
% obs_freq_s1 = histcounts(s1_data, edges_s1);
% obs_freq_s2 = histcounts(s2_data, edges_s2);
% 
% % Degrees of freedom (no of bins - 2, due to mean and std)
% % df = k - 1 - m;  where k is no of bins, m is std + mean
% % is this the same case for normal/lognormal? Clarify
% m = 2; % according to pg 22 wk4 part 5 and pg 124 faber
% df_s1 = (length(edges_s1) - 1) - 1 - m;
% df_s2 = (length(edges_s2) - 1) - 1 - m;
% % Significance level
% alpha = 0.05;
% % Calculate the critical chi-square value for s1
% critical_chi_sq_s1 = chi2inv(1 - alpha, df_s1);
% % Calculate the critical chi-square value for s2
% critical_chi_sq_s2 = chi2inv(1 - alpha, df_s2);
% 
% % NORMAL
% disp('Normal Chi Square')
% 
% % Calculate mean and standard deviation
% mu_norm_s1 = mean(s1_data);
% sigma_norm_s1 = std(s1_data);
% mu_norm_s2 = mean(s2_data);
% sigma_norm_s2 = std(s2_data);
% 
% % Calculate the cumulative probabilities for each bin
% normcdf_s1 = normcdf(edges_s1, mu_norm_s1, sigma_norm_s1);
% normcdf_s2 = normcdf(edges_s2, mu_norm_s2, sigma_norm_s2);
% 
% % Calculate the probabilities for each bin
% normpdf_s1 = diff(normcdf_s1);
% normpdf_s2 = diff(normcdf_s2);
% 
% % Calculate the expected frequencies
% exp_nfreq_s1 = normpdf_s1 * length(s1_data);
% exp_nfreq_s2 = normpdf_s2 * length(s2_data);
% 
% % Perform the Chi-Square test
% chi_sq_norm_s1 = sum((obs_freq_s1 - exp_nfreq_s1).^2 ./ exp_nfreq_s1);
% chi_sq_norm_s2 = sum((obs_freq_s2 - exp_nfreq_s2).^2 ./ exp_nfreq_s2);
% 
% % Compare with the chi-square statistic for s1
% if chi_sq_norm_s1 < critical_chi_sq_s1
%     disp([' Chi-Square Statistic for s1: ' num2str(chi_sq_norm_s1)]);
%     disp([' Critical Chi-Square value for s1 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s1)]);
%     disp(' The null hypothesis cannot be rejected for s1.');
% else
%     disp([' Chi-Square Statistic for s1: ' num2str(chi_sq_norm_s1)]);
%     disp([' Critical Chi-Square value for s1 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s1)]);
%     disp(' The null hypothesis is rejected for s1.');
% end
% 
% % Compare with the chi-square statistic for s2
% if chi_sq_norm_s2 < critical_chi_sq_s2
%     disp([' Chi-Square Statistic for s2: ' num2str(chi_sq_norm_s2)]);
%     disp([' Critical Chi-Square value for s2 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s2)]);
%     disp(' The null hypothesis cannot be rejected for s2.');
% else
%     disp([' Chi-Square Statistic for s2: ' num2str(chi_sq_norm_s2)]);
%     disp([' Critical Chi-Square value for s2 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s2)]);
%     disp(' The null hypothesis is rejected for s2.');
% end
% 
% % Obvserved/expected bar chart
% figure(8);
% % tiledlayout(1,4)
% % nexttile
% bar(edges_s1(1:end-1), [obs_freq_s1', exp_nfreq_s1'], 'grouped', 'BarWidth',zz);
% set(gca, 'TickLabelInterpreter', 'latex');
% ax = gca; ax.FontSize = a;
% ylabel('Frequency', 'Interpreter', 'latex');
% xlabel('Intervals', 'Interpreter', 'latex');
% title('Storey 1 - Normal Model', 'Interpreter', 'latex');
% legend('Observed', 'Expected', 'Location', 'northwest')
% axis square
% xlabel('Intervals'); ylabel('Frequency'); 
% legend('Observed', 'Expected', 'Location', 'southeast','Interpreter', 'latex')
% % Set x-axis ticks and labels
% xticks(edges_s1(1:end-1));
% xticklabels(cellfun(@(a, b) sprintf('%.1f-%.1f', a, b), num2cell(edges_s1(1:end-1)), num2cell(edges_s1(2:end)), 'UniformOutput', false));
% % Show the grid
% grid on;
% 
% figure(9)
% bar(edges_s2(1:end-1), [obs_freq_s2', exp_nfreq_s2'], 'grouped', 'BarWidth',zz);
% set(gca, 'TickLabelInterpreter', 'latex');
% ax = gca; ax.FontSize = a;
% ylabel('Frequency', 'Interpreter', 'latex');
% xlabel('Intervals', 'Interpreter', 'latex');
% title('Storey 2 - Normal Model', 'Interpreter', 'latex');
% legend('Observed', 'Expected', 'Location', 'northwest')
% axis square
% xlabel('Intervals'); ylabel('Frequency'); 
% legend('Observed', 'Expected', 'Location', 'southeast','Interpreter', 'latex')
% % Set x-axis ticks and labels
% xticks(edges_s2(1:end-1));
% xticklabels(cellfun(@(a, b) sprintf('%.1f-%.1f', a, b), num2cell(edges_s2(1:end-1)), num2cell(edges_s2(2:end)), 'UniformOutput', false));
% % Show the grid
% grid on;
% 
% 
% % LOGNORMAL
% disp('Lognormal Chi Square')
% 
% % Lognormal parameters (found using mean/std, see Matlab doc or Wk4 Act 6)
% mu_log_s1 = log(s1_mean^2/(sqrt(s1_std^2 + s1_mean^2))); 
% sigma_log_s1 = sqrt( log( (s1_std^2)/(s1_mean^2) + 1 ) );
% mu_log_s2 = log(s2_mean^2/(sqrt(s2_std^2 + s2_mean^2))); 
% sigma_log_s2 = sqrt( log( (s2_std^2)/(s2_mean^2) + 1 ) );
% 
% % Calculate cdf of each bin
% logncdf_s1 = logncdf(edges_s1, mu_log_s1, sigma_log_s1);
% logncdf_s2 = logncdf(edges_s2, mu_log_s2, sigma_log_s2);
% 
% % Calculate pdf of each bin
% lognpdf_s1 = diff(logncdf_s1);
% lognpdf_s2 = diff(logncdf_s2);
% 
% % Calculate the expected frequencies
% exp_lognfreq_s1 = lognpdf_s1 * length(s1_data);
% exp_lognfreq_s2 = lognpdf_s2 * length(s2_data);
% 
% % Perform the Chi-Square test
% chi_sq_logn_s1 = sum((obs_freq_s1 - exp_lognfreq_s1).^2 ./ exp_lognfreq_s1);
% chi_sq_logn_s2 = sum((obs_freq_s2 - exp_lognfreq_s2).^2 ./ exp_lognfreq_s2);
% 
% % Compare with the chi-square statistic for s1
% if chi_sq_logn_s1 < critical_chi_sq_s1
%     disp([' Chi-Square Statistic for s1: ' num2str(chi_sq_logn_s1)]);
%     disp([' Critical Chi-Square value for s1 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s1)]);
%     disp(' The null hypothesis cannot be rejected for s1.');
% else
%     disp([' Chi-Square Statistic for s1: ' num2str(chi_sq_logn_s1)]);
%     disp([' Critical Chi-Square value for s1 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s1)]);
%     disp(' The null hypothesis is rejected for s1.');
% end
% 
% % Compare with the chi-square statistic for s2
% if chi_sq_logn_s2 < critical_chi_sq_s2
%     disp([' Chi-Square Statistic for s2: ' num2str(chi_sq_logn_s2)]);
%     disp([' Critical Chi-Square value for s2 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s2)]);
%     disp(' The null hypothesis cannot be rejected for s2.');
% else
%     disp([' Chi-Square Statistic for s2: ' num2str(chi_sq_logn_s2)]);
%     disp([' Critical Chi-Square value for s2 at ' num2str(alpha*100) '% significance: ' num2str(critical_chi_sq_s2)]);
%     disp(' The null hypothesis is rejected for s2.');
% end
% 
% % Obvserved/expected bar chart
% 
% figure(10);
% bar(edges_s1(1:end-1), [obs_freq_s1', exp_lognfreq_s1'], 'grouped', 'BarWidth',zz);
% set(gca, 'TickLabelInterpreter', 'latex');
% ax = gca; ax.FontSize = a;
% ylabel('Frequency', 'Interpreter', 'latex');
% xlabel('Intervals', 'Interpreter', 'latex');
% title('Storey 1 - Lognormal Model', 'Interpreter', 'latex');
% legend('Observed', 'Expected', 'Location', 'northwest')
% axis square
% xlabel('Intervals'); ylabel('Frequency'); 
% legend('Observed', 'Expected', 'Location', 'southeast','Interpreter', 'latex')
% % Set x-axis ticks and labels
% xticks(edges_s1(1:end-1));
% xticklabels(cellfun(@(a, b) sprintf('%.1f-%.1f', a, b), num2cell(edges_s1(1:end-1)), num2cell(edges_s1(2:end)), 'UniformOutput', false));
% % Show the grid
% 
% grid on;
% 
% figure(11)
% bar(edges_s2(1:end-1), [obs_freq_s2', exp_lognfreq_s2'], 'grouped', 'BarWidth',zz);
% set(gca, 'TickLabelInterpreter', 'latex');
% ax = gca; ax.FontSize = a;
% ylabel('Frequency', 'Interpreter', 'latex');
% xlabel('Intervals', 'Interpreter', 'latex');
% title('Storey 2 - Lognormal Model', 'Interpreter', 'latex');
% legend('Observed', 'Expected', 'Location', 'northwest')
% axis square
% xlabel('Intervals'); ylabel('Frequency'); 
% legend('Observed', 'Expected', 'Location', 'southeast','Interpreter', 'latex')
% % Set x-axis ticks and labels
% xticks(edges_s2(1:end-1));
% xticklabels(cellfun(@(a, b) sprintf('%.1f-%.1f', a, b), num2cell(edges_s2(1:end-1)), num2cell(edges_s2(2:end)), 'UniformOutput', false));
% % Show the grid
% grid on;
% 
% %% SECTION 1.2.3 Kolmogorov and Smirnov goodness of fit tests.
% 
% % Statistics
% 
% % Could be s2_data as they are the same length
% N = length(s1_data);
% % Creation of the SEED vector: 1,2,...,N
% SEED = 1:1:N;
% % Calculation of the percentiles
% ks_obs = (SEED./(N+1)).';
% 
% % NORMAL
% 
% % Predicted cdf values
% ksnorm_exp_s1 = normcdf(s1_data, s1_mean, s1_std);
% ksnorm_exp_s2 = normcdf(s2_data, s2_mean, s2_std);
% 
% % Errors
% norm_errors_s1 = ks_obs - ksnorm_exp_s1;
% normmax_error_s1 = max(norm_errors_s1);
% norm_errors_s2 = ks_obs - ksnorm_exp_s2;
% normmax_error_s2 = max(norm_errors_s2);
% 
% % Hypothesis
% disp('Normal K-S')
% disp([' For s1 (Sample Size N = ' num2str(N) '): Max error = ' num2str(normmax_error_s1)]);
% disp([' For s2 (Sample Size N = ' num2str(N) '): Max error = ' num2str(normmax_error_s2)]);
% 
% % LOGNORMAL
% 
% % Predicted cdf values
% kslog_exp_s1 = logncdf(s1_data, mu_log_s1, sigma_log_s1);
% kslog_exp_s2 = logncdf(s2_data, mu_log_s2, sigma_log_s2);
% 
% % Errors
% log_errors_s1 = ks_obs - kslog_exp_s1;
% logmax_error_s1 = max(log_errors_s1);
% log_errors_s2 = ks_obs - kslog_exp_s2;
% logmax_error_s2 = max(log_errors_s2);
% 
% % Hypothesis
% disp('Lognormal K-S')
% disp([' For s1 (Sample Size N = ' num2str(N) '): Max error = ' num2str(logmax_error_s1)]);
% disp([' For s2 (Sample Size N = ' num2str(N) '): Max error = ' num2str(logmax_error_s2)]);
% 
% %% 1.3 Sample Likelihood
% 
% disp('Sample Likelihood (greater means better fit)')
% 
% % NORMAL S1
% normal_like_s1 = sampleLike(chi_sq_norm_s1, df_s1);
% % LOGNORMAL S1
% lognormal_like_s1 = sampleLike(chi_sq_logn_s1, df_s1);
% disp([' For s1 the sample likelihood for the normal is ' num2str(normal_like_s1) ', lognormal is ' num2str(lognormal_like_s1)] );
% 
% % NORMAL S2
% normal_like_s2 = sampleLike(chi_sq_norm_s2, df_s2);
% % LOGNORMAL S2
% lognormal_like_s2 = sampleLike(chi_sq_logn_s2, df_s2);
% disp([' For s2 the sample likelihood for the normal is ' num2str(normal_like_s2) ', lognormal is ' num2str(lognormal_like_s2)] );
% 
% function likelihood = sampleLike(cs, df)
%     likelihood = ( (cs^(df/(2-1))) / (2^(df/2) * gamma(df/2)) ) * (exp(-cs-2));
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
