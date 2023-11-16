% Define the PDF function
pdf = @(x) x / 210;

% Create an array of the first 20 real numbers from 1 to 20
x_values = 1:20;

% Calculate the PDF and CDF values for each number DISCRETE IS CUMSUM?
pdf_values = pdf(x_values);
cdf_values = cumsum(pdf_values);

% Plot PDF and CDF
figure;
subplot(1, 2, 1);
bar(x_values, pdf_values, 'b');
xlabel('X');
ylabel('PDF');
title('Probability Density Function (PDF)');

subplot(1, 2, 2);
plot(x_values, cdf_values, 'ro-');
xlabel('X');
ylabel('CDF');
title('Cumulative Distribution Function (CDF)');

% Finding P(X > 17)
cdf_17th = cdf_values(17);  % MATLAB uses 1-based indexing
prob = 1 - cdf_17th;

fprintf('P(X > 17) is %.4f\n', prob);



