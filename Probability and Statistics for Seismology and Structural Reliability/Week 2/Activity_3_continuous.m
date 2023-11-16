% Define the PDF function as a piecewise function:
%   - It calculates 2x for values of x within the interval [0, 1].
%   - It restricts the range or sets the PDF to 0 for values of x outside the interval [0, 1].
pdf = @(x) 2 * x .* (x >= 0 & x <= 1);

% Calculate P(1/4 <= X <= 1/2) CONTINUOUS IS AN INTEGRAL?
probability = integral(pdf, 1/4, 1/2);

fprintf('P(1/4 <= X <= 1/2) is approximately: %.4f\n', probability);

% Calculate the integral of the PDF over the entire range [0, 1]
total_integral = integral(pdf, 0, 1);

fprintf('The integral of the PDF over the entire range [0, 1] is approximately: %.4f\n', total_integral);

% Plot the PDF
x_values = linspace(0, 1, 100);
pdf_values = pdf(x_values);

figure;
plot(x_values, pdf_values, 'b', 'LineWidth', 2);
xlabel('x');
ylabel('PDF');
title('Probability Density Function (PDF)');
grid on;
