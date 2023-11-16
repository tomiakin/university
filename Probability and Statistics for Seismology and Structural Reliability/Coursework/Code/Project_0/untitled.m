clear
clc
close all

% Sample data
data = [9.2, 12.5, 16.8, 18.7, 21.3, 24.6, 28.1, 31.4, 33.7, 36.9, 41.2, 45.6, 49.8];

% Specify the number of bins
num_bins = 7;

% Calculate the data range
data_min = min(data);
data_max = max(data);

% Round the minimum and maximum to integers
data_min = floor(data_min);
data_max = ceil(data_max);

% Calculate the bin width to fit 7 bins
bin_width = (data_max - data_min) / num_bins;

% Calculate the bin edges for 7 bins with integer values
bin_edges = data_min:bin_width:data_max;

% Create the histogram with custom bin edges
figure;
histogram(data, bin_edges);
title('Custom Histogram with Integer Bins');
xlabel('Value');
ylabel('Frequency');
