clear
clc
close all

%% Your data
data = [10, 15, 24, 30, 35];

% Specify the number of bins
numBins = 3;

% Bin the data
[counts, edges] = histcounts(data, numBins);

% Display the counts in each bin
disp('Bin Counts:');
disp(counts);