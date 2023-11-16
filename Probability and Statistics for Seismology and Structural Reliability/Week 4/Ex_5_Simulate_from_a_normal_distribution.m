%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

%% How to simulate from an normal distribution
%  The normal distribution is complitely defined by the mean and the
%  standard deviation
mean_value = 10;
std_value  = 2.5;

% N is the number of simulations
N = 100000;

sample = mean_value + std_value*randn(N,1);

hist(sample,100) % 100 columns
xlabel('Values')
ylabel('Frequency')