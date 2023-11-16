%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

%% How to simulate from an uniform distribution
%  The extrem of the distribution are a and b
a = 1;
b = 10;

% N is the number of simulations
N = 100000;

sample = a + (b-a)*rand(N,1);

hist(sample)
xlabel('Values')
ylabel('Frequency')