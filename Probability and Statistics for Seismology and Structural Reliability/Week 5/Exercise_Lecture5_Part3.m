clear
clc
%% R - Lognormal
muR  = 2300; %kN
COVR = 0.13;

%% D - Normal
muD  = 900; %kN
COVD = 0.10;

%% L - Gumble
muL  = 675; %kN
COVL = 0.25;

%% Step 1: conversion of the available statistics to the proper statistics
%          for the corresponding distributions

% Deviazione standard dei logaritmi
BetaR = sqrt(log(1+COVR^2));

% Logarithm of the median
etaR = log(muR) - 0.5*BetaR^2;

% First parameter for the Gumble
a = pi/sqrt(6)/(COVL*muL);

% Second paramter for the Gumble
u = muL - 0.5777/a;

%% Step 2: generation of the random sampling
n = 100000; % number of simulation to perform

r = exp(etaR+BetaR.*randn(n,1));
d = muD + muD*COVD.*randn(n,1);
l = u - log(-log(rand(n,1)))/a;

g = r-d-l;

%% Output
figure(1)
subplot(131)
hist(r); title('Resistance'); xlabel('Value [kN]'); ylabel('pdf');
subplot(132)
hist(d); title('Dead load'); xlabel('Value [kN]'); ylabel('pdf');
subplot(133)
hist(l); title('Live load'); xlabel('Value [kN]'); ylabel('pdf');

figure(2)
hist(g); title(['Safety margin - P_f = ',num2str(length(find(g<0))./n)]); xlabel('M'); ylabel('pdf');


%% How many extraction to have a given level of probability?
N = 1:1:10000;
Pf = zeros(length(N),1);
for i=1:length(N)

    % checks the effect from 1 siumlation to 10000
    r = exp(etaR+BetaR.*randn(N(i),1));
    d = muD + muD*COVD.*randn(N(i),1);
    l = u - log(-log(rand(N(i),1)))/a;
    
    g = r-d-l;
    
    Pf(i,1) = length(find(g<0))./N(i);
end
figure(3)
plot(Pf); title('Convergence of P_f'); xlabel('number of simulations'); ylabel('P_f');

figure(4)
tiledlayout(1,2)
nexttile
boxplot(storey_1_data)
ylabel('MPa')
ylim([0 50])
nexttile
boxplot(storey_2_data)
ylabel('MPa')
ylim([0 50])
