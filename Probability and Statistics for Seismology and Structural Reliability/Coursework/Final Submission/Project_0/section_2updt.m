clear
clc
close all

%% NUMBER 2.4
% m stands for mean and s for standard dev

% General Parameters
k = 100;
mu_F = 18.93;
sigma_F = 5.68;
h = 7.35;
mu_B = 256.37;
sigma_B = 17.41;

% Column 1
m1_demand = (mu_F*h)/3;
s1_demand = (sigma_F*h)/3;
m1_cap = 0.5 * mu_B;
s1_cap = 0.5 * sigma_B;
z1 = (m1_cap-m1_demand)/ sqrt(((s1_demand^2)+(s1_cap^2)));
B_1 = -z1;
Pf_1 = normcdf(B_1);

% Column 2
m2_demand = (mu_F*h)/2;
s2_demand = (sigma_F*h)/2;
m2_cap = 0.5 * mu_B;
s2_cap = 0.5 * sigma_B;
z2 = (m2_cap-m2_demand)/ sqrt(((s2_demand^2)+(s2_cap^2)));
B_2 = -z2;
Pf_2 = normcdf(B_2);

% Column 3
m3_demand = (mu_F*h)/6;
s3_demand = (sigma_F*h)/6;
m3_cap = 0.2 * mu_B;
s3_cap = 0.2 * sigma_B;
z3 = (m3_cap-m3_demand)/ sqrt(((s3_demand^2)+(s3_cap^2)));
B_3 = -z3;
Pf_3 = normcdf(B_3);

% Probability of failure (parallel system)
Pf = Pf_1*Pf_2*Pf_3;