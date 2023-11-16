%% Solution of the integral numerically
clear
clc
a = 7.9;
s1 = 300:1:2700;
r1 = 100:1:600;

s = repmat(s1,length(r1),1);
r = repmat(r1',1,length(s1));
Pf=zeros(size(s));


tmp = normpdf(r,350,35).* normpdf(s,1500,300).*1*1;
tmp(r.*a-s>=0)=0;
Pf = sum(sum(Pf + tmp));

disp('***************************************************')
disp('Numerical integration of the Pf integral')
disp([' The probability of failure is: ', num2str(Pf)])
disp([' The reliability index \beta is: ', num2str(-norminv(Pf))])
disp('***************************************************')
%% Solution of the integral with the FOSM (First Orfer Second Moment)
beta = (350*a-1500)/sqrt((35^2)*(a^2)+(300^2));
Pf = normcdf(-beta);

disp('***************************************************')
disp('Solution of the integral with the FOSM (First Orfer Second Moment)')
disp([' The probability of failure is: ', num2str(Pf)])
disp([' The reliability index \beta is: ', num2str(beta)])
disp('***************************************************')
%% Solution of the integral via Monte Carlo simulation
Nsimulation = 100000;
rv = 350+35*randn(Nsimulation,1);
sv = 1500+300*randn(Nsimulation,1);
g=rv.*a-sv;
Pf = length(find(g<0))/length(g);
beta = -norminv(Pf);

disp('***************************************************')
disp('Solution of the integral via Monte Carlo simulation')
disp([' The probability of failure is: ', num2str(Pf)])
disp([' The reliability index \beta is: ', num2str(beta)])
disp('***************************************************')