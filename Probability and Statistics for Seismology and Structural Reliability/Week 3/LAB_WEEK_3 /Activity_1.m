clear
clc
close all

%% ACTIVITY 1
% ACquisition of the data
Traffic_gaps=[1.27, 3.65, 3.36, 1.70, 6.90, 31.88, 14.04, 9.78, 1.21, 0.43, 13.08, 0.42, 2.05, 3.55,...
5.20, 4.79, 6.28, 17.60, 1.57, 9.80, 3.05, 2.58, 9.41, 6.79, 0.42, 14.95, 0.71, 6.13,...
4.83, 3.78, 0.36, 17.11, 4.69, 21.91, 10.56, 7.71, 0.98, 2.17, 3.48, 10.13, 3.41,...
14.32, 6.85, 21.57, 1.37, 8.73, 24.87, 12.27, 14.88, 14.08];

% Sorting of the data in ascending order
Sorted_data = sort(Traffic_gaps);

% Quantification of the number of data N
N = length(Traffic_gaps);

% Creation of the SEED vector: 1,2,...,N
SEED = 1:1:N;

% Calculation of the percentiles
Percentiles = SEED./(1+N);

%% Plot of the percentiles against the sorted data
figure % creation of a new figure

subplot(1,2,1) % definition of a panel of plots (1 x 2) - first plot
plot(Sorted_data,Percentiles,'ko','markersize',10,'linewidth',2,'markerfacecolor','y')
xlabel('Sorted data')
ylabel('Percentiles')
axis square
set(gca,'fontsize',16)
hold on
grid on

subplot(1,2,2) % definition of a panel of plots (1 x 2) - second plot
streched_percentiles = -log(1-Percentiles);
plot(Sorted_data,streched_percentiles,'ks','markersize',10,'linewidth',2,'markerfacecolor','g')
xlabel('Sorted data')
ylabel('Percentiles')
axis square
set(gca,'fontsize',16)
hold on
grid on

ft1 = fittype({'x'});
p1 = fit(Sorted_data',streched_percentiles',ft1); %This creates a 'cfit' variable p that is your fitted function
x_fit = linspace(0,max(Sorted_data),100); %x-values to evaluate
y1_fitted = feval(p1, x_fit); %y-values for the evaluated x-values
hold on
plot(x_fit,y1_fitted,'r--','linewidth',2)
legend({'Data';['\lambda = ',num2str(p1.a)]},'location','northwest')

