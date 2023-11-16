%% ACTIVITY 
Data = [0.35, 0.40, 0.41, 0.42, 0.43, 0.48, 0.49, 0.58, 0.68, 0.70, 0.75, 0.87, 0.96];

% Sorting of the data in ascending order
Sorted_data = sort(Data);

% Quantification of the number of data N
N = length(Data);

% Creation of the SEED vector: 1,2,...,N
SEED = 1:1:N;

% Calculation of the percentiles
Percentiles = SEED./(N+1);

%%
subplot(5,1,1:3)
plot(log(Sorted_data),Percentiles,'ko','markersize',10,'linewidth',2,'markerfacecolor','y')
xlabel('log(Sorted data)')
ylabel('Percentiles')
set(gca,'fontsize',16)
hold on
grid on

%% Fit a model
p1 = polyfit(log(Sorted_data)',Percentiles',1); % fit a line to the data
x_fit = linspace(min(log(Sorted_data)),max(log(Sorted_data)),100); %x-values to evaluate
y1_fitted = polyval(p1, x_fit); %y-values for the evaluated x-values
hold on
plot(x_fit,y1_fitted,'r--','linewidth',2)
legend({'Data';'Fit'},'location','northwest')

%% Plot of the errors
subplot(5,1,5)
errors = Percentiles - polyval(p1, log(Sorted_data));
stem(log(Sorted_data),abs(errors),'k','linewidth',2)
xlabel('log(Sorted data)')
ylabel('Errors')
set(gca,'fontsize',16)
hold on
grid on

% Identification of the maximum and overlay with red color
index = find(errors==max(errors));
hold on
plot(log(Sorted_data(index)),errors(index),'rs','linewidth',2)


% Using kstest
[h,p,ks2stat,cv] = kstest((log(Data)-mean(log(Data)))/std(log(Data)),'Alpha',0.02)

%In kstest2, the decision to reject the null hypothesis is based on comparing the p-value p with the significance level Alpha, not by comparing the test statistic ks2stat with a critical value.
x1 = Percentiles;
x2 = polyval(p1, log(Sorted_data));
[h,p,ks2stat] = kstest2(x1,x2,'Alpha',0.02)

