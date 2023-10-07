%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)


%% Plotting example 1
figure(1)
x = 0:0.1:10;   % This defines a vector from 0 to 10 with internaval 0.1
y = sin(pi.*x); % This is a function of x

plot(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
xlabel('x'); %it defines the label of the x axis
ylabel('y'); %it defines the label of the y axis
legend({'Sine'}); %it defines the item of the legend
grid on; %it applies the grid to the plot
xlim([0 20]); %it defines the x limits
ylim([-2 2]); %it defines the y limits

%% Plotting example 2
figure(2)
x = 0:0.1:10;   % This defines a vector from 0 to 10 with internaval 0.1
y1 = sin(pi.*x); % This is a function of x
y2 = cos(pi.*x); % This is a function of x

plot(x,y1,'k','linewidth',3)
hold on; % this command is ESSENTIAL when you have more the one plot
plot(x,y2,':r','linewidth',3); % ':r' means that the plot will be red and dotted line
xlabel('x')
ylabel('y')
legend({'Sine','Cosine'}); %this is when you have more than one series
grid on
xlim([0 20])
ylim([-2 2])

%% Plotting example 3 - scatter plot
figure; % a new figure is created

data = load('MoEdata.txt');
scatter(data(:,1),data(:,2),'bo','linewidth',2,'markerfacecolor','y');
box on
grid on

xlabel('variable 1');
ylabel('variable 2');
axis square
set(gca,'FontSize',16);

%% Plotting example 4 - bar
figure;
data = load('concretestrength.txt');
mean(data)
median(data)
std(data)
skewness(data)
kurtosis(data)
count = hist(data,18:1:34);

subplot(1,2,1)
bar(18:34,count/40)
xlabel('f_{concrete} [MPa]');
ylabel('frequency');
axis square
set(gca,'FontSize',16);
legend({'Frequency'})

countsum = cumsum(count/40);
subplot(1,2,2)
bar(18:34,countsum,'FaceColor',[1 1 0])
xlabel('f_{concrete} [MPa]');
ylabel('CDF');
axis square
set(gca,'FontSize',16);
legend({'Cumulative'})

%% Plotting example 5 - Empricial Cumulative Distribution Function ecdf
figure;
data = load('concretestrength.txt');
[frequency,values]=ecdf(data);
stairs(values,frequency,'linewidth',2)
xlabel('f_{concrete} [MPa]');
ylabel('ECDF');
axis square
set(gca,'FontSize',16);
legend({'Empirical CDF'},'location','northwest')
grid on