%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

mu = 3;        % Definition of the mean
sigma = 1;     % Definition of the standard deviation

%% Plotting Normal distribution
figure(1)

% Subplot for the pdf
subplot(1,2,1)
% Coding of the interval of interest
dx = 0.1;
x = mu-3*sigma:dx:mu+3*sigma;   % This defines a vector with internaval dx

% Preallocation of the memory to speed up for loops
y = normpdf(x,mu,sigma);

plot(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
xlabel('x'); %it defines the label of the x axis
ylabel('f_x'); %it defines the label of the y axis
legend({'Normal pdf'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)
%%

% Subplot for the CDF
subplot(1,2,2)

y= normcdf(x,mu,sigma);

plot(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3

xlabel('x'); %it defines the label of the x axis
ylabel('F_x'); %it defines the label of the y axis
legend({'Normal CDF'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)