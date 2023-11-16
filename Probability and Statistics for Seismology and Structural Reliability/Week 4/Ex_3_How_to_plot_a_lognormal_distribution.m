%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

eta = 3;        % Definition of the median
beta = 1;     % Definition of the logarithimc standard deviation

%% Plotting Normal distribution
figure(1)

% Subplot for the pdf
subplot(1,2,1)
% Coding of the interval of interest
dx = 0.1;
x = 0:dx:exp(log(eta)+3*beta);   % This defines a vector with internaval dx

% Preallocation of the memory to speed up for loops
y = lognpdf(x,log(eta),beta);

plot(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
xlabel('x'); %it defines the label of the x axis
ylabel('f_x'); %it defines the label of the y axis
legend({'LogNormal pdf'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)
xlim([0 40])
%%

% Subplot for the CDF
subplot(1,2,2)

y= normcdf(log(x),log(eta),beta);

plot(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
xlabel('x'); %it defines the label of the x axis
ylabel('F_x'); %it defines the label of the y axis
legend({'LogNormal CDF'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)
xlim([0 40])