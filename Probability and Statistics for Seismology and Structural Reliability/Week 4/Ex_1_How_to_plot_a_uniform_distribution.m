%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

a = 3;     % Definition of the first extreme
b = 6;     % Definition of the second extreme


%% Plotting uniform distribution
figure(1)

% Subplot for the pdf
subplot(1,2,1)
% Coding of the interval of interest
dx = 0.1;
x = a/2:dx:b*1.3;   % This defines a vector with internaval dx

% Preallocation of the memory to speed up for loops
y = zeros(length(x),1); 

% Population of the y vector
for i=1:length(x)
    if x(i)<a
        y(i) = 0;
    elseif x(i) <b
        y(i) = 1/(b-a);
    else
        y(i) = 0;
    end
end

stairs(x,y,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
xlabel('x'); %it defines the label of the x axis
ylabel('f_x'); %it defines the label of the y axis
legend({'Uniform pdf'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)
%%

% Subplot for the CDF
subplot(1,2,2)

plot(x,cumsum(y).*dx,'k','linewidth',3); %it plots y defined on x, color black ('k') and linewidth 3
cumsum = cumsum(y).*dx;
xlabel('x'); %it defines the label of the x axis
ylabel('F_x'); %it defines the label of the y axis
legend({'Uniform CDF'},'location','northwest'); %it defines the item of the legend
grid on; %it applies the grid to the plot
axis square
set(gca,'fontsize',12)