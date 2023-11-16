% THIS IS THE CODE FROM THE EXAMPLE ON WHAT IS MONTE CARLO.. N = 10000, not
% pi will change on each run
Radius = 10; %[m]
Simulations = 10000;
%%
x1=-Radius:0.01:Radius;
x2=Radius:-0.01:-Radius;
y1=sqrt(Radius^2-x1.^2);
y2=-y1;

x=[x1,x2];
y=[y1,y2];

plot(x,y,'r','linewidth',3)
hold on
plot(0,0,'ko','markerfacecolor','y','markersize',14)

axis equal

xlim(1.2*[-Radius Radius])
ylim(1.2*[-Radius Radius])

box off
grid on

xlabel('x - [m]')
ylabel('y - [m]')
%%

hold on
plot([-Radius Radius Radius -Radius -Radius],[Radius Radius -Radius -Radius Radius],'k','linewidth',2) 


%%
xc=-Radius+2*Radius*rand(Simulations,1);
yc=-Radius+2*Radius*rand(Simulations,1);

hold on
plot(xc,yc,'k.')

%%
index = inpolygon(xc,yc,x,y);
hold on
plot(xc(index),yc(index),'g.')
tmp = xc(index);

title(['Area = ',num2str(length(tmp)/Simulations*(2*Radius)^2),' \pi = ',num2str(4*length(tmp)/Simulations*(2*Radius)^2/(2*Radius)^2)])
