%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

P = [10 42]; % Longitude and Latitude of the point

%% How to enter the folder Italy
cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..

plot([Italy_boundaries.X],[Italy_boundaries.Y],'k')
axis equal
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on

%% Is the point inside the polygon?
index = inpolygon(P(:,1),P(:,2),[Italy_boundaries.X],[Italy_boundaries.Y]);

if index == 1
    disp('YES')
    hold on
    plot(P(:,1),P(:,2),'kp','linewidth',1,'markerfacecolor','y','markersize',12)
else
    disp('NO')
    hold on
     plot(P(:,1),P(:,2),'ks','linewidth',1,'markerfacecolor','r','markersize',12)
end