%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

load SHAREver3

cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..

plot([Italy_boundaries.X],[Italy_boundaries.Y],'r','linewidth',3)
axis equal
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on
%% How to plot the spatial distribution of the earhtquakes
hold on
scatter(SHAREver3.Lon,SHAREver3.Lat,0.005*10^-15*10.^(1.5*SHAREver3.Mw+9.05),SHAREver3.Mw,'fill')
xlim([5.22 19.93])
ylim([35.49 47.09])