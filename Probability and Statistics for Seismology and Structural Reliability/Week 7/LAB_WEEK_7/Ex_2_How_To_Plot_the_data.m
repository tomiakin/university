%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

%% How to enter the folder Italy
cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..

plot([Italy_boundaries.X],[Italy_boundaries.Y],'k')
axis equal
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on

%% POINT OF INTEREST

X = [38.111389	15.661944]; % here latitude is 38.11? 

hold on
plot(X(:,2),X(:,1),'ko','markerfacecolor','y','markersize',9,'linewidth',2)

%% Import the shape file of the source zones
cd('SHP')
Zones = shaperead('9_Reggio_zones.shp');
cd ..

hold on
plot([Zones.X],[Zones.Y],'r','linewidth',2)

%% How to make appear a text ID at the centre of the seismic zones
for i=1:length(Zones)
    hold on
    text(Zones(i).LONC,Zones(i).LATC,num2str(i),'fontweight','bold')
end

axis equal