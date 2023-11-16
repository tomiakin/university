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

X = [38.111389	15.661944];

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

%% 
xlim([0.9*min([Zones.X]) 1.1*max([Zones.X])])
ylim([0.99*min([Zones.Y]) 1.01*max([Zones.Y])])

%% Load the catalog
load SHAREver3
hold on
 plot(SHAREver3.Lon,SHAREver3.Lat,'k.')
%  scatter(SHAREver3.Lon,SHAREver3.Lat,0.002*10^-15*10.^(1.5*SHAREver3.Mw+9.05),SHAREver3.Mw,'fill')
%  colorbar

%% Select only the area where there are more earthquakes
Area_to_keep = [2 5 8];

Final_zones = Zones(Area_to_keep);

%% Final figure
figure
cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..

plot([Italy_boundaries.X],[Italy_boundaries.Y],'k')
axis equal
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on

%% POINT OF INTEREST

X = [38.111389	15.661944];

hold on
plot(X(:,2),X(:,1),'ko','markerfacecolor','y','markersize',9,'linewidth',2)

%% Import the shape file of the source zones

hold on
plot([Final_zones.X],[Final_zones.Y],'r','linewidth',2)

%% How to make appear a text ID at the centre of the seismic zones
for i=1:length(Final_zones)
    hold on
    text(Final_zones(i).LONC,Final_zones(i).LATC,num2str(i),'fontweight','bold')
end

axis equal

%% 
xlim([0.9*min([Final_zones.X]) 1.1*max([Final_zones.X])])
ylim([0.99*min([Final_zones.Y]) 1.01*max([Final_zones.Y])])

%% Load the catalog
load SHAREver3
hold on
plot(SHAREver3.Lon,SHAREver3.Lat,'k.')
% scatter(SHAREver3.Lon,SHAREver3.Lat,0.002*10^-15*10.^(1.5*SHAREver3.Mw+9.05),SHAREver3.Mw,'fill')
% colorbar
