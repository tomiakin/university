%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

%% Import the shape file of the source zones
cd('SHP')
Zones = shaperead('9_Reggio_zones.shp');
cd ..

%% Select only the area where there are more earthquakes
Area_to_keep = [1 5 8];

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

minlat = round(min([Final_zones.Y])*100)/100
maxlat = round(max([Final_zones.Y])*100)/100
minlon = round(min([Final_zones.X])*100)/100
maxlon = round(max([Final_zones.X])*100)/100
%% Load the catalog
load SHAREver3
for i=1:length(Final_zones)
    
    index = inpolygon(SHAREver3.Lon,SHAREver3.Lat,Final_zones(i).X,Final_zones(i).Y);
    SUB_CATALOG(i).CATALOG = SHAREver3(index,:);
    hold all
    plot(SUB_CATALOG(i).CATALOG.Lon,SUB_CATALOG(i).CATALOG.Lat,'+')
    % scatter(SHAREver3.Lon,SHAREver3.Lat,0.002*10^-15*10.^(1.5*SHAREver3.Mw+9.05),SHAREver3.Mw,'fill')
    % colorbar
end

save SUB_CATALOG SUB_CATALOG % same as Ex4 but it saves!?
