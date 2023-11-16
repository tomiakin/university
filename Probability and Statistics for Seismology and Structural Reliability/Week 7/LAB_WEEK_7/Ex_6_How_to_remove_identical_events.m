clear
clc

%% Thi script is an example on how to remove identical events from a catalog.
% When you put many catalogs together from different resources (e.g., USGS,
% ISC, etc.) you may have the same earthquake listed in the catalog twice 
% or more times.

% The function to remove events is called Rep_cleaner.m

% This function returns the indexes of the rows of the original catalog
% that correspond to a single event

% The function can work on any subset of the original catalog (in terms of
% columns) to identify identical earthquake events 


%% EXAMPLE: Removal of Identical Events fro a catalog
cd EX_5_SINGLE_CATALOG
load CATALOG
cd ..

% In general, to identify identical elements in a catalog, it is sufficient
% to identify few characteristics: year, month, day, longitude, latitude,
% magnitude.

Year  = [CATALOG.YEAR];  
Month = [CATALOG.MONTH];  
Day   = [CATALOG.DAY];  

Long  = [CATALOG.LONGITUDE];  
Lat   = [CATALOG.LATITUDE];  
Mag   = [CATALOG.MW]; 

% I am now assembling the different columns extracted in the original
% catalog in a temporary support matrix.

temporary_catalog = [Year',Month',Day',round(Long'*100)/100,round(Lat'*100)/100,round(Mag'*10)/10];

% I am now using the function Rep_cleaner.m to obtain the indexes of the
% events witout repetition.

indexes = Rep_cleaner(temporary_catalog);

% Finally, I am going to clean the catalog of duplicates.

CATALOG_CLEAN = CATALOG(indexes);

%%
%%
%% Import the shape file of the source zones
cd('SHP')
Zones = shaperead('9_Reggio_zones.shp');
cd ..

% Select only the area where there are more earthquakes
Area_to_keep = [2 5 8];

Final_zones = Zones(Area_to_keep);

% Final figure
figure
cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..

plot([Italy_boundaries.X],[Italy_boundaries.Y],'k')
axis equal
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on

% POINT OF INTEREST

X = [38.111389	15.661944];

hold on
plot(X(:,2),X(:,1),'ko','markerfacecolor','y','markersize',9,'linewidth',2)

% Import the shape file of the source zones

hold on
plot([Final_zones.X],[Final_zones.Y],'r','linewidth',2)

% How to make appear a text ID at the centre of the seismic zones
for i=1:length(Final_zones)
    hold on
    text(Final_zones(i).LONC,Final_zones(i).LATC,num2str(i),'fontweight','bold')
end

hold on
plot([CATALOG.LONGITUDE],[CATALOG.LATITUDE],'.')
axis equal

% 
xlim([0.9*min([Final_zones.X]) 1.1*max([Final_zones.X])])
ylim([0.99*min([Final_zones.Y]) 1.01*max([Final_zones.Y])])

% Load the catalog