%% INPUT
clear
clc

% Magnitude = 4;  % min magnitude
% 
% Initial_date = '1900-01-01 00:00:00'; %YYYY-MM-DD
% Final_date   = '2020-08-01 00:00:00'; %YYYY-MM-DD
% 
% start_year  = '1964';
% end_year    = '2019';
% start_month = '1';
% end_month   = '12';
% 
% start_day='1';
% start_time='00:00:00';
% end_day='31';
% end_time='00:00:00';
% 
% minlat=36.60;
% maxlat=39.71;
% minlon=14.32;
% maxlon=16.65;

Magnitude = 4;  % min magnitude

Initial_date = '1900-01-01 00:00:00'; %YYYY-MM-DD
Final_date   = '2020-10-01 00:00:00'; %YYYY-MM-DD

start_year  = '1900';
end_year    = '2020';
start_month = '1';
end_month   = '9';

start_day='1';
start_time='00:00:00';
end_day='31';
end_time='00:00:00';

minlat=45.27;
maxlat=46.90;
minlon=9.46;
maxlon=13.50;
%% FIRST CATALOG
%%% all the infor are available at https://earthquake.usgs.gov/fdsnws/event/1/
disp('%%%% USGS CATALOG %%%%')
Rectangle_of_interest = [minlat maxlat minlon maxlon];
[year, month, day, hour, minute, sec, lat_e, long_e, depth, mag, magType] = LoadComCat(datenum(Initial_date),datenum(Final_date),Magnitude,Rectangle_of_interest);

for ii=1:length(year)
CATALOG(ii).text1 = 'USGS';
CATALOG(ii).text2 = 'USGS';
CATALOG(ii).YEAR = year(ii);
CATALOG(ii).MONTH = month(ii);
CATALOG(ii).DAY = day(ii);
CATALOG(ii).HOUR = hour(ii);
CATALOG(ii).MIN = minute(ii);
CATALOG(ii).SEC = sec(ii);
CATALOG(ii).LATITUDE = lat_e(ii);
CATALOG(ii).LONGITUDE = long_e(ii);
CATALOG(ii).DEPTH = depth(ii)/1000;
CATALOG(ii).MW = magnitude_converter(mag(ii),magType(ii));
end

tmp_length = length(year);
disp('%%%% ISC CATALOG %%%%')
URL_ISC_final

save CATALOG CATALOG

%%
%%


%% Import the shape file of the source zones
cd ..
cd('SHP')
Zones = shaperead('9_Reggio_zones.shp');
cd ..
cd('EX_4_SINGLE_CATALOG')

% Select only the area where there are more earthquakes
Area_to_keep = [2 5 8];

Final_zones = Zones(Area_to_keep);

% Final figure
figure
cd ..
cd('Italy')
Italy_boundaries = shaperead('ITA_adm0.shp');
cd ..
cd('EX_4_SINGLE_CATALOG')

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

