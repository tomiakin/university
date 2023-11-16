%% INPUT
% clear
clc

Magnitude = 3;  % min magnitude

Initial_date = '1900-01-01 00:00:00'; %YYYY-MM-DD
Final_date   = '2019-12-31 00:00:00'; %YYYY-MM-DD

start_year  = '1900';
end_year    = '2019';
start_month = '1';
end_month   = '12';

start_day='1';
start_time='00:00:00';
end_day='31';
end_time='00:00:00';

minlat=36.5;
maxlat=40;
minlon=13;
maxlon=18;


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
CATALOG(ii).DEPTH = depth(ii);
CATALOG(ii).MW = magnitude_converter(mag(ii),magType(ii));
end

tmp_length = length(year);
disp('%%%% ISC CATALOG %%%%')
URL_ISC_final

save CATALOG CATALOG