clear
clc
close all

%% Load the catalog
load SHAREver3
figure(1)
hold on

% Import the shape file of the source zones
cd('SHP')
Zones = shaperead('seismic_zones.shp');
cd ..

% Plot all seismic zones in black
plot([Zones.X], [Zones.Y], 'k', 'linewidth', 2)

% Highlight points within selected zones in yellow
highlight_zones = [1 5];
for i = 1:length(Zones)
    in_highlighted_zone = inpolygon(SHAREver3.Lon, SHAREver3.Lat, Zones(i).X, Zones(i).Y);
    
    if ismember(i, highlight_zones)
        plot(SHAREver3.Lon(in_highlighted_zone), SHAREver3.Lat(in_highlighted_zone), 'g.', 'markersize', 5);

    else
        plot(SHAREver3.Lon(in_highlighted_zone), SHAREver3.Lat(in_highlighted_zone), 'k.', 'markersize', 5);
    end
end

% Plot earthquakes


% POINT OF INTEREST
X = [35.36 25.12];
plot(X(:,2), X(:,1), 'ko', 'markerfacecolor', 'y', 'markersize', 9, 'linewidth', 2)

% Plot the source zones in red
plot([Zones.X], [Zones.Y], 'r', 'linewidth', 2)

axis equal
grid on
% plot(SHAREver3.Lon, SHAREver3.Lat, 'k.')
% Display text ID at the center of the seismic zones
for i = 1:length(Zones)
    text(Zones(i).LONC, Zones(i).LATC, num2str(i), 'fontweight', 'bold')
end

box on
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = 11;
ylabel('Latitude [deg]', 'Interpreter', 'latex');
xlabel('Longitude [deg]', 'Interpreter', 'latex');

% Set limits
xlim([0.9 * min([Zones.X]) 1.1 * max([Zones.X])])
ylim([0.99 * min([Zones.Y]) 1.01 * max([Zones.Y])])
