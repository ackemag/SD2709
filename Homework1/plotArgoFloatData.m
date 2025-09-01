% ArgoDataPlot_Map_Validated.m
% SD2709 Underwater Technology
% Plot Latitude vs. Longitude on a Map, along with Temperature, Salinity, and Pressure Profiles from Argo Float Data
% Author: Axel Magnusson
% Date: 1/9/2025

function plotArgoFloatData(filename)

% Load the Argo data from a CSV file with the original column names preserved
%filename = 'PR_PF_4903884.csv';  % Replace with your Argo CSV file name
data = readtable(filename, 'VariableNamingRule', 'preserve');

% Automatically detect and assign columns based on partial name matching
columnNames = data.Properties.VariableNames;
depth_col = columnNames(contains(columnNames, 'PRES', 'IgnoreCase', true));
temp_col = columnNames(contains(columnNames, 'TEMP', 'IgnoreCase', true));
salinity_col = columnNames(contains(columnNames, 'PSAL', 'IgnoreCase', true));
lat_col = columnNames(contains(columnNames, 'LATITUDE', 'IgnoreCase', true));
lon_col = columnNames(contains(columnNames, 'LONGITUDE', 'IgnoreCase', true));

% Extract the data using the detected column names
depth = data.(depth_col{1});
temperature = data.(temp_col{1});
salinity = data.(salinity_col{1});
latitude = data.(lat_col{1});
longitude = data.(lon_col{1});

% Filter out invalid latitude and longitude values
validIndices = latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180;
latitude = latitude(validIndices);
longitude = longitude(validIndices);

% Plot parameters
markersize  = 12;
fontsize    = 14;
linewidth   = 2;

% Create a figure for the map
figure(1);
clf;

% Plot Latitude vs. Longitude on a map
geoplot(latitude, longitude, '-x', 'MarkerSize', markersize, 'LineWidth',linewidth);
geobasemap('colorterrain');  % Choose a base map (e.g., 'colorterrain', 'satellite', etc.)
title('Latitude vs. Longitude on a Map', 'FontSize', fontsize);

% Save the map plot
saveas(gcf, 'Argo_Map_LatLon.png');

% Create a figure for the profiles
figure(2);
clf;

% Subplot 1: Temperature vs. Depth
subplot(1, 3, 1);
plot(temperature, depth, 'o', 'MarkerSize', markersize);
set(gca, 'YDir', 'reverse');
xlabel('Temperature (°C)', 'FontSize', fontsize);
ylabel('Depth (m)', 'FontSize', fontsize);
title('Temperature vs. Depth', 'FontSize', fontsize);
grid on;

% Subplot 2: Salinity vs. Depth
subplot(1, 3, 2);
plot(salinity, depth, 'o', 'MarkerSize', markersize);
set(gca, 'YDir', 'reverse');
xlabel('Salinity (PSU)', 'FontSize', fontsize);
title('Salinity vs. Depth', 'FontSize', fontsize);
grid on;

% Subplot 3: Pressure vs. Depth
subplot(1, 3, 3);
plot(depth, depth, 'o', 'MarkerSize', markersize);
set(gca, 'YDir', 'reverse');
xlabel('Pressure (dbar)', 'FontSize', fontsize);
title('Pressure vs. Depth', 'FontSize', fontsize);
grid on;

% Save the profiles plot
saveas(gcf, 'Argo_Profile_Temp_Sal_Pressure.png');

disp('Argo plots saved as "Argo_Map_LatLon.png" and "Argo_Profile_Temp_Sal_Pressure.png".');
