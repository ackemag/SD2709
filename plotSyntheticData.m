function plotSyntheticData()
% plotSyntheticData.m
% Function to plot synthetic ocean data (temperature, salinity, pressure) vs. depth.
% Author: Axel Magnusson   
% Date: 1/9/2025

% Define example depth, salinity, and temperature for synthetic data
depth = 0:2000;  % Depth from 0 to 2000 meters

% Ensure salinity and temperature have the same length as depth
salinity = 5:5:35;  % Salinity levels from 5 to 35 parts per thousand
temperature = atan((-depth + 500) / 100) * 7 + 11;  % Temperature profile as a function of depth
pressure = 1000 * 10 * depth;  % Pressure in dbar

% Plot parameters
linewidth = 2;
fontsize = 10;

% Create figure for the profiles
figure;
clf;

% -------------------------------------------------------------------------
% Subplot 1: Temperature vs. Depth
subplot(1, 4, 1);
plot(temperature, depth, 'LineWidth', linewidth);
xlabel('Water temp (°C)', 'FontSize', fontsize);
ylabel('Water depth (m)', 'FontSize', fontsize);
set(gca, 'YDir', 'reverse');  % Invert the y-axis to have positive depth down
grid on;
title('Temp', 'FontSize', fontsize);

% -------------------------------------------------------------------------
% Subplot 2: Pressure vs. Depth
subplot(1, 4, 2);
plot(pressure / 101325, depth, 'LineWidth', linewidth);  % Convert pressure to Bar
xlabel('Pressure (Bar)', 'FontSize', fontsize);
ylabel('Water depth (m)', 'FontSize', fontsize);
set(gca, 'YDir', 'reverse');  % Invert the y-axis to have positive depth down
grid on;
title('Pres', 'FontSize', fontsize);

% -------------------------------------------------------------------------
% Subplot 3: Density vs. Depth for different salinities
subplot(1, 4, 3);
hold on;
for sal = salinity
    [rho_w, ~, ~] = calcDensityWater(temperature, sal, pressure);
    plot(rho_w, depth, 'LineWidth', linewidth);
end
xlabel('Water density (kg/m^3)', 'FontSize', fontsize);
ylabel('Water depth (m)', 'FontSize', fontsize);
set(gca, 'YDir', 'reverse');  % Invert the y-axis to have positive depth down
grid on;
legend(arrayfun(@num2str, salinity, 'UniformOutput', false), 'Location', 'best', 'FontSize', fontsize);
title('Dens vs Sal', 'FontSize', fontsize);

% -------------------------------------------------------------------------
% Subplot 4: Density components vs. Depth
subplot(1, 4, 4);
[rho_w_pres, rho_w_temp, rho_w_sal] = calcDensityWater(temperature, salinity(end), pressure);
plot(rho_w_temp, depth, 'LineWidth', linewidth, 'DisplayName', 'Temp');
hold on;
plot(rho_w_sal, depth, 'LineWidth', linewidth, 'DisplayName', 'Salinity');
plot(rho_w_pres, depth, 'LineWidth', linewidth, 'DisplayName', 'Pressure');
xlabel('Water density (kg/m^3)', 'FontSize', fontsize);
ylabel('Water depth (m)', 'FontSize', fontsize);
set(gca, 'YDir', 'reverse');  % Invert the y-axis to have positive depth down
grid on;
legend('show', 'FontSize', fontsize);
title('Dens', 'FontSize', fontsize);

% Save the figure
saveas(gcf, 'SyntheticData_Profile.png');

end
