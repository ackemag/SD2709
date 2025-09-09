% Homework 2 example script
% SD2709 Underwater Technology
% Homework: Module 2
% Sensors and Underwater Acoustics
% Author: [Axel Magnusson]
% Date: [2025-09-09]

% MATLAB version of the Python script
% Constants (Task 2 values)
v = 1500;           % Speed of sound in water in m/s
f = 37.5e3;         % Frequency in Hz
lambda = v / f;     % Wavelength in meters
k = 2 * pi / lambda; % Wave number

% Array design for Task 2
d = lambda/2;       % Element spacing
N = 24;             % Number of elements
D = (N-1)*d;        % Aperture length

% Angle range for beam pattern
theta = linspace(-pi/2, pi/2, 2001);  % Angle from -pi/2 to pi/2


% Beam pattern (uniform linear array, broadside)
psi = k * d * sin(theta);
beam_pattern = ones(size(theta));
idx = abs(psi) > 1e-12;
beam_pattern(idx) = sin(N*psi(idx)/2) ./ (N*sin(psi(idx)/2));

% Normalizing the beam pattern and converting to dB
beam_pattern_normalized = abs(beam_pattern) / max(abs(beam_pattern));
normalized_power_db = 20 * log10(max(beam_pattern_normalized, 1e-6)); % floor avoids -Inf
peak_power_db = max(normalized_power_db);
half_power_db = peak_power_db - 3;

% Plotting the beam pattern
figure('Position', [100, 100, 1000, 600]);
plot(rad2deg(theta), normalized_power_db, 'LineWidth', 2);
title(sprintf('Uniform Linear Array  f = %.1f kHz, N = %d, d = %.3f m', f/1e3, N, d));
xlabel('Angle (degrees)');
ylabel('Normalized Magnitude (dB)');
hold on;
yline(peak_power_db, 'r--', 'Peak Power', 'LabelVerticalAlignment', 'bottom');
yline(half_power_db, 'g--', 'Half Power (-3 dB)', 'LabelVerticalAlignment', 'bottom');

grid on;
ylim([-40, 1]);   % Display down to -40 dB
xlim([-90, 90]);  % Limit the x-axis from -90 to 90 degrees
hold off;
