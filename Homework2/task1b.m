% Task 1b: Rescue pinger range calculation
% Switch DI from 0 to 13.8 for Task2 2
SL = 160.5;      % Source Level [dB re 1 uPa @ 1m]
NL = 41.5;       % Noise Level [dB]
DI = 0;          % Directivity Index [dB]
DT = 12;         % Detection Threshold [dB]

alpha = 10;      % Absorption coefficient [dB/km]
TLmax = SL - (NL - DI) - DT;   % Maximum allowed transmission loss

% Define equation for transmission loss
TLfun = @(r) 20*log10(r) + alpha*(r/1000);  % r in meters

% Solve numerically (fzero finds root of TLfun(r) - TLmax = 0)
range_guess = 3000;                         % initial guess [m]
r_solution = fzero(@(r) TLfun(r) - TLmax, range_guess);

fprintf('Maximum detection distance ≈ %.2f km\n', r_solution/1000);

