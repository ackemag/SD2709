% Flow regime calculator with critical Reynolds number
clear; clc;
% Switch parameters for other vehicles (currently for SCOUT)
% Parameters
U = linspace(0, 4*0.514444, 50);   % speeds (m/s)
L = 1.3;                       % characteristic length (m)
nu = 1.00e-6;                 % kinematic viscosity (m^2/s)

% Calculate Reynolds number
Re = (U .* L) ./ nu;

% Flow regimes (based on hydromechanics notes)
regime = strings(size(U));
for i = 1:length(U)
    if Re(i) < 2e5
        regime(i) = "Laminar";
    elseif Re(i) < 5e5
        regime(i) = "Transition";
    else
        regime(i) = "Turbulent";
    end
end

% Display table
T = table(U', Re', regime', 'VariableNames', {'Speed_m_per_s','Reynolds','Regime'});
disp(T)

% Plot
figure('Color','b');
plot(U, Re, 'b-', 'LineWidth', 1.5); hold on;
yline(2e5, 'r--', 'Re = 2e5 (Laminar→Transition)');
yline(5e5, 'g--', 'Re = 5e5 (Transition→Turbulent)');
xlabel('Speed U (m/s)');
ylabel('Reynolds number');
title('Flow Regime vs Speed');
grid on;
