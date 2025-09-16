function [drag_forces, Cd_total, f, Cf, Cd_body] = calc_drag_force(U, vehicle_params, rho_water, nu_water)
    % Calculate the drag force for a given AUV at different speeds and slenderness ratios.
    
    length = vehicle_params.length;
    diameter = vehicle_params.diameter;

    % Step 1: Calculate slenderness ratio (f = length/diameter)
    f = length / diameter;

    % Step 2: Calculate Reynolds number (Re = ?)
    % Use the formula for Reynolds number here
    Re = (U .* length) ./ nu_water;

    % Step 3: Calculate skin friction coefficient (Cf)
    % Use the ITTC57 correlation here
    Cf = 0.075 ./ ( (log10(Re) - 2).^2 );

    % Step 4: Calculate form factor (1+k) based on slenderness
    % Fill in form factor calculation here
    ratio = diameter / length;                 % d/l
    form_factor = 1 + 1.5*ratio^(3/2) + 7*ratio^3;   % (1 + k)

    % Step 5: Calculate total drag coefficient (Cd_total = ?)
    % Combine Cf with the form factor to get Cd_total
    % Cv = (1+k) * Cf   and   Cd_body = k * Cf
    Cd_total = form_factor .* Cf;          % = (1+k)Cf
    Cd_body  = (form_factor - 1) .* Cf;    % = k*Cf

    % Step 6: Calculate drag forces (drag_force = ?)
    % Use the formula F = 0.5 * rho * U^2 * A * Cd
    R   = diameter/2;
    Aw  = pi*diameter*(length - diameter) + 4*pi*R^2;  % wetted area [m^2]
    
    % Total viscous drag (vectorized over U)
    drag_forces = 0.5 .* rho_water .* (U.^2) .* Aw .* Cd_total;
end
