% calc_rho_water.m

function [rho, rho_temp, rho_sal] = calc_rho_water(temp, salinity, pressure)
    % Calculation of water density.
    % temp in [Celsius]
    % salinity [parts per thousand]
    % pressure [Pa] (absolute)
    % Returns:
    %    rho [kg/m^3]: Water density considering temperature, salinity, and pressure
    %    rho_temp [kg/m^3]: Water density considering temperature only
    %    rho_sal [kg/m^3]: Water density considering temperature and salinity

    % ---------------------------------------------------------
    % Effect of temperature (Millero and Poisson, 1981)
    % ---------------------------------------------------------
    rho = (999.842594 + 6.793952e-2 .* temp - 9.095290e-3 .* temp.^2 + ...
           1.001685e-4 .* temp.^3 - 1.120083e-6 .* temp.^4 + ...
           6.536336e-9 .* temp.^5);
    rho_temp = rho;

    % ---------------------------------------------------------
    % Effects of salinity
    % (International One Atmosphere Equation)
    % ---------------------------------------------------------
    A = (8.24493e-1 - 4.0899e-3 .* temp + 7.6438e-5 .* temp.^2 - ...
         8.2467e-7 .* temp.^3 + 5.3875e-9 .* temp.^4);
    B = -5.72466e-3 + 1.0227e-4 .* temp - 1.6546e-6 .* temp.^2;
    C = 4.8314e-4;
    rho = rho + A .* salinity + B .* (salinity.^1.5) + C .* (salinity.^2);
    rho_sal = rho;

    % ---------------------------------------------------------
    % Effects of hydrostatic pressure
    % ---------------------------------------------------------
    EB = 2.151e9;  % [Pa] Typical bulk elastic modulus
    rho = rho ./ (1 - pressure ./ EB);

end
