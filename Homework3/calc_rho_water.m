function [rho, rho_temp, rho_sal] = calc_rho_water(temp, salinity, pressure)
    % Calculate the density of seawater based on temperature, salinity, and pressure.

    % Step 1: Calculate temperature effect on density
    % Fill in the formula for rho_temp
   
    % %---------------------------------------------------------
    % Effect of temperature on water density [FIGURE 2.2]
    % (Millero and Poisson, 1981)
    %---------------------------------------------------------
    rho_temp = 999.842594 ...
         + 6.793952e-2*temp ...
         - 9.095290e-3*temp.^2 ...
         + 1.001685e-4*temp.^3 ...
         - 1.120083e-6*temp.^4 ...
         + 6.536336e-9*temp.^5;   % kg/m^3 (temp in °C)


    % Step 2: Calculate salinity effect on density
    % Use the equations provided in the script for salinity effects
    
    %---------------------------------------------------------
    % Effects of salinity [FIGURE 2.3]
    % (International One Atmosphere Equation)
    %---------------------------------------------------------
    A =  8.24493e-1 ...
       - 4.0899e-3*temp ...
       + 7.6438e-5*temp.^2 ...
       - 8.2467e-7*temp.^3 ...
       + 5.3875e-9*temp.^4;
    
    B = -5.72466e-3 ...
       + 1.0227e-4*temp ...
       - 1.6546e-6*temp.^2;
    
    C = 4.8314e-4;
    
    rho_sal = A.*salinity + B.*(salinity.^1.5) + C.*(salinity.^2);  % kg/m^3
    rho_0   = rho_temp + rho_sal;   % 1-atm density; Step 3 will add pressure

    
    % Step 3: Adjust density for pressure
    % Use the bulk elastic modulus (EB) to calculate the final density

    EB = 2.151e9;      % Pa  (≈ 2.151 GPa per lecture notes)
    
    rho = rho_0 ./ (1 - pressure./EB);   % final density at pressure


end
