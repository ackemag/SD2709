function part_A1(vehicle_params)
    % Part A1: Buoyancy calculation for REMUS AUV at different depths.

    % Step 1: Calculate the volume of the AUV using the block coefficient formula
    L  = vehicle_params.length;      % [m]
    D  = vehicle_params.diameter;    % [m]
    CB = vehicle_params.CB;          % [-]
    R  = D/2;                        % [m]
    
    % Geometric "ideal" spindle volume: cylinder (L-D) + two hemispheres
    V_geom = pi*R^2*(L - D) + (4/3)*pi*R^3;   % [m^3]
    
    % Apply block coefficient
    V = CB * V_geom;                            % displaced volume [m^3]
   
    
    % Step 2: Create depth, temperature, and salinity arrays for analysis
    switch upper(vehicle_params.name)
        case 'SCOUT'
            z_max = 30;    % Gothenburg harbor
        case 'REMUS'
            z_max = 200;   % Part A
        otherwise
            % Fallback
            if isfield(vehicle_params,'max_depth') && ~isempty(vehicle_params.max_depth)
                z_max = vehicle_params.max_depth;
            else
                z_max = 200;
            end
    end

    dz = 1;                 % depth step [m]
    z  = (0:dz:z_max).';    % column vector

    % --- Temperature, salinity, and pressure profiles by vehicle ---
    g = vehicle_params.g;   % [m/s^2]
    
    switch upper(vehicle_params.name)
        case 'SCOUT'   % Gothenburg estuary (brackish, 0–30 m)
            % Anchor points for a mild halocline/thermocline
            z_knots = [0 5 10 15 20 30]';
            T_knots = [14 13 12 11  9  9]';    % [°C]
            S_knots = [18 20 22 23 24 25]';    % [PSU] brackish
    
            % Interpolate to your depth grid
            temp     = interp1(z_knots, T_knots, z, 'pchip', 'extrap');
            salinity = interp1(z_knots, S_knots, z, 'pchip', 'extrap');
    
            % Slightly lower reference density for hydrostatic p in estuary
            rho_ref  = 1015;   % [kg/m^3] representative brackish
    
        case 'REMUS'   % Open-ocean baseline (saltier, deeper)
            % Simple mixed layer + linear thermocline to z_max
            T_surf = 12;                  % [°C] surface
            T_deep = 6;                   % [°C] at z_max
            z_mld  = min(50, max(5, 0.25*z_max));   % [m] robust MLD choice
    
            temp = zeros(size(z));
            temp(z <= z_mld) = T_surf;
            if z_max > z_mld
                temp(z > z_mld) = T_surf + (T_deep - T_surf) .* ...
                                   (z(z>z_mld) - z_mld) ./ (z_max - z_mld);
            else
                temp(:) = T_surf;
            end
    
            salinity = 35*ones(size(z));  % [PSU] oceanic baseline
            rho_ref  = 1025;              % [kg/m^3] typical seawater
    
        otherwise
            % Fallback: oceanic-like defaults
            T_surf = 12; T_deep = 6; z_mld = min(50, max(5, 0.25*z_max));
            temp = zeros(size(z));
            temp(z <= z_mld) = T_surf;
            if z_max > z_mld
                temp(z > z_mld) = T_surf + (T_deep - T_surf) .* ...
                                   (z(z>z_mld) - z_mld) ./ (z_max - z_mld);
            else
                temp(:) = T_surf;
            end
            salinity = 35*ones(size(z));
            rho_ref  = 1025;
    end
    
    % Hydrostatic gauge pressure (p=0 at surface)
    pressure = rho_ref * g * z;   % [Pa]


   
    % Step 3: Calculate density at different depths using calc_rho_water
    [rho, rho_temp, rho_sal] = calc_rho_water(temp, salinity, pressure);

    % Step 4: Calculate buoyancy force and net force at each depth
    g     = vehicle_params.g;         % [m/s^2]
    m     = vehicle_params.mass;      % [kg]
    W     = m * g;                    % Weight in water (approx same as in air) [N]
    
    FB    = rho .* g .* V;            % Buoyant force [N]
    Fnet  = FB - W;                   % + = upward (positive buoyancy), - = downward
    % Step 5: Plot the results
    % --- Step 5: Plot the results ---
    % Neutral-buoyancy depth (if Fnet crosses zero)
    nz = NaN;
    if any(Fnet(1:end-1).*Fnet(2:end) <= 0)
        nz = interp1(Fnet, z, 0, 'linear');
    end
    
    figure('Name','Part A1: Buoyancy vs Depth','Color','b');
    
    subplot(1,3,1);
    plot(rho, z, 'LineWidth',1.5);
    set(gca,'YDir','reverse'); grid on;
    xlabel('\rho [kg/m^3]'); ylabel('Depth z [m]');
    title('Density profile');
    
    subplot(1,3,2);
    plot(FB, z, 'LineWidth',1.5); hold on;
    plot(W*ones(size(z)), z, '--', 'LineWidth',1.2);
    set(gca,'YDir','reverse'); grid on;
    xlabel('Force [N]'); ylabel('Depth z [m]');
    legend('Buoyancy F_B','Weight W','Location','best');
    title('Buoyancy vs weight');
    
    subplot(1,3,3);
    plot(Fnet, z, 'LineWidth',1.5); hold on;
    xline(0,'k:','LineWidth',1.0);
    if ~isnan(nz), yline(nz,'r--','Neutral depth','LabelHorizontalAlignment','left'); end
    set(gca,'YDir','reverse'); grid on;
    xlabel('F_{net} = F_B - W [N]'); ylabel('Depth z [m]');
    title('Net force');

end
