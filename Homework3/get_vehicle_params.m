function vehicle_params = get_vehicle_params(name, rho, g)
% Return physical parameters for a named AUV.
%
%
% Units: length [m], diameter [m], volume [m^3], mass [kg], depth [m],
% energy [kW-hr], endurance [h at m/s].
%
% Inputs:
%   name  - 'REMUS'
%   rho   - water density [kg/m^3]
%   g     - gravity [m/s^2]

    switch upper(name)
        case 'REMUS'
            % --- Appendix 1: REMUS AUV DATA ---
            length_m        = 1.60;      % Vehicle Length [m]
            diameter_m      = 0.19;      % Vehicle Diameter [m]
            CB              = 0.825;     % Block Coefficient [-]
            volume_m3       = 0.0315;    % Volume Displacement [m^3]
            mass_kg         = 37;        % Weight in air (treated as mass) [kg]
            trim_mass_kg    = 1;         % Trim weight in air [kg]
            max_depth_m     = 100;       % Maximum Operating Depth [m]
            energy_kWhr     = 1;         % Energy Capacity [kW-hr]
            endurance_h     = 22;        % Endurance [h] at the speed below
            endurance_speed = 1.5;       % Speed for stated endurance [m/s]
            velocity_range  = [0.25, 2.8];   % [m/s]
           
            % --- Package into struct ---
            vehicle_params.name             = 'REMUS';
            vehicle_params.length           = length_m;
            vehicle_params.diameter         = diameter_m;
            vehicle_params.CB               = CB;
            vehicle_params.volume           = volume_m3;
            vehicle_params.mass             = mass_kg;
            vehicle_params.trim_mass        = trim_mass_kg;
            vehicle_params.max_depth        = max_depth_m;
            vehicle_params.energy_kWhr      = energy_kWhr;
            vehicle_params.endurance_h      = endurance_h;
            vehicle_params.endurance_speed  = endurance_speed;
            vehicle_params.velocity_range   = velocity_range;
            

            % Inputs passed in:
            vehicle_params.rho = rho;
            vehicle_params.g   = g;

            % (Optional derived helpers you might need later; not from the appendix)
            vehicle_params.radius          = diameter_m/2;
            vehicle_params.frontal_area    = pi*(diameter_m/2)^2;  % for drag, if needed
        case 'SCOUT'
            length_m        = 1.30;       % Vehicle Length [m]
            diameter_m      = 0.16;       % Vehicle Diameter [m]
            CB              = 0.85;       % Block Coefficient [-]
            mass_kg         = 21.5;       % Weight in air (treated as mass) [kg]
            trim_mass_kg    = 0.5;        % Trim weight in air [kg]
            max_depth_m     = 30;         % Maximum Operating Depth [m]
            energy_kWhr     = 0.6;        % Energy Capacity [kW-hr]
            endurance_h     = 10;         % Endurance [h] at the speed below
            endurance_speed = 1.5;        % Speed for stated endurance [m/s]
            velocity_range  = [2, 4]*0.514444;   % knots converted to [m/s]
            
            % --- Package into struct ---
            vehicle_params.name             = 'SCOUT';
            vehicle_params.length           = length_m;
            vehicle_params.diameter         = diameter_m;
            vehicle_params.CB               = CB;
            vehicle_params.mass             = mass_kg;
            vehicle_params.trim_mass        = trim_mass_kg;
            vehicle_params.max_depth        = max_depth_m;
            vehicle_params.energy_kWhr      = energy_kWhr;
            vehicle_params.endurance_h      = endurance_h;
            vehicle_params.endurance_speed  = endurance_speed;
            vehicle_params.velocity_range   = velocity_range;
            
            % Inputs passed in:
            vehicle_params.rho = rho;
            vehicle_params.g   = g;

        
        otherwise
            error('Unknown vehicle name: %s', name);
    end
end
