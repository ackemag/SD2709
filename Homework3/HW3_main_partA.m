function HW3_main_partA(vehicle_choice)
    % Main script for HW3 parts A1, A2, A3
    % Choose SCOUT for part B
    
    rho_water = 1025;  % kg/m^3
    g = 9.81;  % gravitational acceleration (m/s^2)

      % ---- Vehicle selection: 1=REMUS, 2=SCOUT ----
    if nargin < 1 || isempty(vehicle_choice)
        fprintf('Select vehicle:\n  1) REMUS\n  2) SCOUT\n');
        vehicle_choice = input('Enter 1 or 2: ');
    end

    switch vehicle_choice
        case 1
            vehicle_name = 'REMUS';
        case 2
            vehicle_name = 'SCOUT';
        otherwise
            error('Invalid selection. Use 1 for REMUS or 2 for SCOUT.');
    end

    % Get vehicle parameters
    vehicle_params = get_vehicle_params(vehicle_name, rho_water, g);
    
    % Run each part
    part_A1(vehicle_params);
    part_A2(vehicle_params);
    part_A3(vehicle_params);
end
