function part_A2(vehicle_params)
    % Part A2: Hydromechanics and Drag Calculation for the AUV using slenderness ratio analysis.

    % Use vehicle-specific speed range if available
    if isfield(vehicle_params, 'velocity_range') && numel(vehicle_params.velocity_range) == 2
        Umin = vehicle_params.velocity_range(1);
        Umax = vehicle_params.velocity_range(2);
    else
        if isfield(vehicle_params,'name') && strcmpi(vehicle_params.name,'SCOUT')
            Umin = 2 * 0.514444;   % 2 kn -> m/s
            Umax = 4 * 0.514444;   % 4 kn -> m/s
        else
            Umin = 0.25;           % REMUS baseline from HW3
            Umax = 2.80;
        end
    end

    % Build speed vector
    U = linspace(min(Umin,Umax), max(Umin,Umax), 50);

    % Step 1: Define slenderness ratios (range from 12 to 2)
    f_vec = linspace(12, 2, 25)';                 % slenderness sweep (long -> stubby)
    CB    = vehicle_params.CB;
    
    % Pick the fixed displacement volume:
    L0 = vehicle_params.length;
    D0 = vehicle_params.diameter;
    R0 = D0/2;
    
    V_geom0 = pi*R0^2*(L0 - D0) + (4/3)*pi*R0^3; % cylinder + two hemispheres
    V_calc0 = CB * V_geom0;                       % computed baseline volume
    
    if isfield(vehicle_params,'volume_disp') && ~isempty(vehicle_params.volume_disp)
        Vfix = vehicle_params.volume_disp;        % use REMUS given (0.0315 m^3)
    elseif isfield(vehicle_params,'volume') && ~isempty(vehicle_params.volume)
        Vfix = vehicle_params.volume;             % use precomputed if present
    else
        Vfix = V_calc0;                           % compute for SCOUT
    end
    
    % From Appendix 2 with L = f*D:
    %   V = CB * π * D^3 * (3f - 1) / 12  => solve for D(f):
    Df = ((12*Vfix) ./ (CB*pi*(3*f_vec - 1))).^(1/3);   % diameter for each f
    Lf = f_vec .* Df;                                   % length for each f

    
    % Step 2: For each slenderness ratio, calculate drag forces and coefficients
    % Call calc_drag_force for each updated vehicle
    rho = vehicle_params.rho;   % water density [kg/m^3]
    nu  = 1.19e-6;              % kinematic viscosity [m^2/s] (from appendix / notes)
    
    nF = length(f_vec);
    nU = length(U);
    
    drag_matrix   = zeros(nU, nF);   % drag vs U,f
    Cd_matrix     = zeros(nU, nF);   % Cd_total vs U,f
    
    for j = 1:nF
        % Build a temporary vehicle with updated length and diameter
        veh_tmp          = vehicle_params;
        veh_tmp.length   = Lf(j);
        veh_tmp.diameter = Df(j);
    
        % Call your drag function
        [drag_forces, Cd_total, f, Cf, Cd_body] = ...
            calc_drag_force(U, veh_tmp, rho, nu);
    
        drag_matrix(:,j) = drag_forces(:);
        Cd_matrix(:,j)   = Cd_total(:);
    end
    % Step 3: Plot the results
    % (a) Drag vs speed for selected slenderness ratios
    figure('Name','Part A2: Drag vs Speed','Color','b');
    hold on; grid on;
    plot_indices = round(linspace(1, nF, 4));   % pick 4 representative f values
    for k = plot_indices
        plot(U, drag_matrix(:,k), 'LineWidth',1.5, ...
            'DisplayName', sprintf('f = %.1f', f_vec(k)));
    end
    xlabel('Speed U [m/s]');
    ylabel('Drag Force [N]');
    title('Drag force vs speed for different slenderness ratios');
    legend show;
    
    % (b) Surface plot: drag vs U and f
    figure('Name','Part A2: Drag Surface','Color','b');
    surf(f_vec, U, drag_matrix, 'EdgeColor','none');
    xlabel('Slenderness ratio f = L/D');
    ylabel('Speed U [m/s]');
    zlabel('Drag Force [N]');
    title('Drag force as function of speed and slenderness ratio');
    colorbar;
    view(135,30);
end
