function part_A3(vehicle_params, nu_water, rho_water, U)
% Thrust, power, range, endurance vs speed for chosen AUV.
% Matlab co-pilot was used in helping create the code, mainly for plots

    % -------- Speed vector from vehicle definition --------
    if nargin < 4 || isempty(U)
        if isfield(vehicle_params,'velocity_range') && numel(vehicle_params.velocity_range)==2
            U = linspace(vehicle_params.velocity_range(1), ...
                         vehicle_params.velocity_range(2), 20);
        else
            if isfield(vehicle_params,'name') && strcmpi(vehicle_params.name,'SCOUT')
                U = linspace(2,4,20) * 0.514444;  % 2–4 kn → m/s
            else
                U = linspace(0.25, 2.8, 20);      % REMUS baseline
            end
        end
    end
    U = U(:)';  % row vector

    % -------- Fluid properties (used for drag) --------
    if nargin < 2 || isempty(nu_water),  nu_water  = 1.19e-6; end
    if nargin < 3 || isempty(rho_water), rho_water = vehicle_params.rho; end

    % -------- Drag -> Thrust --------
    [drag_N, Cd_total, f, Cf, Cd_body] = calc_drag_force(U, vehicle_params, rho_water, nu_water); %#ok<ASGLU>
    T = drag_N;   % steady speed: thrust = drag

    % -------- Propulsive efficiency, hotel load, battery energy --------
    eta_min = 0.30; 
    eta_max = 0.73;

    if isfield(vehicle_params,'eta_P') && ~isempty(vehicle_params.eta_P)
        eta_P = min(max(vehicle_params.eta_P, eta_min), eta_max);
    else
        % Use a mid-range default if not provided
        eta_P = 0.54;  % around the middle of 0.30–0.73
    end


    % hotel load
    if isfield(vehicle_params,'hotel_W') && ~isempty(vehicle_params.hotel_W)
        P_hotel = vehicle_params.hotel_W;
    else
        if isfield(vehicle_params,'name') && strcmpi(vehicle_params.name,'SCOUT')
            P_hotel = 35;   % SCOUT baseline (from literature review)
        else
            P_hotel = 77;   % REMUS baseline (from lecture slides)
        end
    end

    % battery energy
    if isfield(vehicle_params,'energy_kWhr') && ~isempty(vehicle_params.energy_kWhr)
        E_bat_J = vehicle_params.energy_kWhr * 3.6e6;
    elseif isfield(vehicle_params,'energy_Wh') && ~isempty(vehicle_params.energy_Wh)
        E_bat_J = vehicle_params.energy_Wh * 3600;
    elseif isfield(vehicle_params,'energy_J') && ~isempty(vehicle_params.energy_J)
        E_bat_J = vehicle_params.energy_J;
    else
        if isfield(vehicle_params,'name') && strcmpi(vehicle_params.name,'SCOUT')
            E_bat_J = 0.6 * 3.6e6;   % 0.6 kWh
        else
            E_bat_J = 1.0 * 3.6e6;   % 1.0 kWh
        end
    end

    % -------- Power, endurance, range --------
    P_prop = (T .* U) ./ max(eta_P, eps);   % [W]
    P_tot  = P_prop + P_hotel;                 % [W]

    endurance_s = E_bat_J ./ P_tot;            % [s]
    endurance_h = endurance_s / 3600;          % [h]
    range_km    = (U .* endurance_s) / 1000;   % [km]

    
    % -------- Plots (η sensitivity) --------
    eta_vec = [0.3, 0.4, 0.5, 0.6, 0.7];
    colors = lines(numel(eta_vec));

    % (1) Thrust vs Speed
    figure('Name','Thrust vs Speed (η sensitivity)','Color','b'); hold on; grid on;
    plot(U, T, 'r-', 'LineWidth',1.5, 'DisplayName','T = D (ideal)');
    for k = 1:numel(eta_vec)
        eta_k = eta_vec(k);
        plot(U, T./eta_k, '--', 'Color',colors(k,:), 'LineWidth',1.5, ...
            'DisplayName',sprintf('Effective thrust (η = %.1f)',eta_k));
    end
    xlabel('Speed U [m/s]'); ylabel('Thrust [N]');
    title(sprintf('Thrust vs Speed (%s)', iff_name(vehicle_params)));
    legend('Location','northwest');

    % (2) Range vs Speed
    figure('Name','Range vs Speed (η sensitivity)','Color','b'); hold on; grid on;
    for k = 1:numel(eta_vec)
        eta_k = eta_vec(k);
        P_prop_k = (T .* U) ./ eta_k;
        P_tot_k  = P_prop_k + P_hotel;
        range_km_k = (U .* (E_bat_J ./ P_tot_k)) / 1000;

        plot(U, range_km_k, '-', 'Color',colors(k,:), 'LineWidth',1.5, ...
            'DisplayName',sprintf('η = %.1f',eta_k));

        [max_range, idx] = max(range_km_k);
        plot(U(idx), max_range, 'ro', 'MarkerSize',8, 'LineWidth',1.2);
        text(U(idx), max_range, sprintf('  Max @ %.2f m/s', U(idx)), ...
            'VerticalAlignment','bottom','FontSize',8);
    end
    xlabel('Speed U [m/s]'); ylabel('Range [km]');
    title('Range vs Speed for different propeller efficiencies');
    legend('Location','best');

    % (3) Endurance vs Speed
    figure('Name','Endurance vs Speed (η sensitivity)','Color','b'); hold on; grid on;
    for k = 1:numel(eta_vec)
        eta_k = eta_vec(k);
        P_prop_k = (T .* U) ./ eta_k;
        P_tot_k  = P_prop_k + P_hotel;
        endurance_h_k = (E_bat_J ./ P_tot_k) / 3600;

        plot(U, endurance_h_k, '-', 'Color',colors(k,:), 'LineWidth',1.5, ...
            'DisplayName',sprintf('η = %.1f',eta_k));

        [max_end, idx] = max(endurance_h_k);
        plot(U(idx), max_end, 'ro', 'MarkerSize',8, 'LineWidth',1.2);
        text(U(idx), max_end, sprintf('  Max @ %.2f m/s', U(idx)), ...
            'VerticalAlignment','bottom','FontSize',8);
    end
    xlabel('Speed U [m/s]'); ylabel('Endurance [h]');
    title('Endurance vs Speed for different propeller efficiencies');
    legend('Location','best');

    % (4) Combined Range & Endurance vs Speed
    figure('Name','Range & Endurance vs Speed (η sensitivity)','Color','b'); hold on; grid on;
    for k = 1:numel(eta_vec)
        eta_k = eta_vec(k);
        P_prop_k = (T .* U) ./ eta_k;
        P_tot_k  = P_prop_k + P_hotel;
        endurance_h_k = (E_bat_J ./ P_tot_k) / 3600;
        range_km_k    = (U .* (E_bat_J ./ P_tot_k)) / 1000;

        yyaxis left
        plot(U, range_km_k, '-', 'Color',colors(k,:), 'LineWidth',1.5, ...
            'DisplayName',sprintf('Range (η = %.1f)',eta_k));
        ylabel('Range [km]');

        yyaxis right
        plot(U, endurance_h_k, '--', 'Color',colors(k,:), 'LineWidth',1.5, ...
            'HandleVisibility','off');
        ylabel('Endurance [h]');
    end
    xlabel('Speed U [m/s]');
    title('Range (solid) and Endurance (dashed) vs Speed for different η');
    yyaxis left

    
end

% ---- tiny helper for safe name printing ----
function nm = iff_name(vp)
    if isstruct(vp) && isfield(vp,'name') && ~isempty(vp.name)
        nm = vp.name;
    else
        nm = 'AUV';
    end
end
