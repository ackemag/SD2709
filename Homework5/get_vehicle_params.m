function [a, b, c, d, aofs, cofs, n, theta, l, Vol, B, W, m, Ixx, Iyy, Izz, cg, cb, Ud] = get_vehicle_params(name, rho, g)
    if strcmp(name, 'REMUS')
        % Vehicle parameters REMUS
        a = 191e-3;  % Nose length
        b = 654e-3;  % Mid-body length
        c = 541e-3;  % Tail length
        d = 191e-3;  % Diameter
        
        aofs = 16.5e-3;  % Truncated nose
        cofs = 36.8e-3;  % Truncated tail
        n = 2;  % Exponential coefficient
        theta = 4.36e-1;  % Tail angle

        l = a + b + c - aofs - cofs;  % Total length

        Vol = 3.15e-2;  % Hull volume in m^3
        B = Vol * rho * g;  % Buoyancy force in N

        W = 299;  % Weight in N
        m = W / g;  % Mass in kg

        % Moments of inertia
        Ixx = 1.77e-1;
        Iyy = 3.45;
        Izz = 3.45;

        % Center of gravity and center of buoyancy
        cg = [611e-3, 0, 19.6e-3];  % Center of gravity
        cb = [611e-3, 0, 0];  % Center of buoyancy

        % Design speed
        Ud = 1.5;  % m/s

    elseif strcmp(name, 'AutoSub3')
        % Vehicle parameters AutoSub3
        a = 1.0;  % Nose length
        b = 3.5;  % Mid-body length
        c = 2.5;  % Tail length
        d = 0.9;  % Diameter
        
        aofs = 0;  % Truncated nose
        cofs = 0;  % Truncated tail
        n = 2;  % Exponential coefficient
        theta = 4.36e-1;  % Tail angle

        l = a + b + c;  % Total length

        Vol = 7.96;  % Hull volume in m^3
        B = Vol * rho * g;  % Buoyancy force in N

        W = 7672.02;  % Weight in N
        m = W / g;  % Mass in kg

        % Moments of inertia
        Ixx = 515;
        Iyy = 7368;
        Izz = 7368;

        % Center of gravity and center of buoyancy
        cg = [3.0, 0, 0];  % Center of gravity
        cb = [3.0, 0, 0];  % Center of buoyancy

        % Design speed
        Ud = 2.5;  % m/s

    % Add your vehicle here by copying the above elseif block and changing the parameters
    elseif strcmp(name, 'SCOUT')
        % ========================== Vehicle parameters SCOUT ==========================
        % Design: Compact harbor AUV for CTD monitoring (Göteborg harbor)
        % References: REMUS-100, Graal Tech X300, Lo-MARVE, MIT Hydrofoil notes

        % --- Geometry ---
        L = 1.30;              % [m] overall length (similar to REMUS-100: 1.33 m)
        D = 0.16;              % [m] max diameter (scaled smaller than REMUS-100: 0.19 m)
        CB = 0.85;             % [-] block coefficient (typical streamlined torpedo value)

        a = 0.20*L;            % [m] nose length (20% of L)
        b = 0.60*L;            % [m] mid-body length (60% of L)
        c = 0.20*L;            % [m] tail length (20% of L)
        d = D;                 % [m] diameter

        aofs = 0.0;            % [m] truncated nose (none)
        cofs = 0.0;            % [m] truncated tail (none)
        n = 2;                 % [-] Myring exponent (2 = ellipsoidal)
        theta = 4.36e-1;       % [rad] tail angle (same convention as REMUS)

        l = a + b + c - aofs - cofs;   % [m] modeled total length

        % --- Hydrostatics ---
        Vol = CB * pi * (D/2)^2 * L;  % [m^3] displaced volume
        B   = Vol * rho * g;          % [N] buoyancy force

        m = 22.0;                     % [kg] total mass (portable scale, cf. X300 AUV)
        W = m * g;                    % [N] weight

        % --- Inertia (slender-body approximations) ---
        R   = D/2;
        Ixx = 0.5 * m * R^2;                 % [kg·m^2] roll
        Iyy = (m/12) * (L^2 + 3*R^2);        % [kg·m^2] pitch
        Izz = Iyy;                           % [kg·m^2] yaw

        % --- Centers ---
        x_cg = 0.60 * L;                     % [m] CG ~60% L from nose
        z_cg = +0.015;                       % [m] CG slightly above CB (BG=0.015 m)
        cg = [x_cg, 0, z_cg];                % [m] center of gravity
        cb = [x_cg, 0, 0.000];               % [m] center of buoyancy on centerline

        % --- Design / operating speed ---
        Ud = 1.5;                            % [m/s] design speed (typical for survey AUVs)


    % elseif strcmp(name, 'YourVehicle')
    %   a = ...;  % Nose length
    %   b = ...;  % Mid-body length
    %   % Set other parameters for your custom vehicle
    else
        error('Unknown vehicle: %s', name);
    end
end
