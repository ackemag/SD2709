function [x_combined, r_combined] = myringShape(a, aofs, b, c, cofs, n, theta, d, dx)
    % Total length of the vehicle
    l = a + b + c - aofs - cofs;
    
    % Define x based on dx
    x = 0:dx:l;
    
    % Nose section
    r_nose = (d / 2) * (1 - ((x + aofs - a) / a).^2).^(1 / n);
    nose_ind = find(x <= a);
    r_nose = r_nose(nose_ind);
    x_nose = x(nose_ind);
    
    % Nose geometric properties
    As_nose = calc_surface_area(x_nose, r_nose);
    Ap_nose = calc_plane_area(x_nose, r_nose);
    V_nose = calc_volume(x_nose, r_nose);
    CB_nose = calc_cb(x_nose, r_nose);
    
    % Tail section
    r_tail = (d / 2) - (3 * d / (2 * (l - a - b)^2) - tan(theta) / (l - a - b)) * (x - a - b).^2 + ...
             (d / (l - a - b)^3 - tan(theta) / (l - a - b)^2) * (x - a - b).^3;
    tail_ind = find(x >= (a + b));
    r_tail = r_tail(tail_ind);
    x_tail = x(tail_ind);
    
    % Tail geometric properties
    As_tail = calc_surface_area(x_tail, r_tail);
    Ap_tail = calc_plane_area(x_tail, r_tail);
    V_tail = calc_volume(x_tail, r_tail);
    CB_tail = calc_cb(x_tail, r_tail);
    
    % Midsection
    r_mids = [d / 2, d / 2];
    x_mids = [a, a + b];
    
    % Midsection geometric properties
    As_mids = pi * d * b;
    Ap_mids = d * b;
    V_mids = pi * (d / 2)^2 * b;
    CB_mids = b / 2;

    % Total properties
    As = As_mids + As_nose + As_tail;
    Ap = Ap_mids + Ap_nose + Ap_tail;
    Af = pi * (d / 2)^2;
    V = V_mids + V_nose + V_tail;
    cb = (V_nose * CB_nose + V_mids * (a + CB_mids) + V_tail * (a + b + CB_tail)) / V;
    
    % Concatenate all sections
    r_combined = [r_nose, r_mids, r_tail];
    x_combined = [x_nose, x_mids, x_tail];
end

% Additional geometric calculation functions

function As = calc_surface_area(x, r)
    dx = diff(x);
    r = r(1:end-1);
    dA = pi * 2 * r .* dx;
    As = sum(dA);
end

function A = calc_plane_area(x, r)
    dx = diff(x);
    r = r(1:end-1);
    dA = 2 * r .* dx;
    A = sum(dA);
end

function V = calc_volume(x, r)
    dx = diff(x);
    r = r(1:end-1);
    dV = pi * r.^2 .* dx;
    V = sum(dV);
end

function CB = calc_cb(x, r)
    dx = diff(x);
    x = x(1:end-1);
    r = r(1:end-1);
    dV = pi * r.^2 .* dx;
    V = sum(dV);
    CB = sum(x .* dV) / V;
end
