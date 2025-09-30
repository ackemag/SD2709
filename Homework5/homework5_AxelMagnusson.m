% SD2709 – Homework 5: Maneuvering & Stability
% Axel Magnusson
clc; clear; close all;

%% -----------------------------------------------------------------------------
% Inputs (pick REMUS or SCOUT)
rho  = 1025;            % seawater density [kg/m^3] (1015 SCOUT, 1025 REMUS)
g    = 9.81;            % gravity [m/s^2]
name = 'REMUS';         % 'REMUS' or 'SCOUT'
N    = 150;             % axial stations for strip integrals

% Vehicle data
[a,b,c,d,aofs,cofs,n,theta,L,Vol,B,W,m,Ixx,Iyy,Izz,cg_org,cb_org,Ud] = ...
    get_vehicle_params(name, rho, g);

%% -----------------------------------------------------------------------------
% Myring hull discretization (x from nose -> tail), then shift to CG frame
dx_nom = L/N;
[x_raw, r_raw] = myringShape(a,aofs,b,c,cofs,n,theta,d,dx_nom);
x = x_raw - cg_org(1);                        % x relative to CG
d_local = 2*r_raw;
A_local = (pi/4)*d_local.^2;
m_a = rho * A_local;

% Strip integrals
I0 = trapz(x,               m_a);
I1 = trapz(x,      x     .* m_a);
I2 = trapz(x, (x.^2) .*    m_a);

U = Ud;   % design speed

%% -----------------------------------------------------------------------------
% Task 1: Hydrodynamic derivatives (hull + fins)
Yv_h = -U*I0;   Yr_h =  U*I1;   Nv_h =  U*I1;   Nr_h = -U*I2;
Zw_h = -U*I0;   Zq_h =  U*I1;   Mw_h =  U*I1;   Mq_h = -U*I2;

S_fin_single = 6.65e-3;   S_pair = 2*S_fin_single;
a3D = 3.12;               xT = abs(-0.638);

Yv_f =  0.5*rho*U*a3D*S_pair;   Yr_f = -0.5*rho*U*a3D*S_pair*xT;
Nv_f =  Yv_f*xT;                Nr_f =  Yr_f*xT;
Zw_f =  0.5*rho*U*a3D*S_pair;   Zq_f = -0.5*rho*U*a3D*S_pair*xT;
Mw_f =  Zw_f*xT;                Mq_f =  Zq_f*xT;

Yv = Yv_h + Yv_f;   Yr = Yr_h + Yr_f;   Nv = Nv_h + Nv_f;   Nr = Nr_h + Nr_f;
Zw = Zw_h + Zw_f;   Zq = Zq_h + Zq_f;   Mw = Mw_h + Mw_f;   Mq = Mq_h + Mq_f;

%% -----------------------------------------------------------------------------
% Plot 1: Myring hull with CG and CB
figure;
plot(x_raw,  r_raw,'r','LineWidth',2); hold on;
plot(x_raw, -r_raw,'r','LineWidth',2);
plot(cg_org(1), -cg_org(3),'w+','MarkerSize',8,'LineWidth',1.5);
plot(cb_org(1),  cb_org(3),'bo','MarkerSize',6,'LineWidth',1.5);
xlabel('x (m)'); ylabel('radius (m)');
title(sprintf('%s hull (Myring shape)', upper(name)));
legend('Hull','Hull','CG','CB','Location','best');
axis equal; grid on;

%% -----------------------------------------------------------------------------
% Task 2: Lateral stability
Ydotv = -I0;  Ndotr = -I2;
m_eff = m - Ydotv;
Iz_eff = Izz - Ndotr;
Kconv = Yr - m*U;

A = [Yv/m_eff, Kconv/m_eff;
     Nv/Iz_eff, Nr/Iz_eff];
trA = trace(A); detA = det(A); eigA = eig(A);

%% -----------------------------------------------------------------------------
% Task 3: Critical forward velocity
BG = abs(cg_org(3)-cb_org(3));
delta_deg = 10; delta = deg2rad(delta_deg);
Vc = sqrt( (2*m*g*BG) / (rho * a3D * delta * S_pair * xT) );

%% -----------------------------------------------------------------------------
% Task 4: Turning radii
Y_delta = 0.5*rho*U^2*a3D*S_pair;
N_delta = xT*Y_delta;
den = Yv*Nr - Kconv*Nv;

delta_r_deg = [5 10 15];
delta_r = deg2rad(delta_r_deg);
R_out = zeros(size(delta_r));

for k = 1:numel(delta_r)
    r_ss = delta_r(k) * (Nv*Y_delta - N_delta*Yv) / den;
    R_out(k) = U / abs(r_ss);
end

% Plot 2: Turning circles
theta = linspace(0,2*pi,400);
figure; hold on; axis equal; grid on;
for k = 1:numel(R_out)
    R = R_out(k);
    xC = R*sin(theta); 
    yC = R*(1-cos(theta));
    plot(xC,yC,'LineWidth',2, ...
        'DisplayName',sprintf('\\delta_r = %d° → R = %.1f m',delta_r_deg(k),R));
end
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Turning circles at U = %.2f m/s', U));
legend('Location','best'); box on;

%% -----------------------------------------------------------------------------
% Print results
fprintf('Vehicle: %s\n', upper(name));
fprintf('I0 = %.3f, I1 = %.4f, I2 = %.5f\n', I0,I1,I2);
fprintf('Derivatives (total): Yv=%.2f, Yr=%.2f, Nv=%.2f, Nr=%.2f\n',Yv,Yr,Nv,Nr);
fprintf('                    Zw=%.2f, Zq=%.2f, Mw=%.2f, Mq=%.2f\n',Zw,Zq,Mw,Mq);
fprintf('\nLateral stability: trace=%.3e, det=%.3e,\n', trA,detA);
fprintf('Critical speed Vc = %.3f m/s (δ_e = %.1f°)\n', Vc, delta_deg);
for k=1:numel(delta_r_deg)
    fprintf('Turning radius delta_r=%2d° → R=%.2f m\n', delta_r_deg(k), R_out(k));
end
