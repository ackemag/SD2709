%% Task 1a — Atlantic sound-speed profile using GSW (TEOS-10) with WOA23
% Matlab copilot was used in helping construct this code
% Requires:
%   - GSW Oceanographic Toolbox on MATLAB path
%   - Files: woa23_decav_t00_01.nc, woa23_decav_s00_01.nc

clear; clc;

% -------- USER SETTINGS --------
tgt_lat = 30;            % deg N
tgt_lon = -45;           % deg E (negative = West)
z_out   = (0:100:8000)'; % depths for output (m)
save_plot = true;

tempFile = 'woa23_decav_t00_01.nc';
saltFile = 'woa23_decav_s00_01.nc';

% -------- READ COORDS --------
lat = ncread(tempFile,'lat');      % [nlat]
lon = ncread(tempFile,'lon');      % [nlon] (often 0..360)
z   = ncread(tempFile,'depth');    % WOA standard depths (~0..5500 m)

% Map target lon into file's convention
if all(lon >= 0), tgt_lon_file = mod(tgt_lon,360); else, tgt_lon_file = tgt_lon; end

[~,ilat] = min(abs(lat - tgt_lat));
[~,ilon] = min(abs(lon - tgt_lon_file));

% -------- READ T & S COLUMNS --------
% WOA variable names:
%   t_an: in-situ temperature (deg C)
%   s_an: Practical Salinity (unitless PSS-78)
T3 = ncread(tempFile,'t_an');  % [lon x lat x depth]
S3 = ncread(saltFile,'s_an');  % [lon x lat x depth]

t = squeeze(T3(ilon, ilat, :));  % deg C (in-situ)
SP = squeeze(S3(ilon, ilat, :)); % Practical Salinity (PSS-78)

% Clean fill values if any (very negative placeholders)
t(t < -9e3) = NaN; SP(SP < -9e3) = NaN;
if any(isnan(t)),  t  = fillmissing(t,'linear');  end
if any(isnan(SP)), SP = fillmissing(SP,'linear'); end

% -------- EXTEND TO 8000 m --------
z_max = max(z);
t_ext  = interp1(z, t,  min(z_out,z_max),  'linear','extrap');   % hold below last level
SP_ext = interp1(z, SP, min(z_out,z_max), 'linear','extrap');

% -------- PRESSURE FROM DEPTH (GSW expects z negative downward) --------
z_gsw = -z_out;                          % meters, negative below sea surface
p = gsw_p_from_z(z_gsw, tgt_lat);        % dbar

% -------- CONVERT SP -> SA -> CT (TEOS-10) --------
% Need lon/lat for Absolute Salinity (spatial anomalies)
% GSW expects lon in 0..360
lon_gsw = mod(tgt_lon,360);
lat_gsw = tgt_lat;

SA = gsw_SA_from_SP(SP_ext, p, lon_gsw*ones(size(p)), lat_gsw*ones(size(p)));
CT = gsw_CT_from_t(SA, t_ext, p);

% -------- SOUND SPEED (two equivalent ways) --------
% 1) From SA, CT, p:
c = gsw_sound_speed(SA, CT, p);                 % m/s
% (Alternatively) 2) From SA, in-situ t, p:
% c = gsw_sound_speed_t_exact(SA, t_ext, p);

% -------- FIND SOFAR-LIKE MIN --------
[cmin, idxmin] = min(c);
zmin = z_out(idxmin);

%% --- PLOTTING: SA & CT on top, Sound Speed underneath (2x2 layout) ---
f = figure('Color','b');
t = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

% --- Top-left: Absolute Salinity (SA) ---
ax1 = nexttile(t,1);
plot(SA, z_out, 'LineWidth',1.6);
set(ax1,'YDir','reverse'); grid(ax1,'on');
xlabel('Absolute Salinity, SA (g/kg)');
ylabel('Depth (m)');
title('Absolute Salinity (TEOS-10)');

% --- Top-right: Conservative Temperature (CT) ---
ax2 = nexttile(t,2);
plot(CT, z_out, 'LineWidth',1.6);
set(ax2,'YDir','reverse'); grid(ax2,'on');
xlabel('Conservative Temperature, CT (°C)');
ylabel('Depth (m)');
title('Conservative Temperature (TEOS-10)');

% --- Bottom: Sound speed (span two columns) ---
ax3 = nexttile(t,[1 2]); % span across both columns
plot(c, z_out, 'LineWidth',1.8);
set(ax3,'YDir','reverse'); grid(ax3,'on');
xlabel('Sound speed, c (m/s)');
ylabel('Depth (m)');
title(sprintf('GSW Sound-Speed Profile @ (%.1f^{\\circ}N, %.1f^{\\circ}E)', tgt_lat, tgt_lon));

% Mark SOFAR-like minimum if you computed cmin/zmin earlier
hold(ax3,'on');
[cmin, idxmin] = min(c);
zmin = z_out(idxmin);
plot(ax3, cmin, zmin, 'o', 'MarkerSize',6, 'MarkerFaceColor',[0.2 0.2 0.2]);
text(ax3, cmin+2, zmin, sprintf('c_{min}=%.1f m/s @ %dm', cmin, round(zmin)), ...
     'VerticalAlignment','middle');

% Optional: save
exportgraphics(f,'Atlantic_SA_CT_SSP.png','Resolution',200);
