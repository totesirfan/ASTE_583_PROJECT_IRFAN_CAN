% Phase1_Main_RefTraj.m
% MASTER SCRIPT FOR PHASE 1
% - Propagates Reference Trajectory (Sun-Centered EMO2000)
% - Applies Lunar Targeting Maneuver (LTM)
% - Verifies Lunar Flyby
% - Generates 3 Plots: Heliocentric, Geocentric, Ground Track

clear; clc; close all;

%% 1. Setup & Constants
init_project();
const = lib_constants();

t0 = cspice_str2et(const.epoch_utc_str);

% Time Definitions
t_LTM = cspice_str2et(const.LTM.date_utc);
tf    = t0 + (251 * const.day2sec);

%% 2. Propagation (Sun-Centered EMO2000)
% State: [r(3); v(3); k_SRP(1); bias(1)]
X0 = [const.X0_ref; 1.0; 0.0];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

fprintf('1/4 Propagating Segment 1 (Detection -> LTM)...\n');
[T1, X1] = ode45(@(t,x) lib_dynamics(t, x, const), [t0 t_LTM], X0, options);

% Apply LTM Maneuver
state_at_LTM = X1(end, :)';
state_at_LTM(4:6) = state_at_LTM(4:6) + const.LTM.dV;
fprintf('    * Applied LTM Delta-V: [%.4f, %.4f, %.4f] km/s\n', const.LTM.dV);

fprintf('2/4 Propagating Segment 2 (LTM -> Flyby)...\n');
[T2, X2] = ode45(@(t,x) lib_dynamics(t, x, const), [t_LTM tf], state_at_LTM, options);

% Combine Trajectory
T_vec = [T1; T2];
X_sc_sun = [X1; X2];
n_steps = length(T_vec);

%% 3. Coordinate Transformations (For Plotting)
fprintf('3/4 Processing Coordinates (Earth, Moon, Ground Track)...\n');

r_earth_sun = zeros(n_steps, 3);
r_sc_earth  = zeros(n_steps, 3);
r_moon_earth = zeros(n_steps, 3);
lats = zeros(n_steps, 1);
lons = zeros(n_steps, 1);

for i = 1:n_steps
    t = T_vec(i);
    r_sc_now = X_sc_sun(i, 1:3)';
    
    % A. Earth w.r.t Sun (J2000 -> EMO2000)
    st_earth = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_J2000 = st_earth(1:3);
    r_E_Sun_EMO = const.R_EME_EMO * r_E_Sun_J2000;
    r_earth_sun(i, :) = r_E_Sun_EMO';
    
    % B. Spacecraft w.r.t Earth (EMO2000)
    r_sc_earth_emo = r_sc_now - r_E_Sun_EMO;
    r_sc_earth(i, :) = r_sc_earth_emo';
    
    % C. Moon w.r.t Earth (J2000 -> EMO2000)
    st_moon = cspice_spkezr('MOON', t, 'J2000', 'NONE', 'EARTH');
    r_moon_J2000 = st_moon(1:3);
    r_moon_EMO = const.R_EME_EMO * r_moon_J2000;
    r_moon_earth(i, :) = r_moon_EMO';
    
    % D. Ground Track Calculations (ECF)
    % 1. Rotate SC from EMO2000 -> ECI (J2000)
    r_sc_eci = const.R_EME_EMO' * r_sc_earth_emo;
    
    % 2. Rotate ECI -> ECF (Body Fixed)
    % FIX: Use Updated GST Calculation relative to Detection Epoch
    phi_G = const.phi_G_detect + const.we * (t - const.t_detect_et);
    
    R_ECI_ECF = [cos(phi_G), sin(phi_G), 0;
                -sin(phi_G), cos(phi_G), 0;
                 0,          0,          1];
    
    r_sc_ecf = R_ECI_ECF * r_sc_eci;
    
    % 3. Cartesian -> Spherical (Lat/Lon)
    [lon, lat, ~] = cart2sph(r_sc_ecf(1), r_sc_ecf(2), r_sc_ecf(3));
    lats(i) = lat * const.rad2deg;
    lons(i) = lon * const.rad2deg;
end

%% 4. Visualization
fprintf('4/4 Generating Plots...\n');

% --- FIGURE 1: HELIOCENTRIC VIEW ---
figure('Name', 'Heliocentric View', 'Color', 'w', 'Position', [50 50 600 500]);
hold on; grid on; axis equal;
plot3(0, 0, 0, 'o', 'Color', [1 0.5 0], 'MarkerSize', 12, 'MarkerFaceColor', 'y'); % Sun
plot3(r_earth_sun(:,1), r_earth_sun(:,2), r_earth_sun(:,3), 'b--', 'LineWidth', 1); % Earth Orbit
plot3(r_earth_sun(1,1), r_earth_sun(1,2), r_earth_sun(1,3), 'bo', 'MarkerFaceColor', 'b'); % Earth Start
plot3(X_sc_sun(:,1), X_sc_sun(:,2), X_sc_sun(:,3), 'm-', 'LineWidth', 1.5); % SC
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Heliocentric Trajectory (Sun-Centered EMO2000)');
legend('Sun', 'Earth Orbit', 'Earth Start', 'SC Trajectory');
view(2);

% --- FIGURE 2: GEOCENTRIC VIEW (FLYBY) ---
figure('Name', 'Geocentric View', 'Color', 'w', 'Position', [660 50 600 500]);
hold on; grid on; axis equal;
[xE,yE,zE] = sphere(20);
surf(xE*const.R_E, yE*const.R_E, zE*const.R_E, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none'); % Earth
plot3(r_moon_earth(:,1), r_moon_earth(:,2), r_moon_earth(:,3), 'k:', 'LineWidth', 0.5); % Moon Orbit
plot3(r_sc_earth(:,1), r_sc_earth(:,2), r_sc_earth(:,3), 'm-', 'LineWidth', 2); % SC

% Highlight LTM
ltm_idx = length(T1);
plot3(r_sc_earth(ltm_idx,1), r_sc_earth(ltm_idx,2), r_sc_earth(ltm_idx,3), 'r^', 'MarkerFaceColor', 'r');

% Closest Approach Calculation
dist_to_moon = sqrt(sum((r_sc_earth - r_moon_earth).^2, 2));
[min_dist, idx_ca] = min(dist_to_moon);
time_ca = T_vec(idx_ca);
utc_ca = cspice_et2utc(time_ca, 'C', 3);

% Plot Moon & SC at Closest Approach
plot3(r_moon_earth(idx_ca,1), r_moon_earth(idx_ca,2), r_moon_earth(idx_ca,3), 'o', ...
    'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);
plot3(r_sc_earth(idx_ca,1), r_sc_earth(idx_ca,2), r_sc_earth(idx_ca,3), 'rp', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'm');

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Geocentric View & Lunar Flyby');
legend('Earth', 'Moon Orbit', 'SC Trajectory', 'LTM', 'Moon @ CA', 'SC @ CA');

fprintf('\n--- LUNAR ENCOUNTER STATISTICS ---\n');
fprintf('Closest Approach Distance: %.2f km\n', min_dist);
fprintf('Altitude above Surface:    %.2f km\n', min_dist - 1737.4);
fprintf('Time of Closest Approach:  %s (UTC)\n', utc_ca);

% --- FIGURE 3: GROUND TRACK (DOTS) ---
figure('Name', 'Ground Track', 'Color', 'w', 'Position', [100 100 1000 600]);
hold on; grid on;

% Map Outline
try
    load coastlines; 
    plot(coastlon, coastlat, 'k', 'Color', [0.7 0.7 0.7]);
catch
end

% Trajectory Dots
plot(lons, lats, 'b.', 'MarkerSize', 2);

% Stations
st_colors = {'ro', 'mo', 'go', 'ko'}; 
for k = 1:4
    s = const.stations(k);
    plot(s.lon*const.rad2deg, s.lat*const.rad2deg, st_colors{k}, 'MarkerFaceColor', st_colors{k}(1), 'MarkerSize', 8);
    text(s.lon*const.rad2deg+3, s.lat*const.rad2deg+3, s.name, 'FontSize', 10, 'FontWeight', 'bold');
end

% Start & LTM Labels
plot(lons(1), lats(1), 'bs', 'MarkerFaceColor', 'b');
text(lons(1)+3, lats(1), 'Start', 'Color', 'b', 'FontWeight', 'bold');

plot(lons(ltm_idx), lats(ltm_idx), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
text(lons(ltm_idx)+3, lats(ltm_idx), 'LTM', 'Color', 'r', 'FontWeight', 'bold');

xlim([-180 180]); ylim([-90 90]);
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
title('Lunar Trailblazer Ground Track (Reference - Updated Earth Rotation)');