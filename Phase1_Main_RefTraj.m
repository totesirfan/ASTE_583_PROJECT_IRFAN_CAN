% run_reference_trajectory.m
% Generates the nominal trajectory with LTM maneuver and detailed visualization.
% Displays Sun, Earth, and Moon relative geometries.

clear; clc; close all;

%% 1. Setup & Constants
const = get_constants();

% Ensure kernels are loaded
try
    t0 = cspice_str2et(const.epoch_utc_str);
catch
    setup_project();
    t0 = cspice_str2et(const.epoch_utc_str);
end

% Time Definitions
t_LTM = cspice_str2et(const.LTM.date_utc);
tf    = t0 + (251 * const.day2sec);

%% 2. Propagation (Sun-Centered EMO2000)
X0 = [const.X0_ref; const.k_SRP_0; 0];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

fprintf('1/4 Propagating Segment 1 (Detection -> LTM)...\n');
[T1, X1] = ode45(@(t,x) equations_of_motion(t, x, const), [t0 t_LTM], X0, options);

% Apply LTM Maneuver
state_at_LTM = X1(end, :)';
state_at_LTM(4:6) = state_at_LTM(4:6) + const.LTM.dV;
fprintf('    * Applied LTM Delta-V: [%.4f, %.4f, %.4f] km/s\n', const.LTM.dV);

fprintf('2/4 Propagating Segment 2 (LTM -> Flyby)...\n');
[T2, X2] = ode45(@(t,x) equations_of_motion(t, x, const), [t_LTM tf], state_at_LTM, options);

% Combine Trajectory
T_vec = [T1; T2];
X_sc_sun = [X1; X2];
n_steps = length(T_vec);

%% 3. Coordinate Transformations (Post-Processing)
fprintf('3/4 Retrieving Ephemerides (Earth, Moon)...\n');

r_earth_sun = zeros(n_steps, 3);
r_sc_earth  = zeros(n_steps, 3);
r_moon_earth = zeros(n_steps, 3);

for i = 1:n_steps
    t = T_vec(i);
    
    % A. Earth w.r.t Sun (J2000 -> EMO2000)
    st_earth = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_J2000 = st_earth(1:3);
    r_E_Sun_EMO = const.R_EME_EMO * r_E_Sun_J2000;
    r_earth_sun(i, :) = r_E_Sun_EMO';
    
    % B. Spacecraft w.r.t Earth (EMO2000)
    % r_sc_earth = r_sc_sun - r_earth_sun
    r_sc_earth(i, :) = X_sc_sun(i, 1:3) - r_earth_sun(i, :);
    
    % C. Moon w.r.t Earth (J2000 -> EMO2000)
    st_moon = cspice_spkezr('MOON', t, 'J2000', 'NONE', 'EARTH');
    r_moon_J2000 = st_moon(1:3);
    r_moon_EMO = const.R_EME_EMO * r_moon_J2000;
    r_moon_earth(i, :) = r_moon_EMO';
end



%% 4. Visualization
fprintf('4/4 Generating Plots...\n');

% --- FIGURE 1: HELIOCENTRIC VIEW (Sun, Earth Orbit, SC) ---
figure('Name', 'Heliocentric View', 'Color', 'w', 'Position', [100 100 800 600]);
hold on; grid on; axis equal;

% Plot Sun
plot3(0, 0, 0, 'o', 'Color', [1 0.5 0], 'MarkerSize', 12, 'MarkerFaceColor', 'y');

% Plot Earth Orbit
plot3(r_earth_sun(:,1), r_earth_sun(:,2), r_earth_sun(:,3), 'b--', 'LineWidth', 1);
% Plot Earth at start/end
plot3(r_earth_sun(1,1), r_earth_sun(1,2), r_earth_sun(1,3), 'bo', 'MarkerFaceColor', 'b');
text(r_earth_sun(1,1), r_earth_sun(1,2), r_earth_sun(1,3), '  Earth (Start)');

% Plot Spacecraft
plot3(X_sc_sun(:,1), X_sc_sun(:,2), X_sc_sun(:,3), 'm-', 'LineWidth', 1.5);

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Heliocentric Trajectory (Sun-Centered EMO2000)');
legend('Sun', 'Earth Orbit', 'Earth Pos', 'SC Trajectory');
view(2); % Top-down view

% --- FIGURE 2: GEOCENTRIC VIEW (Earth, Moon, Flyby) ---
figure('Name', 'Geocentric View (Moon Flyby)', 'Color', 'w', 'Position', [150 150 800 600]);
hold on; grid on; axis equal;

% Plot Earth
[xE,yE,zE] = sphere(20);
surf(xE*const.R_E, yE*const.R_E, zE*const.R_E, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none');

% Plot Moon Orbit
plot3(r_moon_earth(:,1), r_moon_earth(:,2), r_moon_earth(:,3), 'k:', 'LineWidth', 0.5);
% Plot Moon at Encounter (End)
plot3(r_moon_earth(end,1), r_moon_earth(end,2), r_moon_earth(end,3), 'o', ...
    'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);
text(r_moon_earth(end,1), r_moon_earth(end,2), r_moon_earth(end,3), '  Moon (Encounter)');

% Plot Spacecraft Trajectory
plot3(r_sc_earth(:,1), r_sc_earth(:,2), r_sc_earth(:,3), 'm-', 'LineWidth', 2);

% Plot Sun Direction (Scaled Vector)
sun_vec = -r_earth_sun(1, :) / norm(r_earth_sun(1, :)) * 100000; % 100,000 km length vector
quiver3(0, 0, 0, sun_vec(1), sun_vec(2), sun_vec(3), 'Color', [1 0.5 0], 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(sun_vec(1), sun_vec(2), sun_vec(3), '  To Sun');

% Highlight LTM Point
ltm_idx = length(T1);
plot3(r_sc_earth(ltm_idx,1), r_sc_earth(ltm_idx,2), r_sc_earth(ltm_idx,3), 'r^', 'MarkerFaceColor', 'r');
text(r_sc_earth(ltm_idx,1), r_sc_earth(ltm_idx,2), r_sc_earth(ltm_idx,3), '  LTM Burn');

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Geocentric Trajectory & Lunar Encounter (Earth-Centered EMO2000)');
legend('Earth', 'Moon Orbit', 'Moon', 'SC Trajectory', 'Sun Direction', 'LTM');
%Zoom in on the flyby area if needed
xlim([-4e5 4e5]); ylim([-4e5 4e5]);

% Calculate final distance to Moon center
r_sc_final = r_sc_earth(end, :);
r_moon_final = r_moon_earth(end, :);
dist_to_center = norm(r_sc_final - r_moon_final);
dist_to_surface = dist_to_center - 1737.4; % Moon radius ~1737 km

fprintf('\n--- LUNAR ENCOUNTER STATISTICS ---\n');
fprintf('Distance to Moon Center:  %.2f km\n', dist_to_center);
fprintf('Altitude above Surface:   %.2f km\n', dist_to_surface);

if dist_to_surface > 0
    fprintf('Result: SUCCESSFUL FLYBY (No Impact)\n');
else
    fprintf('Result: IMPACT DETECTED\n');
end
%% 5. Closest Approach Analysis
% Calculate distance to Moon for every time step
dist_to_moon = zeros(n_steps, 1);
for i = 1:n_steps
    % We already have r_sc_earth and r_moon_earth from the previous loop
    dist_to_moon(i) = norm(r_sc_earth(i, :) - r_moon_earth(i, :));
end

% Find the minimum
[min_dist, idx] = min(dist_to_moon);
time_of_closest_approach = T_vec(idx);

% Convert ET to UTC string for readability
utc_ca = cspice_et2utc(time_of_closest_approach, 'C', 3);

fprintf('\n==============================================\n');
fprintf('      LUNAR FLYBY ANALYSIS\n');
fprintf('==============================================\n');
fprintf('Closest Approach Distance: %.2f km\n', min_dist);
fprintf('Altitude above Surface:    %.2f km\n', min_dist - 1737.4);
fprintf('Time of Closest Approach:  %s (UTC)\n', utc_ca);
fprintf('==============================================\n');

% Add the Closest Approach point to your Geocentric Plot
figure(findobj('Name', 'Geocentric View (Moon Flyby)')); % Switch to figure
hold on;
plot3(r_sc_earth(idx,1), r_sc_earth(idx,2), r_sc_earth(idx,3), 'rp', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'm');
text(r_sc_earth(idx,1), r_sc_earth(idx,2), r_sc_earth(idx,3), ...
    sprintf('  Closest Approach: %.0f km', min_dist));
