% Phase1_Main_RefTraj.m
% MASTER SCRIPT: Reference Trajectory, LTM, Flyby, and Plots.
clear; clc; close all;

%% 1. Setup & Propagation
init_project();
const = lib_constants();

% Times & Options
t0    = cspice_str2et(const.epoch_utc_str);
t_LTM = cspice_str2et(const.LTM.date_utc);
tf    = t0 + (251 * const.day2sec);
opts  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

fprintf('1/4 Propagating (Detection -> LTM -> Flyby)...\n');
% Segment 1: Detection -> LTM
[T1, X1] = ode45(@(t,x) lib_dynamics(t, x, const), [t0 t_LTM], [const.X0_ref; 1; 0], opts);

% Apply LTM
X_LTM = X1(end, :)';
X_LTM(4:6) = X_LTM(4:6) + const.LTM.dV;
fprintf('    * LTM Delta-V Applied: [%.4f, %.4f, %.4f] km/s\n', const.LTM.dV);

% Segment 2: LTM -> Flyby
[T2, X2] = ode45(@(t,x) lib_dynamics(t, x, const), [t_LTM tf], X_LTM, opts);
T_vec = [T1; T2]; X_sc_sun = [X1; X2];

%% 3. Coordinate Transformations (Vectorized)
fprintf('3/4 Processing Coordinates...\n');

% Ephemeris (Vectorized) - Note the T_vec' transpose correction
st_earth = cspice_spkezr('EARTH', T_vec', 'J2000', 'NONE', 'SUN');
st_moon  = cspice_spkezr('MOON',  T_vec', 'J2000', 'NONE', 'EARTH');
r_E_Sun_EMO = (const.R_EME_EMO * st_earth(1:3,:))'; 
r_moon_EMO  = (const.R_EME_EMO * st_moon(1:3,:))';

% Relative Vectors
r_sc_earth = X_sc_sun(:, 1:3) - r_E_Sun_EMO;

% Ground Track (ECF Rotation)
r_sc_eci = (const.R_EME_EMO' * r_sc_earth')'; 
phi_G = const.phi_G_detect + const.we * (T_vec - const.t_detect_et);
c = cos(phi_G); s = sin(phi_G);
r_sc_ecf_x =  c .* r_sc_eci(:,1) + s .* r_sc_eci(:,2);
r_sc_ecf_y = -s .* r_sc_eci(:,1) + c .* r_sc_eci(:,2);
[lons, lats, ~] = cart2sph(r_sc_ecf_x, r_sc_ecf_y, r_sc_eci(:,3));
lats = lats * const.rad2deg; lons = lons * const.rad2deg;

% Stats
[min_dist, idx_ca] = min(sqrt(sum((r_sc_earth - r_moon_EMO).^2, 2)));
utc_ca = cspice_et2utc(T_vec(idx_ca), 'C', 3);
ltm_idx = length(T1);

%% 4. Visualization
fprintf('4/4 Generating Plots...\n');

% --- Styles ---
c_traj = [0.85, 0.33, 0.10]; c_earth = [0, 0.45, 0.74]; c_moon = [0.5, 0.5, 0.5];
mk_sty = {'MarkerEdgeColor','k', 'HandleVisibility','off'};
txt_sty = {'FontSize',9,'FontWeight','bold','Color','k','BackgroundColor','w','EdgeColor','k','Margin',1,'VerticalAlignment','bottom','HorizontalAlignment','left'};

% --- FIG 1: HELIOCENTRIC ---
setup_fig('Heliocentric View', 'Heliocentric Trajectory');
rectangle('Position',[-696340,-696340,1.4e6,1.4e6],'Curvature',[1 1],'FaceColor',[1 0.8 0],'EdgeColor','none');
plot(0,0,'o','Color',[1 0.5 0],'DisplayName','Sun');
plot(r_E_Sun_EMO(:,1), r_E_Sun_EMO(:,2), '--', 'Color', c_earth, 'DisplayName', 'Earth Orbit');
plot(X_sc_sun(:,1), X_sc_sun(:,2), '-', 'Color', c_traj, 'LineWidth', 1.5, 'DisplayName', 'SC Trajectory');
legend('Location','northeast');

% --- FIG 2: GEOCENTRIC ---
setup_fig('Geocentric View', 'Geocentric View & Lunar Flyby');
rectangle('Position',[-const.R_E,-const.R_E,2*const.R_E,2*const.R_E],'Curvature',[1 1],'FaceColor',c_earth,'EdgeColor','none');
plot(NaN,NaN,'o','MarkerFaceColor',c_earth,'Color','none','DisplayName','Earth'); 
plot(r_moon_EMO(:,1), r_moon_EMO(:,2), ':', 'Color', c_moon, 'LineWidth', 1, 'DisplayName', 'Moon Orbit');
plot(r_sc_earth(:,1), r_sc_earth(:,2), '-', 'Color', c_traj, 'LineWidth', 1.5, 'DisplayName', 'SC Trajectory');

% Markers (LTM, Moon@CA, SC@CA)
plot(r_sc_earth(ltm_idx,1), r_sc_earth(ltm_idx,2), '^', 'MarkerFaceColor','r', 'Color','k', 'MarkerSize',8, 'DisplayName','LTM');
plot(r_moon_EMO(idx_ca,1), r_moon_EMO(idx_ca,2), 'o', 'MarkerFaceColor',[0.7 0.7 0.7], 'Color','k', 'DisplayName','Moon @ CA');
plot(r_sc_earth(idx_ca,1), r_sc_earth(idx_ca,2), 'p', 'MarkerFaceColor','m', 'Color','k', 'MarkerSize',12, 'DisplayName','SC @ CA');
plot([r_sc_earth(idx_ca,1) r_moon_EMO(idx_ca,1)], [r_sc_earth(idx_ca,2) r_moon_EMO(idx_ca,2)], 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
xlim([-4.5e5, 4.5e5]); ylim([-4.5e5, 4.5e5]); legend('Location','northeast');

% --- FIG 3: GROUND TRACK ---
figure('Name','Ground Track','Color','w','Position',[100 100 1000 600]); ax=axes; hold on;
try load coastlines; plot(coastlon, coastlat, 'k', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility','off'); catch; end

% Trajectory (Blue->Purple->Red)
n_c = 256; cmap_cust = [linspace(0,1,n_c)', zeros(n_c,1), linspace(1,0,n_c)'];
scatter(lons, lats, 6, (T_vec-t0)/const.day2sec, 'filled', 'DisplayName', 'Ground Track');
colormap(ax, cmap_cust); cb=colorbar; cb.Label.String='Days Since Detection';

% Markers Loop (Stations + Events)
% Format: {Lon, Lat, Marker, Color, Label}
markers = {
    const.stations(1).lon, const.stations(1).lat, 's', 'r', 'Goldstone';
    const.stations(2).lon, const.stations(2).lat, 'd', 'm', 'Canberra';
    const.stations(3).lon, const.stations(3).lat, '^', 'g', 'Madrid';
    const.stations(4).lon, const.stations(4).lat, 'p', 'k', 'Antarctica';
    lons(1), lats(1), 's', 'c', 'Start';
    lons(ltm_idx), lats(ltm_idx), '^', 'r', 'LTM';
    lons(idx_ca), lats(idx_ca), 'p', 'm', 'CA'
};

for i=1:size(markers,1)
    % Convert Station Radians to Degrees if needed (Events are already degrees)
    if i <= 4, ln = markers{i,1}*const.rad2deg; lt = markers{i,2}*const.rad2deg; else, ln = markers{i,1}; lt = markers{i,2}; end
    plot(ln, lt, markers{i,3}, 'MarkerFaceColor', markers{i,4}, 'MarkerSize', 9, mk_sty{:});
    text(ln+3, lt+2, markers{i,5}, txt_sty{:});
end

set(gca, 'YDir', 'normal'); xlim([-180 180]); ylim([-90 90]);
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); title('Lunar Trailblazer Ground Track'); grid on;

fprintf('\n--- STATS ---\nCA Dist: %.2f km\nAlt: %.2f km\nTime: %s\n', min_dist, min_dist-1737.4, utc_ca);

% --- Local Helper Function to Replace Anonymous Function ---
function setup_fig(name, title_str)
    figure('Name',name,'Color','w','Position',[50 50 600 500]);
    hold on; grid on; axis equal; view(2);
    title(title_str,'Interpreter','none');
    xlabel('X EMO (km)','Interpreter','none'); 
    ylabel('Y EMO (km)','Interpreter','none');
end
