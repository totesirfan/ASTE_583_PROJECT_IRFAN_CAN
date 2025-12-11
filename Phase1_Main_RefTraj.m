% Phase1_Main_RefTraj.m
% Master script: Reference trajectory, LTM, flyby, and plots (incl. 30-min ground track).

clear; clc; close all;

%% User-adjustable time steps
dt_state_sec = 10;    % [s] output step for state history
dt_track_sec = 600;   % [s] ground-track sampling step (30 min)

%% 1. Setup & propagation
init_project();
const = lib_constants();

% Times & ODE options
t0    = cspice_str2et(const.epoch_utc_str);   % detection epoch (ET)
t_LTM = cspice_str2et(const.LTM.date_utc);    % LTM epoch (ET)
tf    = t0 + 251*const.day2sec;               % final time (ET)
opts  = odeset('RelTol',1e-12,'AbsTol',1e-12);

fprintf('1/4 Propagating (Detection -> LTM -> Flyby)...\n');

% Build output grids for ode45
t1_grid = (t0:dt_state_sec:t_LTM).';
if t1_grid(end) < t_LTM
    t1_grid(end+1,1) = t_LTM;
end

t2_grid = (t_LTM:dt_state_sec:tf).';
if t2_grid(end) < tf
    t2_grid(end+1,1) = tf;
end

% Segment 1: detection -> LTM (state: [r; v; k_SRP; bias])
X0_ref = [const.X0_ref; 1; 0];  % [r; v; k_SRP; bias]
[T1, X1] = ode45(@(t,x) lib_dynamics(t,x,const), t1_grid, X0_ref, opts);

% Apply LTM dV at LTM epoch
X_LTM = X1(end,:)';
X_LTM(4:6) = X_LTM(4:6) + const.LTM.dV(:);
fprintf('    * LTM ΔV applied: [%.4f, %.4f, %.4f] km/s\n', const.LTM.dV);

% Segment 2: LTM -> tf
[T2, X2] = ode45(@(t,x) lib_dynamics(t,x,const), t2_grid, X_LTM, opts);

% Concatenate time & state history (Sun-centered EMO2000)
T_vec    = [T1; T2];   % ET times, user-controlled spacing
X_sc_sun = [X1; X2];   % [r_sc_sun; v_sc_sun; params]
ltm_idx  = numel(T1);  % index of LTM in concatenated arrays

%% 2. Coordinate transformations (Earth/Moon, ground track)
fprintf('2/4 Processing coordinates...\n');

% Ephemeris: Earth wrt Sun (J2000, then EMO) & Moon wrt Earth (J2000->EMO)
st_earth = cspice_spkezr('EARTH', T_vec', 'J2000', 'NONE', 'SUN');
st_moon  = cspice_spkezr('MOON',  T_vec', 'J2000', 'NONE', 'EARTH');

r_E_Sun_EMO = (const.R_EME_EMO * st_earth(1:3,:)).';   % km, EMO frame
r_moon_EMO  = (const.R_EME_EMO * st_moon(1:3,:)).';    % km, EMO frame

% SC relative to Earth in EMO
r_sc_sun   = X_sc_sun(:,1:3);
r_sc_earth = r_sc_sun - r_E_Sun_EMO;

% Ground track: EMO -> ECI -> ECF -> lat/lon
r_sc_eci = (const.R_EME_EMO' * r_sc_earth.').';        % EMO->ECI

phi_G = const.phi_G_detect + const.we*(T_vec - const.t_detect_et);  % GMST angle
c = cos(phi_G); s = sin(phi_G);

% Rotate ECI -> ECF (about z-axis)
r_ecf_x =  c.*r_sc_eci(:,1) + s.*r_sc_eci(:,2);
r_ecf_y = -s.*r_sc_eci(:,1) + c.*r_sc_eci(:,2);
r_ecf_z =  r_sc_eci(:,3);

[lons, lats, ~] = cart2sph(r_ecf_x, r_ecf_y, r_ecf_z);
lons = lons*const.rad2deg;
lats = lats*const.rad2deg;

% Moon closest approach (geocentric)
d_vec              = r_sc_earth - r_moon_EMO;
dist               = sqrt(sum(d_vec.^2,2));
[min_dist, idx_ca] = min(dist);
utc_ca             = cspice_et2utc(T_vec(idx_ca), 'C', 3);
dt_ca_days         = (T_vec(idx_ca) - t0)/const.day2sec;
alt_ca             = min_dist - 1737.4;   % km above lunar mean radius

%% 3. Ground-track sampling (using dt_track_sec)
fprintf('3/4 Sampling ground track every %.1f minutes...\n', dt_track_sec/60);

[T_unique, ia] = unique(T_vec,'stable');   % enforce unique sample points
lons_u         = lons(ia);
lats_u         = lats(ia);

T_30    = (T_unique(1):dt_track_sec:T_unique(end)).';    % ET grid
lons_30 = interp1(T_unique, lons_u, T_30, 'linear');     % deg
lats_30 = interp1(T_unique, lats_u, T_30, 'linear');     % deg
days_30 = (T_30 - t0)/const.day2sec;                     % days since detection

%% 4. Visualization
fprintf('4/4 Generating plots...\n');

% --- Color palette (aligned with prelim plots) ---
c_traj       = [0.90 0.40 0.10];   % SC trajectory (warm orange)
c_earth      = [0.00 0.45 0.74];   % Earth / Earth orbit (deep blue)
c_moon       = [0.40 0.40 0.70];   % Moon orbit (slate blue)

c_goldstone  = [0.89 0.10 0.11];   % Station 1
c_canberra   = [0.23 0.49 0.77];   % Station 2
c_madrid     = [0.20 0.63 0.17];   % Station 3
c_antarctica = [0.60 0.31 0.64];   % Station 4

c_start      = [0.00 0.60 0.60];   % Start of ground track (teal)
c_LTM        = [0.80 0.00 0.00];   % LTM marker
c_CA         = [0.70 0.20 0.70];   % CA marker

mk_sty  = {'MarkerEdgeColor','k','HandleVisibility','off'};
txt_sty = {'FontSize',9,'FontWeight','bold','Color','k', ...
           'BackgroundColor','w','EdgeColor','k','Margin',1, ...
           'VerticalAlignment','bottom','HorizontalAlignment','left'};

%% FIG 1 – Heliocentric view
setup_fig('Heliocentric View','Heliocentric Trajectory', ...
          'X_{Sun-EMO} (km)','Y_{Sun-EMO} (km)');

% Sun (to scale)
rectangle('Position',[-696340 -696340 1.4e6 1.4e6], ...
          'Curvature',[1 1],'FaceColor',[1 0.8 0],'EdgeColor','none');
plot(0,0,'o','Color',[1 0.5 0],'DisplayName','Sun');

% Earth orbit + spacecraft trajectory
plot(r_E_Sun_EMO(:,1), r_E_Sun_EMO(:,2),'--', ...
     'Color',c_earth,'DisplayName','Earth Orbit');
plot(r_sc_sun(:,1),   r_sc_sun(:,2),  '-', ...
     'Color',c_traj,'LineWidth',1.5,'DisplayName','SC Trajectory');

% Heliocentric LTM & CA
plot(r_sc_sun(ltm_idx,1), r_sc_sun(ltm_idx,2), '^', ...
     'MarkerFaceColor',c_LTM, ...
     'MarkerEdgeColor','k', ...
     'MarkerSize',8, ...
     'DisplayName','LTM');

plot(r_sc_sun(idx_ca,1),  r_sc_sun(idx_ca,2),  'p', ...
     'MarkerFaceColor',c_CA, ...
     'MarkerEdgeColor',[0.3 0 0.3], ...
     'MarkerSize',10, ...
     'DisplayName','CA');

legend('Location','northeast');

%% FIG 2 – Geocentric + lunar flyby
setup_fig('Geocentric View','Geocentric View & Lunar Flyby', ...
          'X_{Earth-EMO} (km)','Y_{Earth-EMO} (km)');

% Earth to scale
rectangle('Position',[-const.R_E -const.R_E 2*const.R_E 2*const.R_E], ...
          'Curvature',[1 1],'FaceColor',c_earth,'EdgeColor','none');
plot(NaN,NaN,'o','MarkerFaceColor',c_earth,'Color','none', ...
     'DisplayName','Earth');

% Moon orbit + spacecraft trajectory
plot(r_moon_EMO(:,1), r_moon_EMO(:,2),':', ...
     'Color',c_moon,'LineWidth',1,'DisplayName','Moon Orbit');
plot(r_sc_earth(:,1), r_sc_earth(:,2),'-', ...
     'Color',c_traj,'LineWidth',1.5,'DisplayName','SC Trajectory');

% Moon @ CA, SC @ CA (geocentric)
plot(r_moon_EMO(idx_ca,1),  r_moon_EMO(idx_ca,2),  'o', ...
     'MarkerFaceColor',[0.8 0.8 0.8], ...
     'MarkerEdgeColor','k', ...
     'MarkerSize',7, ...
     'DisplayName','Moon @ CA');

plot(r_sc_earth(idx_ca,1),  r_sc_earth(idx_ca,2),  'p', ...
     'MarkerFaceColor',c_CA, ...
     'MarkerEdgeColor',[0.3 0 0.3], ...
     'LineStyle','none', ...
     'LineWidth',0.5, ...
     'MarkerSize',7, ...
     'DisplayName','SC @ CA');

% Line from SC to Moon at CA
plot([r_sc_earth(idx_ca,1) r_moon_EMO(idx_ca,1)], ...
     [r_sc_earth(idx_ca,2) r_moon_EMO(idx_ca,2)], ...
     'k-','LineWidth',0.5,'HandleVisibility','off');

xlim([-4.5e5 4.5e5]);
ylim([-4.5e5 4.5e5]);
legend('Location','northeast');

%% FIG 3 – Ground track (dt_track_sec samples)
figure('Name','Ground Track','Color','w','Position',[100 100 1000 600]);
ax = axes; hold on;

% Map outline
try
    load coastlines;
    plot(coastlon, coastlat,'Color',[0.6 0.6 0.6], ...
         'LineWidth',0.5,'HandleVisibility','off');
catch
end

% Color map for time progression
n_c = 256;
cmap_cust = [linspace(0,1,n_c)', zeros(n_c,1), linspace(1,0,n_c)'];

% Ground-track points
scatter(lons_30, lats_30, 1, days_30, 'filled', 'DisplayName','Ground Track');
colormap(ax,cmap_cust);
cb = colorbar; 
cb.Label.String = 'Days Since Detection';

% Stations + trajectory events
markers = {
    const.stations(1).lon, const.stations(1).lat, 's', c_goldstone,  'Goldstone';
    const.stations(2).lon, const.stations(2).lat, 'd', c_canberra,   'Canberra';
    const.stations(3).lon, const.stations(3).lat, '^', c_madrid,     'Madrid';
    const.stations(4).lon, const.stations(4).lat, 'p', c_antarctica, 'Antarctica';
    lons(1),        lats(1),        's', c_start,  'Start';
    lons(ltm_idx),  lats(ltm_idx),  '^', c_LTM,    'LTM';
    lons(idx_ca),   lats(idx_ca),   'p', c_CA,     'CA';
};

for i = 1:size(markers,1)
    if i <= 4
        ln = markers{i,1}*const.rad2deg;
        lt = markers{i,2}*const.rad2deg;
    else
        ln = markers{i,1};
        lt = markers{i,2};
    end
    plot(ln, lt, markers{i,3}, ...
         'MarkerFaceColor',markers{i,4}, ...
         'MarkerEdgeColor','k', ...
         'MarkerSize',9, ...
         mk_sty{:});
    text(ln+3, lt+2, markers{i,5}, txt_sty{:});
end

set(gca,'YDir','normal');
xlim([-180 180]); ylim([-90 90]);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title(sprintf('Lunar Trailblazer Ground Track (%.1f min Samples)', dt_track_sec/60));
grid on;

%% FIG 4 – 3D Moon flyby at closest approach (textured Moon)
figure('Name','3D Moon Flyby at CA','Color','w','Position',[200 100 800 700]);
hold on; grid on; axis equal; view(3);

% Use Moon-centered coordinates at CA for clean local geometry
rM_ca      = r_moon_EMO(idx_ca,:);      % Moon position at CA (Earth-EMO)
r_moon_MC  = r_moon_EMO - rM_ca;        % Moon trajectory, Moon-centered
r_sc_MC    = r_sc_earth - rM_ca;        % SC trajectory, Moon-centered
sc_ca_MC   = r_sc_MC(idx_ca,:);         % SC at CA in Moon-centered frame

% Subset near CA (e.g., where SC-Moon distance < 20 R_moon)
R_moon = 1737.4;    % km
if isfield(const,'R_Moon')
    R_moon = const.R_Moon;
end
idx_win = dist < 20*R_moon;
if ~any(idx_win)
    idx_win = true(size(dist));   % fallback: plot all
end

% Plot Moon trajectory and spacecraft trajectory (near CA)
plot3(r_moon_MC(idx_win,1), r_moon_MC(idx_win,2), r_moon_MC(idx_win,3), ...
      '-','Color',c_moon,'LineWidth',1.2,'DisplayName','Moon Trajectory');
plot3(r_sc_MC(idx_win,1),   r_sc_MC(idx_win,2),   r_sc_MC(idx_win,3), ...
      '-','Color',c_traj,'LineWidth',1.5,'DisplayName','SC Trajectory');

% Textured Moon model at origin of this frame
[XS, YS, ZS] = sphere(256);
XS = R_moon * XS;
YS = R_moon * YS;
ZS = R_moon * ZS;

has_tex = false;
try
    moon_tex = imread('lroc_color_poles_2k.tif');
    has_tex  = true;
catch
    warning('Moon texture file lroc_color_poles_2k.tif not found. Using gray sphere.');
end

if has_tex
    surface(XS, YS, ZS, ...
            'FaceColor','texturemap', ...
            'CData',moon_tex, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');
else
    surf(XS, YS, ZS, ...
         'FaceColor',[0.8 0.8 0.8], ...
         'EdgeColor','none', ...
         'HandleVisibility','off');
end

shading interp;
camlight headlight;
lighting gouraud;

% SC @ CA marker (smaller)
plot3(sc_ca_MC(1), sc_ca_MC(2), sc_ca_MC(3), 'p', ...
      'MarkerFaceColor',c_CA, ...
      'MarkerEdgeColor',[0.3 0 0.3], ...
      'MarkerSize',5, ...
      'DisplayName','SC @ CA');

% CA altitude annotation – close to SC, offset along radial direction
r_hat = sc_ca_MC / norm(sc_ca_MC);
text_pos = sc_ca_MC + 0.3*R_moon * r_hat;   % small radial offset
text(text_pos(1), text_pos(2), text_pos(3), ...
     sprintf('CA Altitude = %.1f km', alt_ca), ...
     'FontSize',10, ...
     'FontWeight','bold', ...
     'HorizontalAlignment','center');

% Zoom around the Moon (~±4000 km box)
lim = 2500;   % km
xlim([-lim lim]);
ylim([-lim lim]);
zlim([-lim lim]);

xlabel('X_{Moon-EMO, CA} (km)','Interpreter','tex');
ylabel('Y_{Moon-EMO, CA} (km)','Interpreter','tex');
zlabel('Z_{Moon-EMO, CA} (km)','Interpreter','tex');

title(sprintf(['Lunar Flyby Geometry at Closest Approach\n', ...
               't = %.3f Days from t0, %s UTC'], ...
               dt_ca_days, utc_ca), ...
      'Interpreter','none');

legend('Location','bestoutside');


%% Stats
fprintf('\n--- STATS ---\n');
fprintf('Closest approach distance: %.2f km\n', min_dist);
fprintf('Closest approach altitude: %.2f km\n', alt_ca);
fprintf('Closest approach UTC     : %s\n', utc_ca);

%% Local helper: figure setup
function setup_fig(name, title_str, xlab, ylab)
    figure('Name',name,'Color','w','Position',[50 50 600 500]);
    hold on; grid on; axis equal; view(2);
    title(title_str,'Interpreter','none');
    xlabel(xlab,'Interpreter','tex');
    ylabel(ylab,'Interpreter','tex');
end
