function plot_ground_track()
    % PLOT_GROUND_TRACK Generates the ground track map for the Reference Trajectory
    % Satisfies "Reference Trajectory" requirement 
    
    clear; clc; close all;
    const = get_constants();
    
    % 1. Regenerate Reference Trajectory (Identical to run_reference_trajectory)
    try
        t0 = cspice_str2et(const.epoch_utc_str);
    catch
        setup_project();
        t0 = cspice_str2et(const.epoch_utc_str);
    end
    t_LTM = cspice_str2et(const.LTM.date_utc);
    tf    = t0 + (251 * const.day2sec);
    
    X0 = [const.X0_ref; const.k_SRP_0; 0];
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    
    [T1, X1] = ode45(@(t,x) equations_of_motion(t, x, const), [t0 t_LTM], X0, options);
    
    % Apply LTM
    state_at_LTM = X1(end, :)';
    state_at_LTM(4:6) = state_at_LTM(4:6) + const.LTM.dV;
    
    [T2, X2] = ode45(@(t,x) equations_of_motion(t, x, const), [t_LTM tf], state_at_LTM, options);
    
    % Combine
    T = [T1; T2];
    X_sun = [X1; X2];
    
    % 2. Convert to Body-Fixed (Lat/Lon)
    fprintf('Calculating Ground Track...\n');
    n = length(T);
    lats = zeros(n, 1);
    lons = zeros(n, 1);
    
    for i = 1:n
        t = T(i);
        r_sc_sun = X_sun(i, 1:3)';
        
        % Get Earth State (Sun-Centered J2000)
        st_earth = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
        r_E_Sun_J2000 = st_earth(1:3);
        
        % Rotate Earth State to EMO2000 to match SC
        r_E_Sun_EMO = const.R_EME_EMO * r_E_Sun_J2000;
        
        % SC w.r.t Earth (In EMO2000)
        r_sc_earth_emo = r_sc_sun - r_E_Sun_EMO;
        
        % Rotate EMO2000 -> ECI (J2000)
        % (Inverse of R_EME_EMO is transpose)
        r_sc_eci = const.R_EME_EMO' * r_sc_earth_emo;
        
        % Rotate ECI -> ECF (Body Fixed)
        % GST Calculation
        phi_G_J2000 = 280.46061837 * const.deg2rad;
        phi_G = phi_G_J2000 + const.we * t;
        
        R_ECI_ECF = [cos(phi_G), sin(phi_G), 0;
                    -sin(phi_G), cos(phi_G), 0;
                     0,          0,          1];
        
        r_sc_ecf = R_ECI_ECF * r_sc_eci;
        
        % Cartesian -> Spherical
        [lon, lat, ~] = cart2sph(r_sc_ecf(1), r_sc_ecf(2), r_sc_ecf(3));
        lats(i) = lat * const.rad2deg;
        lons(i) = lon * const.rad2deg;
    end
    
    % 3. Plotting
    figure('Name', 'Ground Track', 'Color', 'w', 'Position', [100 100 1000 600]);
    hold on; grid on;
    
    % Draw Map Outline (Approximate)
    load coastlines; % Built-in MATLAB dataset
    plot(coastlon, coastlat, 'k', 'Color', [0.7 0.7 0.7]);
    
    % Plot Trajectory
    % We use scatter for clarity as lines wrap around 180/-180
    plot(lons, lats, 'b.', 'MarkerSize', 2);
    
    % Plot Stations
    st_colors = {'ro', 'mo', 'go', 'ko'}; % Gold, Can, Mad, Ant
    for k = 1:4
        s = const.stations(k);
        slat = s.lat * const.rad2deg;
        slon = s.lon * const.rad2deg;
        plot(slon, slat, st_colors{k}, 'MarkerFaceColor', st_colors{k}(1), 'MarkerSize', 8);
        text(slon+3, slat+3, s.name, 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    xlim([-180 180]); ylim([-90 90]);
    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
    title('Lunar Trailblazer Ground Track (Post-Detection)');
    
    % Add Start/End labels
    plot(lons(1), lats(1), 'bs', 'MarkerFaceColor', 'b');
    text(lons(1)+3, lats(1), 'Start', 'Color', 'b');
    
    % LTM Location
    ltm_idx = length(T1);
    plot(lons(ltm_idx), lats(ltm_idx), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    text(lons(ltm_idx)+3, lats(ltm_idx), 'LTM', 'Color', 'r');
end