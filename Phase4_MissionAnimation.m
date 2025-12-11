function Phase4_MissionAnimation()
% PHASE4_MISSIONANIMATION
% - 60 FPS: Ultra-smooth playback.
% - COMPRESSION: High quality MPEG-4 (Quality 100).
% - DYNAMICS: Stride adjusted for fluid 60fps timing.
% - VISUALS: Improved 3D plot aesthetics and robust fixed legends.

    clear; clc; close all;
    fprintf('=== GENERATING 60FPS MISSION DASHBOARD ===\n');

    %% 0. USER SETTINGS
    save_video = false;      
    video_filename = 'Lunar_Trailblazer_60fps.mp4';
    target_fps = 60;        % <--- HIGH FRAME RATE (Standard Smooth)
    % target_fps = 120;     % <--- OPTIONAL: Very High Frame Rate (Large File)
    video_quality = 50;    % Max Quality

    %% 1. Setup & Constants
    try, init_project(); catch, warning('Run init_project first!'); end
    const = lib_constants();
    R_Moon = 1737.4; % km
    
    % Colors
    c_back   = 'k'; c_axis = 'w';              
    c_traj   = [0.90 0.40 0.10]; % Orange
    c_earth  = [0.00 0.45 0.74]; % Earth Blue
    c_moon   = [0.70 0.70 0.70]; % Moon Grey
    c_sun    = [1.00 0.84 0.00]; % Sun Gold
    c_LTM    = [0.80 0.00 0.00]; c_LCM    = [0.00 0.80 0.80]; c_CA = [0.70 0.20 0.70];
    c_gold   = [0.89 0.10 0.11]; c_can    = [0.23 0.49 0.77];
    c_mad    = [0.20 0.63 0.17]; c_ant    = [0.60 0.31 0.64];
    
    % Texture
    moon_tex = cat(3, uint8([150 150; 150 150]), uint8([150 150; 150 150]), uint8([150 150; 150 150]));
    tex_filename = 'lroc_color_poles_2k.tif'; 
    if exist(tex_filename, 'file')
        try
            [img, map] = imread(tex_filename);
            if ~isempty(map), img = ind2rgb(img, map); end
            if size(img, 3) == 1, img = repmat(img, [1 1 3]); end
            moon_tex = img;
        catch, end
    end

    %% 2. Propagation
    t0 = cspice_str2et(const.epoch_utc_str);
    t_LTM = cspice_str2et(const.LTM.date_utc);
    t_LCM = cspice_str2et(const.LCM.date_utc);
    tf = t0 + 251*const.day2sec; 
    
    dt_step = 60; 
    t_vec = (t0 : dt_step : tf).';
    
    fprintf('Propagating trajectory...\n');
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    X0_ref = [const.X0_ref; 1; 0];
    
    t1_vec = t_vec(t_vec <= t_LTM); if isempty(t1_vec), t1_vec = [t0; t_LTM]; end
    [~, X1] = ode45(@(t,x) lib_dynamics(t,x,const), t1_vec, X0_ref, opts);
    X_LTM_State = X1(end,:)'; X_LTM_State(4:6) = X_LTM_State(4:6) + const.LTM.dV(:);
    t2_vec = t_vec(t_vec > t_LTM);
    [~, X2] = ode45(@(t,x) lib_dynamics(t,x,const), [t1_vec(end); t2_vec], X_LTM_State, opts);
    X2 = X2(2:end, :); X_all = [X1; X2];
    
    %% 3. Coordinate Calculations
    fprintf('Calculating frames...\n');
    st_earth = cspice_spkezr('EARTH', t_vec', 'J2000', 'NONE', 'SUN');
    r_earth_sun  = (const.R_EME_EMO * st_earth(1:3,:)).'; 
    st_moon_earth  = cspice_spkezr('MOON',  t_vec', 'J2000', 'NONE', 'EARTH'); 
    r_moon_earth  = (const.R_EME_EMO * st_moon_earth(1:3,:)).';
    
    r_sc_sun_raw = X_all(:,1:3);
    r_sc_earth = r_sc_sun_raw - r_earth_sun;
    r_moon_sun = r_earth_sun + r_moon_earth;
    r_sc_moon  = r_sc_earth - r_moon_earth;  
    
    dist_sc_moon = sqrt(sum(r_sc_moon.^2, 2));
    [min_dist, idx_CA] = min(dist_sc_moon);
    ca_altitude = min_dist - R_Moon; 
    fprintf('  * CA Altitude: %.2f km\n', ca_altitude);
    
    [~, idx_LTM] = min(abs(t_vec - t_LTM)); 
    [~, idx_LCM] = min(abs(t_vec - t_LCM));
    
    r_sc_eci = (const.R_EME_EMO' * r_sc_earth.').';
    phi_G = const.phi_G_detect + const.we*(t_vec - t0);
    c = cos(phi_G); s = sin(phi_G);
    x = c.*r_sc_eci(:,1) + s.*r_sc_eci(:,2);
    y = -s.*r_sc_eci(:,1) + c.*r_sc_eci(:,2);
    z = r_sc_eci(:,3);
    [lon_rad, lat_rad, ~] = cart2sph(x, y, z);
    lon_deg = lon_rad * const.rad2deg;
    lat_deg = lat_rad * const.rad2deg;

    %% 4. Setup Figure (STABLE 1080p)
    fig = figure('Name','Mission Dashboard 3D', ...
                 'Color',c_back, ...
                 'Position', [50, 50, 1920, 1080], ...  
                 'InvertHardcopy','off');
    
    drawnow; 
    pause(1.0); % Stabilization pause

    if save_video
        fprintf('Setting up Video Writer (%s)...\n', video_filename);
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = target_fps;
        v.Quality = video_quality; 
        open(v);
    end

    % --- LEFT: 3D TRAJECTORY ---
    ax1 = axes('Position', [0.05 0.25 0.40 0.65], 'Color',c_back, 'XColor',c_axis, 'YColor',c_axis, 'ZColor',c_axis); 
    hold on; grid on; axis equal; view(3);
    ax1.GridColor = [0.3 0.3 0.3]; ax1.GridAlpha = 0.5;
    
    % --- LEGEND FIX: Create Dummy Handles for Legend ---
    % These plots are invisible (NaN data) but define the legend entry styles
    L_sun   = plot3(nan,nan,nan, 'o', 'Color', c_sun, 'MarkerFaceColor', c_sun, 'MarkerSize', 8, 'DisplayName', 'Sun');
    L_earth = plot3(nan,nan,nan, 'o', 'Color', c_earth, 'MarkerFaceColor', c_earth, 'MarkerSize', 6, 'DisplayName', 'Earth');
    L_moon  = plot3(nan,nan,nan, 'o', 'Color', 'none', 'MarkerFaceColor', c_moon, 'MarkerEdgeColor', 'w', 'MarkerSize', 6, 'DisplayName', 'Moon');
    L_sc    = plot3(nan,nan,nan, '-', 'Color', c_traj, 'LineWidth', 1.5, 'DisplayName', 'Spacecraft');
    L_LTM   = plot3(nan,nan,nan, '^', 'MarkerFaceColor', c_LTM, 'MarkerEdgeColor', 'w', 'DisplayName', 'LTM');
    L_LCM   = plot3(nan,nan,nan, 'd', 'MarkerFaceColor', c_LCM, 'MarkerEdgeColor', 'k', 'DisplayName', 'LCM');
    
    % Add Legend to 3D Plot
    legend([L_sun, L_earth, L_moon, L_sc, L_LTM, L_LCM], 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none', 'Location', 'northeast');

    % Actual Objects (Dynamic)
    h_sun      = plot3(0,0,0, 'o', 'Color',c_sun, 'MarkerFaceColor',c_sun, 'MarkerSize',8, 'HandleVisibility','off'); 
    h_earth_mk = plot3(0,0,0, 'o', 'Color',c_earth, 'MarkerFaceColor',c_earth, 'MarkerSize',6, 'HandleVisibility','off');
    h_trail_sc_helio    = plot3(nan,nan,nan, '-', 'Color',c_traj, 'LineWidth',1.5, 'HandleVisibility','off');
    h_trail_earth_helio = plot3(nan,nan,nan, '--', 'Color',c_earth, 'LineWidth',1, 'HandleVisibility','off');
    h_trail_sc_zoom     = plot3(nan,nan,nan, '-', 'Color',c_traj, 'LineWidth',2.0, 'Visible', 'off', 'HandleVisibility','off');
    
    [sx,sy,sz] = sphere(50); 
    h_moon_body = surface(sx*R_Moon, sy*R_Moon, sz*R_Moon, ...
        'FaceColor','texturemap', 'CData', moon_tex, ...
        'EdgeColor','none', 'Visible','off', ...
        'FaceLighting','gouraud', 'SpecularStrength',0.1, 'AmbientStrength', 0.8, 'HandleVisibility','off'); 
    
    h_sc_3d    = plot3(0,0,0, 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor','none', 'MarkerSize',6, 'HandleVisibility','off');
    h_moon_mk  = plot3(0,0,0, 'o', 'MarkerFaceColor',c_moon, 'MarkerEdgeColor','none', 'MarkerSize',4, 'HandleVisibility','off'); 
    h_mk_LTM = plot3(0,0,0, '^', 'MarkerFaceColor',c_LTM, 'MarkerEdgeColor','w', 'MarkerSize',8, 'Visible','off', 'HandleVisibility','off');
    h_mk_LCM = plot3(0,0,0, 'd', 'MarkerFaceColor',c_LCM, 'MarkerEdgeColor','k', 'MarkerSize',8, 'Visible','off', 'HandleVisibility','off');
    
    h_light = camlight(ax1, 'headlight'); 
    material dull; 
    title('Heliocentric Transfer','Color',c_axis);

    % --- RIGHT: Ground Track ---
    ax2 = axes('Position', [0.52 0.25 0.45 0.65], 'Color',c_back, 'XColor',c_axis, 'YColor',c_axis); 
    hold on; grid on; axis equal;
    ax2.GridColor = [0.3 0.3 0.3]; ax2.GridAlpha = 0.5;
    try load coastlines; plot(coastlon, coastlat, 'Color',[0.5 0.5 0.5]); catch; end
    
    st_colors = {c_gold, c_can, c_mad, c_ant};
    for i=1:4
        plot(const.stations(i).lon*const.rad2deg, const.stations(i).lat*const.rad2deg, ...
            's', 'MarkerFaceColor', st_colors{i}, 'MarkerEdgeColor','w', 'MarkerSize',8);
        text(const.stations(i).lon*const.rad2deg+3, const.stations(i).lat*const.rad2deg+3, ...
            const.stations(i).name, 'FontSize',9, 'Color',c_axis, 'Interpreter','none');
    end
    
    h_map_LTM = plot(lon_deg(idx_LTM), lat_deg(idx_LTM), '^', 'MarkerFaceColor',c_LTM, 'MarkerEdgeColor','w', 'MarkerSize',10, 'Visible','off');
    h_txt_LTM = text(lon_deg(idx_LTM)+3, lat_deg(idx_LTM), 'LTM', 'Color',c_LTM, 'FontSize',10, 'FontWeight','bold', 'Visible','off');
    h_map_LCM = plot(lon_deg(idx_LCM), lat_deg(idx_LCM), 'd', 'MarkerFaceColor',c_LCM, 'MarkerEdgeColor','k', 'MarkerSize',10, 'Visible','off');
    h_txt_LCM = text(lon_deg(idx_LCM)+3, lat_deg(idx_LCM), 'LCM', 'Color',c_LCM, 'FontSize',10, 'FontWeight','bold', 'Visible','off');
    
    h_trail_ground = plot(nan,nan, '-', 'Color',c_traj, 'LineWidth',1.5);
    h_sc_ground = plot(nan,nan, 'p', 'MarkerFaceColor',c_traj, 'MarkerEdgeColor','w', 'MarkerSize',12);
    
    xlim([-180 180]); ylim([-90 90]);
    title('Ground Track (Last 6 Hours)','Color',c_axis);

    % Main Title
    htitle = sgtitle('Init...', 'Color',c_axis, 'FontSize',14, 'FontWeight','bold');

    %% 5. Animation Loop
    fprintf('Starting animation... (Recording: %s)\n', video_filename);
    n_total = length(t_vec);
    mode = 'HELIO'; 
    k = 1;
    tail_length = 180; 

    while k <= n_total
        if ~isvalid(fig), break; end
        
        curr_dist = dist_sc_moon(k);
        curr_t = t_vec(k);
        
        % 60 FPS SCALING:
        if curr_dist > 1e6
            stride = 250; 
        else
            stride = max(1, ceil(250 * (curr_dist / 1e6)));
        end
        
        slice_helio = 1:stride:k; 
        start_idx = max(1, k - tail_length);
        slice_ground = start_idx:ceil(stride/5):k; 
        
        if curr_dist < 800000 
            % ZOOM MODE
            if strcmp(mode, 'HELIO')
                mode = 'MOON';
                set([h_sun, h_earth_mk, h_trail_sc_helio, h_trail_earth_helio, h_moon_mk], 'Visible', 'off');
                set([h_trail_sc_zoom, h_moon_body], 'Visible', 'on');
                set(h_moon_body, 'XData', sx*R_Moon, 'YData', sy*R_Moon, 'ZData', sz*R_Moon);
            end
            
            slice_zoom = max(1, k-2000):stride:k; 
            set(h_trail_sc_zoom, 'XData', r_sc_moon(slice_zoom,1), 'YData', r_sc_moon(slice_zoom,2), 'ZData', r_sc_moon(slice_zoom,3));
            sc_pos = r_sc_moon(k, :);
            set(h_sc_3d, 'XData', sc_pos(1), 'YData', sc_pos(2), 'ZData', sc_pos(3));
            
            win = max(4000, 1.5 * curr_dist); 
            xlim(ax1, [-win, win]); ylim(ax1, [-win, win]); zlim(ax1, [-win, win]);
            set(h_sc_3d, 'MarkerSize', max(2, 6 * (curr_dist / 100000))); 
            camlight(h_light, 'headlight'); 
            
            title(ax1, sprintf('FLYBY: Altitude %.1f km', curr_dist - R_Moon), 'Color', c_CA);
            
        else
            % CRUISE MODE
            if strcmp(mode, 'MOON')
                mode = 'HELIO';
                set([h_sun, h_earth_mk, h_trail_sc_helio, h_trail_earth_helio, h_moon_mk], 'Visible', 'on');
                set([h_trail_sc_zoom, h_moon_body], 'Visible', 'off');
                xlim(ax1, 'auto'); ylim(ax1, 'auto'); zlim(ax1, 'auto');
                set(h_sc_3d, 'MarkerSize', 6);
                title(ax1, 'Heliocentric Transfer', 'Color', c_axis);
            end
            
            set(h_trail_sc_helio, 'XData', r_sc_sun_raw(slice_helio,1), 'YData', r_sc_sun_raw(slice_helio,2), 'ZData', r_sc_sun_raw(slice_helio,3));
            set(h_trail_earth_helio, 'XData', r_earth_sun(slice_helio,1), 'YData', r_earth_sun(slice_helio,2), 'ZData', r_earth_sun(slice_helio,3));
            set(h_sc_3d, 'XData', r_sc_sun_raw(k,1), 'YData', r_sc_sun_raw(k,2), 'ZData', r_sc_sun_raw(k,3));
            set(h_earth_mk, 'XData', r_earth_sun(k,1), 'YData', r_earth_sun(k,2), 'ZData', r_earth_sun(k,3));
            set(h_moon_mk, 'XData', r_moon_sun(k,1), 'YData', r_moon_sun(k,2), 'ZData', r_moon_sun(k,3));
            set(h_mk_LTM, 'XData', r_sc_sun_raw(idx_LTM,1), 'YData', r_sc_sun_raw(idx_LTM,2), 'ZData', r_sc_sun_raw(idx_LTM,3));
            set(h_mk_LCM, 'XData', r_sc_sun_raw(idx_LCM,1), 'YData', r_sc_sun_raw(idx_LCM,2), 'ZData', r_sc_sun_raw(idx_LCM,3));
        end
        
        set(h_trail_ground, 'XData', lon_deg(slice_ground), 'YData', lat_deg(slice_ground));
        set(h_sc_ground, 'XData', lon_deg(k), 'YData', lat_deg(k));
        
        if curr_t >= t_LTM, set([h_mk_LTM, h_map_LTM, h_txt_LTM], 'Visible', 'on'); end
        if curr_t >= t_LCM, set([h_mk_LCM, h_map_LCM, h_txt_LCM], 'Visible', 'on'); end
        
        days = (curr_t - t0)/const.day2sec;
        set(htitle, 'String', sprintf('MET: %.1f days | Dist to Moon: %.0f km', days, curr_dist));
        
        drawnow;
        
        if save_video
            frame = getframe(fig);
            writeVideo(v, frame);
        end
        
        if k < idx_CA && (k + stride) >= idx_CA
            k = idx_CA;
            sc_pos = r_sc_moon(k, :);
            set(h_sc_3d, 'XData', sc_pos(1), 'YData', sc_pos(2), 'ZData', sc_pos(3));
            set(h_moon_body, 'Visible', 'on');
            camlight(h_light, 'headlight');
            
            % FINAL TITLE UPDATE
            title(ax1, sprintf('CLOSEST APPROACH: Altitude %.1f km', min_dist - R_Moon), 'Color', c_CA);
            
            drawnow;
            if save_video
                % Freeze on final frame for 2 seconds (120 frames)
                for z = 1:120, writeVideo(v, getframe(fig)); end
            end
            fprintf('Flyby Reached.\n');
            break; 
        end
        
        k = k + stride;
    end
    
    if save_video
        close(v);
        fprintf('Video saved to %s\n', video_filename);
    end
end