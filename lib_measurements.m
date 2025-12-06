function [Y_comp, H] = lib_measurements(t, X, station_id, const)
    % LIB_MEASUREMENTS Computes computed observations G(X) and the H-matrix.
    % Frame: Sun-Centered EMO2000
    % Updates: 
    %   - Uses accurate GST Epoch from lib_constants (Handout aligned)
    %   - Handles Doppler Bias for t < 6 days
    %   - Supports 7-state (Prelim) and 10-state (Final) estimation
    
    % --- 1. Unpack State ---
    r_sc = X(1:3);
    v_sc = X(4:6);
    
    % Handle variable state sizes (7, 8, or 10 elements)
    % k_SRP is usually index 7
    if length(X) >= 7
        k_SRP = X(7); 
    else
        k_SRP = 1.0; 
    end
    
    % Bias is usually index 8
    if length(X) >= 8
        bias = X(8); 
    else
        bias = 0.0; 
    end
    
    % Detect if estimating Station 4 coordinates (State size 10)
    % If X has 10 elements AND we are processing Station 4, use the estimated Lat/Lon
    est_station = (length(X) >= 10) && (station_id == 4);
    
    % --- 2. Ephemeris & Rotation Matrices ---
    st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(1:3);
    v_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(4:6);
    
    % --- GST Calculation (Handout Specific) ---
    % We use the exact epoch and angle defined in lib_constants
    % phi_G(t) = phi_G0 + we * (t - t0)
    phi_G = const.phi_G_detect + const.we * (t - const.t_detect_et);
    
    c = cos(phi_G); s = sin(phi_G);
    
    % Rotation ECF -> ECI (Z-axis rotation)
    R_ECF_ECI = [c, -s, 0; 
                 s,  c, 0; 
                 0,  0, 1];
    
    % Rotation ECF -> EMO (Inertial Integration Frame)
    R_ECF_EMO = const.R_EME_EMO * R_ECF_ECI;
    
    % --- 3. Station State ---
    if est_station
        % Use Estimated Lat/Lon from State Vector
        lat = X(9);
        lon = X(10);
    else
        % Use Nominal Lat/Lon from Constants
        stat = const.stations(station_id);
        lat = stat.lat;
        lon = stat.lon;
    end
    
    cl = cos(lat); sl = sin(lat);
    cd = cos(lon); sd = sin(lon);
    
    % Station Position (ECF)
    r_sta_ecf = const.R_E * [cl*cd; cl*sd; sl];
    
    % Station Velocity (ECF) - Due to Earth Rotation
    w_vec_ecf = [0; 0; const.we];
    v_sta_ecf = cross(w_vec_ecf, r_sta_ecf);
    
    % Transform to Inertial Sun-Centered Frame
    r_sta_emo_earth = R_ECF_EMO * r_sta_ecf;
    v_sta_emo_earth = R_ECF_EMO * v_sta_ecf;
    
    r_sta = r_E_Sun_EMO + r_sta_emo_earth;
    v_sta = v_E_Sun_EMO + v_sta_emo_earth;
    
    % --- 4. Compute Observations ---
    rho_vec = r_sc - r_sta;
    rho = norm(rho_vec);
    u_rho = rho_vec / rho;
    
    rel_vel = v_sc - v_sta;
    range_rate = dot(rho_vec, rel_vel) / rho;
    
    % Apply Range Rate Bias (First 6 days only)
    % DCO for Bias is 6 days after detection
    t_bias_end = const.t_detect_et + (6 * 86400);
    
    if t < t_bias_end
        y_bias = bias; 
        H_b = 1; % Partial derivative d(Doppler)/d(bias)
    else
        y_bias = 0; 
        H_b = 0;
    end
    
    Y_comp = [rho; range_rate + y_bias];
    
    % --- 5. H Matrix (State Partials) ---
    % Partials w.r.t Spacecraft State [r, v]
    H_r = u_rho';
    H_v = zeros(1,3);
    
    % Doppler Partials
    % d(rho_dot)/dr = (v_rel' - rho_dot*u_rho') / rho
    H_rr_r = (rel_vel' - range_rate * u_rho') / rho;
    % d(rho_dot)/dv = u_rho'
    H_rr_v = u_rho';
    
    % Base H matrix [r(3), v(3), k_SRP(1), bias(1)]
    % Note: k_SRP partial is 0 for measurement geometry
    H = [H_r,    H_v,    0, 0; 
         H_rr_r, H_rr_v, 0, H_b];
     
    % --- 6. H Matrix (Station Partials) ---
    if est_station
        % Derivative of Station Position (ECF) w.r.t Lat/Lon
        dr_dlat = const.R_E * [-sl*cd; -sl*sd; cl];
        dr_dlon = const.R_E * [-cl*sd;  cl*cd;  0];
        
        % Transform Partials to Inertial Frame
        dr_I_dlat = R_ECF_EMO * dr_dlat;
        dr_I_dlon = R_ECF_EMO * dr_dlon;
        
        % Derivative of Station Velocity (Inertial) w.r.t Lat/Lon
        % v_I = R * (w x r_ecf) -> dv_I = R * (w x dr_ecf)
        dv_I_dlat = R_ECF_EMO * cross(w_vec_ecf, dr_dlat);
        dv_I_dlon = R_ECF_EMO * cross(w_vec_ecf, dr_dlon);

        % Range Partials (Chain Rule)
        % d(rho)/d(sta) = -u_rho
        H_rho_lat = -u_rho' * dr_I_dlat;
        H_rho_lon = -u_rho' * dr_I_dlon;
        
        % Range Rate Partials (Chain Rule)
        d_rr_dr_sta = -H_rr_r; % Partial w.r.t station pos is negative of partial w.r.t satellite pos
        d_rr_dv_sta = -u_rho'; % Partial w.r.t station vel is negative of partial w.r.t satellite vel
        
        H_rr_lat = d_rr_dr_sta * dr_I_dlat + d_rr_dv_sta * dv_I_dlat;
        H_rr_lon = d_rr_dr_sta * dr_I_dlon + d_rr_dv_sta * dv_I_dlon;
        
        % Append Station Partials to H (Cols 9 and 10)
        H = [H, [H_rho_lat; H_rr_lat], [H_rho_lon; H_rr_lon]];
        
    elseif length(X) >= 10
        % Padding for non-estimated stations (so H remains 2x10)
        H = [H, zeros(2,2)];
    end
end