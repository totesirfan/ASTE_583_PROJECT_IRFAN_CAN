function [Y_comp, H] = lib_measurements(t, X, station_id, const)
    % LIB_MEASUREMENTS Computes G(X) and H-matrix.
    % Frame: Sun-Centered EMO2000
    
    % --- 1. Unpack State (Adaptive) ---
    r_sc = X(1:3);
    v_sc = X(4:6);
    k_SRP = X(7);  % Not used for measurement, but part of state
    bias  = X(8);
    
    % Check if we are estimating Station 4 coordinates
    est_station = (length(X) >= 10) && (station_id == 4);
    
    % --- 2. Get Earth & Rotation Matrices ---
    st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(1:3);
    v_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(4:6);
    
    % GST Calculation
    phi_G = const.phi_G_J2000 + const.we * t;
    c = cos(phi_G); s = sin(phi_G);
    R_ECF_ECI = [c, -s, 0; s, c, 0; 0, 0, 1];
    
    % Combined Rotation (ECF -> EMO)
    % EMO = R_EME_EMO * R_ECI(EME)_ECF
    R_ECF_EMO = const.R_EME_EMO * R_ECF_ECI;
    
    % --- 3. Station State ---
    % If estimating Station 4, use values from State Vector X
    % Otherwise, use values from Constants struct
    if est_station
        lat = X(9);
        lon = X(10);
    else
        stat = const.stations(station_id);
        lat = stat.lat;
        lon = stat.lon;
    end
    
    cl = cos(lat); sl = sin(lat);
    cd = cos(lon); sd = sin(lon);
    
    % Station in ECF
    r_sta_ecf = const.R_E * [cl*cd; cl*sd; sl];
    
    % Velocity in ECF (Earth Rotation Only)
    w_vec_ecf = [0; 0; const.we];
    v_sta_ecf = cross(w_vec_ecf, r_sta_ecf);
    
    % Rotate to EMO2000
    r_sta_emo_earth = R_ECF_EMO * r_sta_ecf;
    v_sta_emo_earth = R_ECF_EMO * v_sta_ecf;
    
    % Shift Origin to Sun
    r_sta = r_E_Sun_EMO + r_sta_emo_earth;
    v_sta = v_E_Sun_EMO + v_sta_emo_earth;
    
    % --- 4. Compute Observations ---
    rho_vec = r_sc - r_sta;
    rho = norm(rho_vec);
    u_rho = rho_vec / rho;
    
    rel_vel = v_sc - v_sta;
    range_rate = dot(rho_vec, rel_vel) / rho;
    
    % Bias Logic (Switch at Day 6)
    t_detect = cspice_str2et('2025 DEC 01 00:00:00.00'); 
    if (t - t_detect) < (6 * 86400)
        y_bias = bias; H_b = 1;
    else
        y_bias = 0; H_b = 0;
    end
    
    Y_comp = [rho; range_rate + y_bias];
    
    % --- 5. H Matrix (Standard) ---
    % Partials w.r.t State [r, v, k, b]
    H_r = u_rho';
    H_v = zeros(1,3);
    H_rr_r = (rel_vel' - range_rate * u_rho') / rho;
    H_rr_v = u_rho';
    
    % Base H matrix (2x8)
    H = [H_r, H_v, 0, 0; 
         H_rr_r, H_rr_v, 0, H_b];
     
    % --- 6. H Matrix (Station Partials) ---
    if est_station
        % We need partials of r_sta_ecf w.r.t lat/lon
        % r_ecf = R [cl*cd; cl*sd; sl]
        
        dr_dlat = const.R_E * [-sl*cd; -sl*sd; cl];
        dr_dlon = const.R_E * [-cl*sd;  cl*cd;  0];
        
        % Rotate derivatives to Inertial Frame (EMO)
        % dr_sta_I / dlat = R_ECF_EMO * dr_dlat
        dr_I_dlat = R_ECF_EMO * dr_dlat;
        dr_I_dlon = R_ECF_EMO * dr_dlon;
        
        % Velocity Partials (Chain rule on v = w x r)
        % dv/dlat = w x (dr/dlat)
        % We need w in Inertial Frame for the cross product with Inertial r? 
        % Easier: transform v_ecf partials.
        % v_ecf = [ -we*y; we*x; 0 ]
        % dv_dlat = [ -we * dy/dlat; we * dx/dlat; 0 ]
        
        % Let's stick to Inertial: v_I = R * (w_ecf x r_ecf)
        % But R is not dependent on lat/lon. w_ecf is const.
        % So dv_I_dlat = R * (w_ecf x dr_ecf_dlat)
        dv_I_dlat = R_ECF_EMO * cross(w_vec_ecf, dr_dlat);
        dv_I_dlon = R_ECF_EMO * cross(w_vec_ecf, dr_dlon);

        % RANGE PARTIALS
        % rho = ||r_sc - r_sta||
        % d_rho/d_sta = -u_rho'
        H_rho_lat = -u_rho' * dr_I_dlat;
        H_rho_lon = -u_rho' * dr_I_dlon;
        
        % RANGE RATE PARTIALS
        % rr = dot(rho_vec, v_sc - v_sta) / rho
        % d_rr / d_r_sta = -H_rr_r (symmetric logic)
        % d_rr / d_v_sta = -H_rr_v = -u_rho'
        
        % Chain Rule:
        % d_rr/dlat = (d_rr/d_r_sta * d_r_sta/dlat) + (d_rr/d_v_sta * d_v_sta/dlat)
        
        d_rr_dr_sta = -H_rr_r; % (1x3)
        d_rr_dv_sta = -u_rho'; % (1x3)
        
        H_rr_lat = d_rr_dr_sta * dr_I_dlat + d_rr_dv_sta * dv_I_dlat;
        H_rr_lon = d_rr_dr_sta * dr_I_dlon + d_rr_dv_sta * dv_I_dlon;
        
        % Append to H Matrix
        H = [H, [H_rho_lat; H_rr_lat], [H_rho_lon; H_rr_lon]];
        
    elseif length(X) >= 10
        % If X has 10 elements but we are not tracking Station 4 (e.g. Goldstone)
        % The partials for lat/lon are 0
        H = [H, zeros(2,2)];
    end
end
