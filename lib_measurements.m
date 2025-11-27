function [Y_comp, H] = lib_measurements(t, X, station_id, const)
    % LIB_MEASUREMENTS Computes G(X) and H-matrix.
    % Frame: Sun-Centered EMO2000
    
    r_sc = X(1:3);
    v_sc = X(4:6);
    bias = X(8);
    
    % 1. Get Earth State (Rotate J2000 -> EMO2000)
    st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(1:3);
    v_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(4:6);
    
    % 2. Station State (Chain: ECF -> ECI -> EMO)
    % GST relative to J2000
    phi_G = const.phi_G_J2000 + const.we * t;
    c = cos(phi_G); s = sin(phi_G);
    R_ECF_ECI = [c, -s, 0; s, c, 0; 0, 0, 1];
    
    % Station in ECF
    stat = const.stations(station_id);
    cl = cos(stat.lat); sl = sin(stat.lat);
    cd = cos(stat.lon); sd = sin(stat.lon);
    r_sta_ecf = const.R_E * [cl*cd; cl*sd; sl];
    v_sta_ecf = cross([0;0;const.we], r_sta_ecf);
    
    % Chain Rotation to EMO2000
    % r_emo = R_EME_EMO * R_ECF_ECI * r_ecf
    R_ECF_EMO = const.R_EME_EMO * R_ECF_ECI;
    
    r_sta_emo_earth = R_ECF_EMO * r_sta_ecf;
    v_sta_emo_earth = R_ECF_EMO * v_sta_ecf;
    
    % Shift Origin to Sun
    r_sta = r_E_Sun_EMO + r_sta_emo_earth;
    v_sta = v_E_Sun_EMO + v_sta_emo_earth;
    
    % 3. Compute Measurement
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
    
    % 4. H Matrix
    H_r = u_rho';
    H_v = zeros(1,3);
    H_rr_r = (rel_vel' - range_rate * u_rho') / rho;
    H_rr_v = u_rho';
    
    H = [H_r, H_v, 0, 0; 
         H_rr_r, H_rr_v, 0, H_b];
end