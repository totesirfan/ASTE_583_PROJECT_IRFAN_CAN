function dX = equations_of_motion(t, X, const)
    % EQUATIONS_OF_MOTION Sun-Centered EMO2000 Frame
    
    r_sc = X(1:3);
    v_sc = X(4:6);
    k_SRP = X(7);
    
    % --- 1. Get Earth State (Standardize to EMO2000) ---
    % SPICE gives J2000 (EME). We rotate it.
    try
        st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
        r_E_EMO = const.R_EME_EMO * st_earth_J2000(1:3); % Rotate Position
        % Velocity not needed for force model, but good practice
    catch
        error('MICE toolkit not loaded.');
    end
    
    % --- 2. Forces in EMO2000 ---
    r_mag = norm(r_sc);
    
    % Sun Gravity
    a_sun = -const.mu_S * r_sc / r_mag^3;
    
    % Earth Gravity (Third Body)
    % Vectors must be in same frame (EMO2000)
    r_rel = r_sc - r_E_EMO; % SC relative to Earth
    a_earth = -const.mu_E * ( r_rel / norm(r_rel)^3 + r_E_EMO / norm(r_E_EMO)^3 );
    
    % SRP
    AU_km = 149597870.7;
    P_dist = const.P_SRP * (AU_km / r_mag)^2;
    a_srp_mag = (k_SRP * P_dist * (1 + const.rho_r) * const.A_sc / const.m_sc) / 1000;
    a_srp = a_srp_mag * (r_sc / r_mag);
    
    dX = [v_sc; a_sun + a_earth + a_srp; 0; 0];
end
