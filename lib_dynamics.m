function dX = lib_dynamics(t, X, const)
    % LIB_DYNAMICS Computes state derivative in Sun-Centered EMO2000.
    
    r_sc = X(1:3);
    v_sc = X(4:6);
    k_SRP = X(7);
    
    % 1. Get Earth State (Rotate J2000 -> EMO2000)
    % Uses cspice_spkezr for high precision
    st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_EMO = const.R_EME_EMO * st_earth_J2000(1:3);
    
    % 2. Forces
    r_mag = norm(r_sc);
    
    % Sun Gravity
    a_sun = -const.mu_S * r_sc / r_mag^3;
    
    % Earth Gravity (Third Body)
    r_rel = r_sc - r_E_EMO; 
    a_earth = -const.mu_E * ( r_rel / norm(r_rel)^3 + r_E_EMO / norm(r_E_EMO)^3 );
    
    % SRP
    P_dist = const.P_SRP * (const.AU / r_mag)^2;
    a_srp_mag = (k_SRP * P_dist * (1 + const.rho_r) * const.A_sc / const.m_sc) / 1000;
    a_srp = a_srp_mag * (r_sc / r_mag);
    
    dX = [v_sc; a_sun + a_earth + a_srp; 0; 0];
end