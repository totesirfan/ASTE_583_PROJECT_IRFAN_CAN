function dX_out = lib_dynamics(t, X, const)
    % LIB_DYNAMICS Physics engine for Lunar Trailblazer.
    % Handles both:
    %   1. State Propagation (Phase 1): X is 7x1 to 10x1
    %   2. State + STM Propagation (Phase 2): X is 110x1
    %
    % Input:
    %   t : Ephemeris Time (seconds past J2000)
    %   X : State Vector (size N)
    %       - if N=7-10: [r; v; params...]
    %       - if N=110 : [r; v; params...; STM(100)]
    %
    % Output:
    %   dX_out : Derivative vector (size N)

    %% 1. Parse Input State
    n_total = length(X);
    
    % Determine mode based on vector length
    % Phase 2 state is 10 (state) + 100 (STM) = 110 elements
    calc_stm = (n_total > 20); 
    
    % Extract core state (Always the first 7-10 elements)
    if calc_stm
        n_state = 10; % Fixed size for Phase 2 augmented state
        state = X(1:n_state);
        Phi = reshape(X(n_state+1:end), n_state, n_state);
    else
        n_state = n_total;
        state = X;
    end
    
    r_sc  = state(1:3);
    v_sc  = state(4:6);
    k_SRP = state(7);
    % bias, lat, lon are state(8:10) but don't affect dynamics directly
    
    %% 2. Ephemeris (Sun-Centered EMO2000)
    % Get Earth position
    st_earth = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_EMO = const.R_EME_EMO * st_earth(1:3);
    
    %% 3. Compute Forces (Physics)
    % Magnitudes and Relative Vectors
    r_mag = norm(r_sc);
    r_rel = r_sc - r_E_EMO; % SC relative to Earth
    r_rel_mag = norm(r_rel);
    r_E_mag = norm(r_E_EMO);
    
    % A. Sun Gravity (Point Mass)
    a_sun = -const.mu_S * r_sc / r_mag^3;
    
    % B. Earth Gravity (Third Body Perturbation)
    % a_3rd = -mu_3 * ( r_rel/|r_rel|^3 + r_3/|r_3|^3 )
    a_earth = -const.mu_E * (r_rel / r_rel_mag^3 + r_E_EMO / r_E_mag^3);
    
    % C. Solar Radiation Pressure (Cannonball)
    % P_SRP scaled by inverse square of distance from Sun
    % C_srp factor includes Area, Mass, Reflectivity
    C_srp_base = const.P_SRP * const.AU^2 * (1 + const.rho_r) * const.A_sc / const.m_sc / 1000; % km/s^2 units
    a_srp = (k_SRP * C_srp_base) * r_sc / r_mag^3;
    
    % Total Acceleration
    total_acc = a_sun + a_earth + a_srp;
    
    %% 4. State Derivatives
    dX_state = zeros(n_state, 1);
    dX_state(1:3) = v_sc;       % dr/dt = v
    dX_state(4:6) = total_acc;  % dv/dt = a
    % Parameters (k_SRP, bias, etc.) are constant -> derivatives are 0
    
    %% 5. STM Derivatives (Only if requested)
    if calc_stm
        % Initialize Jacobian A (10x10)
        % A = [ 0   I   0 ]
        %     [ G   0   S ]
        %     [ 0   0   0 ]
        A = zeros(10, 10);
        A(1:3, 4:6) = eye(3); 
        
        I3 = eye(3);
        
        % --- Gravity Gradients (G = da/dr) ---
        % Sun: mu/r^5 * (3*r*r' - r^2*I)
        G_sun = (const.mu_S / r_mag^5) * (3 * (r_sc * r_sc') - r_mag^2 * I3);
        
        % Earth: mu/r_rel^5 * (3*r_rel*r_rel' - r_rel^2*I)
        G_earth = (const.mu_E / r_rel_mag^5) * (3 * (r_rel * r_rel') - r_rel_mag^2 * I3);
        
        % SRP Gradient
        % deriv of (C * r / r^3) -> C/r^5 * (r^2*I - 3*r*r')
        % Note: Signs are opposite of gravity gradient
        G_srp = (k_SRP * C_srp_base / r_mag^5) * (r_mag^2 * I3 - 3 * (r_sc * r_sc'));
        
        % Total Position Partial
        A(4:6, 1:3) = G_sun + G_earth + G_srp;
        
        % --- Parameter Sensitivities (S = da/dp) ---
        % Sensitivity to k_SRP: a_srp / k_SRP
        A(4:6, 7) = a_srp / k_SRP;
        
        % Compute STM Derivative: dPhi = A * Phi
        dPhi = A * Phi;
        
        % Output: Stacked [State_Deriv; STM_Deriv_Flattened]
        dX_out = [dX_state; dPhi(:)];
        
    else
        % Phase 1 Mode: Return only state derivatives
        dX_out = dX_state;
    end
end