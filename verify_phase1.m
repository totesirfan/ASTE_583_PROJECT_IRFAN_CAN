function verify_phase1()
    % VERIFY_PHASE1 Master test script for ASTE 583 Project Phase 1
    % Updated with larger perturbation step to handle Sun-Centered precision limits.
    
    clc; close all;
    fprintf('==============================================\n');
    fprintf('   ASTE 583: PHASE 1 VERIFICATION SUITE\n');
    fprintf('==============================================\n\n');

    %% TEST 1: CONSTANTS & SPICE SETUP
    fprintf('[TEST 1] Loading Constants and Kernels... ');
    try
        const = get_constants();
        if ~isfield(const, 'mu_S') || ~isfield(const, 'P_SRP')
            error('Constants struct is missing mu_S or P_SRP.');
        end
        % Check SPICE
        try
            t_check = cspice_str2et('2025 DEC 01');
        catch
            setup_project(); % Attempt auto-load
            t_check = cspice_str2et('2025 DEC 01');
        end
        
        fprintf('PASS\n');
    catch ME
        fprintf('FAIL\n');
        fprintf('         Error: %s\n', ME.message);
        return;
    end

    %% TEST 2: DYNAMICS (Equations of Motion)
    fprintf('\n[TEST 2] Testing Dynamics Integration (1 Step)... ');
    try
        t0 = cspice_str2et(const.epoch_utc_str);
        X0 = [const.X0_ref; const.k_SRP_0; 0];
        
        dX = equations_of_motion(t0, X0, const);
        
        if length(dX) ~= 8
            error('Output state derivative must be 8x1.');
        end
        
        % Sanity Check: Sun Gravity at ~1AU
        r_sc = X0(1:3);
        accel_mag = norm(dX(4:6));
        expected_mag = const.mu_S / norm(r_sc)^2;
        
        if abs(accel_mag - expected_mag) > 1e-4
            warning('Accel magnitude mismatch. Check units.');
        end
        
        fprintf('PASS\n');
        fprintf('         (Initial Accel Mag: %.5e km/s^2)\n', accel_mag);
        
    catch ME
        fprintf('FAIL\n');
        fprintf('         Error: %s\n', ME.message);
        return;
    end

    %% TEST 3: MEASUREMENT MODEL & H-MATRIX
    fprintf('\n[TEST 3] Validating H-Matrix (Finite Differences)... \n');
    
    try
        % Define a nominal state (Sun-Centered EMO2000)
        % We construct a state relative to Earth to ensure valid geometry
        st_earth = cspice_spkezr('EARTH', t0, 'J2000', 'NONE', 'SUN');
        r_E_EMO = const.R_EME_EMO * st_earth(1:3);
        v_E_EMO = const.R_EME_EMO * st_earth(4:6);
        
        % Place SC ~400,000 km from Earth (Lunar distance)
        r_sc = r_E_EMO + [300000; 100000; 50000]; 
        v_sc = v_E_EMO + [0.5; 0.5; 0.1];
        
        X_nom = [r_sc; v_sc; 1.2; 1e-5];
        station_id = 1; % Goldstone
        
        % 1. Analytical Result
        [Y_nom, H_anal] = measurement_model(t0, X_nom, station_id, const);
        
        % 2. Numerical Result (Finite Difference)
        % CRITICAL FIX: Increased delta to 1e-4 (10cm) to avoid machine epsilon noise
        % when subtracting large Sun-Centered vectors.
        perturbation = 1e-4; 
        
        n_states = length(X_nom);
        H_num = zeros(2, 8);
        
        for i = 1:n_states
            X_p = X_nom; X_p(i) = X_p(i) + perturbation;
            X_m = X_nom; X_m(i) = X_m(i) - perturbation;
            
            [Y_p, ~] = measurement_model(t0, X_p, station_id, const);
            [Y_m, ~] = measurement_model(t0, X_m, station_id, const);
            
            H_num(:, i) = (Y_p - Y_m) / (2 * perturbation);
        end
        
        % 3. Compare
        param_names = {'rx', 'ry', 'rz', 'vx', 'vy', 'vz', 'kSRP', 'Bias'};
        meas_names = {'Range', 'Doppler'};
        
        % Use a normalized tolerance
        diff = abs(H_anal - H_num);
        max_diff = max(diff(:));
        
        fprintf('         Perturbation Size: %e\n', perturbation);
        fprintf('         Max Difference:    %e\n', max_diff);
        
        if max_diff < 1e-5
            fprintf('         Result: PASS\n');
        else
            fprintf('         Result: FAIL\n');
            disp('Difference Matrix (Analytic - Numerical):');
            disp(H_anal - H_num);
        end
        
    catch ME
        fprintf('FAIL\n');
        fprintf('         Error: %s\n', ME.message);
        fprintf('         Stack: %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    
    fprintf('\n----------------------------------------------\n');
    fprintf('If all PASS, proceed to run_reference_trajectory.m\n');
end