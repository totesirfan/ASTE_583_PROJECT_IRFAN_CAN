function Phase2_Final_Kalman(prior_source)
tic
% Phase2_Final_Kalman.m
% FINAL 10-STATE KALMAN (CKF) OD USING 0–14 DAY DATA
% Optimization: Vectorized interpolation to remove loop overhead.
%
% prior_source (optional): 'kalman' | 'batch' | 'handout'
%   'kalman'  -> use ASTE583_PrelimKalman_Results.mat if available (default)
%   'batch'   -> use ASTE583_PrelimBatch_Results.mat if available
%   'handout' -> use lib_constants() handout a priori only

    if nargin < 1
        prior_source = 'handout';   % default behavior
    end
    prior_source = lower(strtrim(prior_source));

    clearvars -except prior_source; clc; close all;
    fprintf('=== PHASE 2: FINAL 10-STATE KALMAN OD (0–14 DAYS) ===\n');
    
    %% 1. Setup & Constants
    try
        init_project();
    catch
    end
    const   = lib_constants();
    t0_et   = const.t_detect_et;
    n_state = 10;
    
    % Default a priori from handout / lib_constants
    X0_ref = [const.X0_ref; ...
              const.k_SRP_0; ...
              0.0; ...
              const.stations(4).lat; ...
              const.stations(4).lon];
    P0_ref = const.P0_aug;
    used_prior = 'Handout a priori (lib_constants)';
    
    %% 1.1 Select a priori source based on user choice
    switch prior_source
        case 'kalman'
            if exist('ASTE583_PrelimKalman_Results.mat','file')
                S = load('ASTE583_PrelimKalman_Results.mat');
                % Support multiple possible field names
                if isfield(S,'X0_final_kalman')
                    X0_ref = S.X0_final_kalman(:);
                elseif isfield(S,'X0_kalman')
                    X0_ref = S.X0_kalman(:);
                elseif isfield(S,'kf_results') && isfield(S.kf_results,'X0_KF')
                    X0_ref = S.kf_results.X0_KF(:);
                end
                
                if isfield(S,'P0_final_kalman')
                    P0_ref = S.P0_final_kalman;
                elseif isfield(S,'P0_kalman')
                    P0_ref = S.P0_kalman;
                elseif isfield(S,'kf_results') && isfield(S.kf_results,'P0_KF')
                    P0_ref = S.kf_results.P0_KF;
                end
                
                used_prior = 'Preliminary Kalman solution';
            else
                warning('Prelim Kalman file not found. Falling back to handout a priori.');
            end
            
        case 'batch'
            if exist('ASTE583_PrelimBatch_Results.mat','file')
                S = load('ASTE583_PrelimBatch_Results.mat');
                % Support either flat vars or a struct
                if isfield(S,'prelim_results')
                    pr = S.prelim_results;
                    if isfield(pr,'X0_batch')
                        X0_ref = pr.X0_batch(:);
                    end
                    if isfield(pr,'P0_batch')
                        P0_ref = pr.P0_batch;
                    end
                else
                    if isfield(S,'X0_batch')
                        X0_ref = S.X0_batch(:);
                    end
                    if isfield(S,'P0_batch')
                        P0_ref = S.P0_batch;
                    end
                end
                used_prior = 'Preliminary Batch solution';
            else
                warning('Prelim Batch file not found. Falling back to handout a priori.');
            end
            
        case 'handout'
            % Already set above; nothing to change
            used_prior = 'Handout a priori (lib_constants)';
            
        otherwise
            warning('Unknown prior_source "%s". Using handout a priori.', prior_source);
    end
    
    if numel(X0_ref) ~= n_state
        error('A priori state must be 10x1; got %d elements.', numel(X0_ref));
    end
    fprintf('Using %s as a priori for FINAL Kalman.\n', used_prior);
    
    %% 2. Load ALL measurements, then restrict to 0–14 days
    fprintf('Loading measurements (full set)...\n');
    meas = load_project_measurements();   
    t_sec   = meas.time_sec(:);          
    st_id   = meas.station_id(:);
    rng_km  = meas.range_km(:);
    rr_kmps = meas.rr_kmps(:);
    has_range_all = ~isnan(rng_km);
    
    t_days_all = t_sec / const.day2sec;
    t_et_all   = t0_et + t_sec;
    
    t_max_days = 14.0;
    use_idx = (t_days_all >= 0.0) & (t_days_all <= t_max_days);
    
    t_days_all = t_days_all(use_idx);
    t_et_all   = t_et_all(use_idx);
    st_id      = st_id(use_idx);
    rng_km     = rng_km(use_idx);
    rr_kmps    = rr_kmps(use_idx);
    has_range  = has_range_all(use_idx);
    
    % Sort by time
    [t_et_all, sort_idx] = sort(t_et_all);
    t_days_all = t_days_all(sort_idx);
    st_id      = st_id(sort_idx);
    rng_km     = rng_km(sort_idx);
    rr_kmps    = rr_kmps(sort_idx);
    has_range  = has_range(sort_idx);
    
    n_meas      = numel(t_et_all);
    t_last_days = max(t_days_all);
    
    fprintf('Final Kalman OD: using %d measurements up to %.3f days.\n', ...
            n_meas, t_last_days);
            
    %% 3. Measurement noise (km, km/s)
    sigma_r_DSN   = 1e-3;    % km
    sigma_rr_DSN  = 1e-7;    % km/s
    sigma_r_ANT   = 1e-2;    % km
    sigma_rr_ANT  = 1e-6;    % km/s
    
    %% 4. CKF / EKF Settings
    options_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-9);
    max_iter  = 9;        % relinearization iterations
    tol_dx0   = 1e-3;      % convergence threshold on ||ΔX0||
    
    % Storage for final (converged) pre/postfits
    rho_prefit  = NaN(n_meas,1);
    rho_postfit = NaN(n_meas,1);
    rr_prefit   = NaN(n_meas,1);
    rr_postfit  = NaN(n_meas,1);
    
    %% 5. Relinearizing CKF at detection epoch
    X0_curr = X0_ref(:);
    P0_curr = P0_ref;
    
    for iter = 1:max_iter
        fprintf('\n=== FINAL 10-STATE KALMAN ITERATION %d ===\n', iter);
        
        % 5.1 Propagate nominal state + STM from t0 to last measurement
        t_last = t_et_all(end);
        Phi0 = eye(n_state);
        Z0   = [X0_curr; Phi0(:)];
        
        fprintf('  [Iter %d] Propagating reference & STM from t0 to %.3f days...\n', ...
                iter, (t_last - t0_et)/const.day2sec);
                
        [T_prop, Z_prop] = ode45(@(t,Z) lib_dynamics(t, Z, const), ...
                                 [t0_et, t_last], Z0, options_ode);
                                 
        X_hist   = Z_prop(:, 1:n_state);        % nominal state history
        Phi_hist = Z_prop(:, n_state+1:end);    % STM history (flattened)
        
        fprintf('  [Iter %d] Propagation complete. Steps: %d. Pre-interpolating states...\n', ...
                iter, numel(T_prop));
        
        % Vectorized interpolation at all measurement times
        X_meas_all   = interp1(T_prop, X_hist,   t_et_all, 'spline');
        Phi_meas_all = interp1(T_prop, Phi_hist, t_et_all, 'spline');

        % 5.2 Sequential CKF at t0
        fprintf('  [Iter %d] Starting CKF measurement loop over %d measurements...\n', ...
                iter, n_meas);
                
        deltaX0 = zeros(n_state,1);   % deviation at t0
        P0      = P0_curr;            % covariance at t0
        I10 = eye(n_state);
        
        % Per-iteration storage
        rho_prefit_iter  = NaN(n_meas,1);
        rho_postfit_iter = NaN(n_meas,1);
        rr_prefit_iter   = zeros(n_meas,1);
        rr_postfit_iter  = zeros(n_meas,1);
        
        debug_stride = max(1, floor(n_meas / 2));
        
        for k = 1:n_meas
            t_k   = t_et_all(k);
            stat  = st_id(k);
            use_range = has_range(k);
            rho_o = rng_km(k);
            rr_o  = rr_kmps(k);
            
            % Use pre-interpolated state & STM
            X_k_flat   = X_meas_all(k, :)';
            Phi_k_flat = Phi_meas_all(k, :)';
            Phi_k0     = reshape(Phi_k_flat, n_state, n_state);  
            
            % Measurement model at nominal state
            [Y_comp_nom, H_tilde] = lib_measurements(t_k, X_k_flat, stat, const);
            H_full = H_tilde * Phi_k0;   % 2x10
            
            if use_range
                % 2x1 measurement (rho, rr)
                y_obs = [rho_o; rr_o];
                dy    = y_obs - Y_comp_nom;   % prefit residual
                Hk    = H_full;              % 2x10
                
                if stat == 4
                    Rk = diag([sigma_r_ANT^2; sigma_rr_ANT^2]);
                else
                    Rk = diag([sigma_r_DSN^2; sigma_rr_DSN^2]);
                end
                
                rho_prefit_iter(k) = dy(1);
                rr_prefit_iter(k)  = dy(2);
                
                dy_pred = Hk * deltaX0;
                innov   = dy - dy_pred;
                
                S_k = Hk * P0 * Hk' + Rk;
                K0  = (P0 * Hk') / S_k;   % 10x2
                
                deltaX0 = deltaX0 + K0 * innov;
                P0      = (I10 - K0 * Hk) * P0;
                
                dy_post = dy - Hk * deltaX0;
                rho_postfit_iter(k) = dy_post(1);
                rr_postfit_iter(k)  = dy_post(2);
                
            else
                % Range-rate only (1x1)
                dy = rr_o - Y_comp_nom(2);
                Hk = H_full(2,:);      % 1x10
                
                if stat == 4
                    Rk = sigma_rr_ANT^2;
                else
                    Rk = sigma_rr_DSN^2;
                end
                
                rr_prefit_iter(k) = dy;
                
                dy_pred = Hk * deltaX0;
                innov   = dy - dy_pred;
                
                S_k = Hk * (P0 * Hk') + Rk;   % scalar
                K0  = (P0 * Hk') / S_k;       % 10x1
                
                deltaX0 = deltaX0 + K0 * innov;
                P0      = (I10 - K0 * Hk) * P0;
                
                dy_post = dy - Hk * deltaX0;
                rr_postfit_iter(k) = dy_post;
            end
            
            if mod(k, debug_stride) == 0 || k == n_meas
                fprintf('    [Iter %d] Processed %5d / %5d measurements\n', ...
                        iter, k, n_meas);
            end
        end
        
        % 5.3 Update reference at t0 and check convergence
        dX0_norm = norm(deltaX0);
        
        rr_prefit_rms_mm  = sqrt(mean(rr_prefit_iter.^2)) * 1e6;
        rr_postfit_rms_mm = sqrt(mean(rr_postfit_iter.^2)) * 1e6;
        
        rho_idx            = ~isnan(rho_prefit_iter);
        rho_prefit_rms_km  = sqrt(mean(rho_prefit_iter(rho_idx).^2));
        rho_postfit_rms_km = sqrt(mean(rho_postfit_iter(rho_idx).^2));
        
        fprintf('  [Iter %d] ||ΔX0||_2 = %.3e\n', iter, dX0_norm);
        fprintf('  [Iter %d] RR prefit  RMS ≈ %.2f mm/s, postfit RMS ≈ %.2f mm/s\n', ...
                iter, rr_prefit_rms_mm, rr_postfit_rms_mm);
        fprintf('  [Iter %d] R  prefit  RMS ≈ %.3f km,   postfit RMS ≈ %.3f km\n', ...
                iter, rho_prefit_rms_km, rho_postfit_rms_km);
        
        X0_new = X0_curr + deltaX0;
        P0_new = P0;
        
        rho_prefit  = rho_prefit_iter;
        rho_postfit = rho_postfit_iter;
        rr_prefit   = rr_prefit_iter;
        rr_postfit  = rr_postfit_iter;
        
        X0_curr = X0_new;
        P0_curr = P0_new;
        
        if dX0_norm < tol_dx0
            fprintf('  [Iter %d] Converged: ||ΔX0|| < %.1e. Stopping relinearization.\n', ...
                    iter, tol_dx0);
            break;
        end
    end
    
    %% 6. Final solution at detection epoch
    X0_final = X0_curr;
    P0_final = P0_curr;
    
    r0 = X0_final(1:3);
    v0 = X0_final(4:6);
    k_SRP_final = X0_final(7);
    bias_final  = X0_final(8);
    lat4_final  = X0_final(9);
    lon4_final  = X0_final(10);
    
    fprintf('\n=== FINAL 10-STATE KALMAN RESULT (t0: DETECTION EPOCH) ===\n');
    fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', r0(1), r0(2), r0(3));
    fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', v0(1), v0(2), v0(3));
    fprintf('k_SRP:        %.6f\n', k_SRP_final);
    fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', ...
            bias_final, bias_final*1e6);
    fprintf('lat_4 (deg):  %.6f\n', lat4_final * const.rad2deg);
    fprintf('lon_4 (deg):  %.6f\n', lon4_final * const.rad2deg);
    
    %% 7. Save results for Final OD plots & LTM mapping
    X0_final_kalman = X0_final;
    P0_final_kalman = P0_final;
    
    save('ASTE583_FinalKalman_Results.mat', ...
         'X0_final_kalman', 'P0_final_kalman', ...
         't_et_all', 'st_id', 'has_range', ...
         'rho_prefit', 'rho_postfit', ...
         'rr_prefit', 'rr_postfit');
         
    fprintf('\nSaved FINAL Kalman results to ASTE583_FinalKalman_Results.mat\n');
    toc
end
