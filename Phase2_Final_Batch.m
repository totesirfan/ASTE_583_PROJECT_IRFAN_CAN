function Phase2_Final_Batch(prior_source)
% Phase2_Final_Batch.m
% Final 10-state Batch OD using full 0–14 day radiometric arc.
%
%   X0 = [r0(3); v0(3); k_SRP; rho_dot_bias; lat_4; lon_4]  @ detection epoch
%
% Usage:
%   Phase2_Final_Batch;                 % uses 'handout' prior
%   Phase2_Final_Batch('handout');      % same as above
%   Phase2_Final_Batch('batch');        % prior from Prelim Batch result
%   Phase2_Final_Batch('kalman');       % prior from Prelim Kalman result
%
% Output:
%   ASTE583_FinalBatch_Results.mat with:
%       final_batch.X0_prior        (10×1)
%       final_batch.P0_prior        (10×10)
%       final_batch.X0_batch        (10×1 final estimate)
%       final_batch.P0_batch        (10×10 posterior covariance)
%       final_batch.time_days       (n×1 time since detection in days)
%       final_batch.st_id           (n×1 station IDs)
%       final_batch.prefit_rr       (n×1 RR prefit, km/s, final iter)
%       final_batch.postfit_rr      (n×1 RR postfit, km/s, final iter)
%       final_batch.iteration_results (struct array for all iterations)
%
% Dependencies:
%   init_project.m, lib_constants.m, lib_dynamics.m,
%   lib_measurements.m, load_project_measurements.m

    if nargin < 1
        prior_source = 'handout';
    end

    clearvars -except prior_source;
    clc; close all;
    tic;

    fprintf('=== Phase 2: FINAL 10-state Batch OD (0–14 days) ===\n');
    fprintf('Prior source: %s\n', prior_source);

    %% 1. Setup
    init_project();
    C      = lib_constants();
    t0_et  = C.t_detect_et;

    DCO_days = 14.0;                 % full arc
    t_DCO    = DCO_days * C.day2sec;

    %% 2. Measurements (0–14 days)
    meas_full = load_project_measurements();  % 0–14 day dataset

    use      = meas_full.time_sec <= t_DCO;
    t_sec    = meas_full.time_sec(use);
    st_id    = meas_full.station_id(use);
    rng_km   = meas_full.range_km(use);
    rr_kmps  = meas_full.rr_kmps(use);

    % Sort by time (just in case)
    [t_sec, ksort] = sort(t_sec);
    st_id   = st_id(ksort);
    rng_km  = rng_km(ksort);
    rr_kmps = rr_kmps(ksort);

    n_meas   = numel(t_sec);
    time_days= t_sec / C.day2sec;
    t_et_all = t0_et + t_sec;

    fprintf('Final Batch OD: %d meas, 0–%.3f days.\n', n_meas, DCO_days);

    %% 3. A priori state & covariance (10-state) from selected prior
    st4 = C.stations(4);

    switch lower(strtrim(prior_source))
        case 'handout'
            % Handout prior (same as prelim batch start)
            X0_prior = [ C.X0_ref; ...      % 1–6: r0, v0
                         C.k_SRP_0; ...     % 7:   SRP scale
                         0.0;       ...     % 8:   Doppler bias (km/s)
                         st4.lat;   ...     % 9:   lat_4 (rad)
                         st4.lon ];         % 10:  lon_4 (rad)
            P0_prior = C.P0_aug;

        case 'batch'
            % Use Prelim Batch result as prior
            S = load('ASTE583_PrelimBatch_Results.mat');
            pr = S.prelim_results;
            X0_prior = pr.X0_batch(:);
            P0_prior = pr.P0_batch;

        case 'kalman'
            % Use Prelim Kalman result as prior
            if exist('ASTE583_PrelimKalman_Results.mat','file')
                S = load('ASTE583_PrelimKalman_Results.mat');
            else
                S = load('ASTE583_PrelimKalman_Debug.mat');
            end
            kr = S.kf_results;
            X0_prior = kr.X0_KF(:);
            P0_prior = kr.P0_KF;

        otherwise
            error('Unknown prior_source "%s". Use ''handout'', ''batch'', or ''kalman''.', prior_source);
    end

    X0_est = X0_prior;

    % P0_inv via Cholesky
    [L0,p0] = chol(P0_prior,'lower');
    if p0 ~= 0
        error('P0_prior is not positive definite.');
    end
    L0i    = L0 \ eye(size(L0));
    P0_inv = L0i' * L0i;

    %% 4. Measurement noise
    sigma_r_DSN  = 1e-3;   % km   (1 m)
    sigma_rr_DSN = 1e-7;   % km/s (0.1 mm/s)
    sigma_r_ANT  = 1e-2;   % km   (10 m)
    sigma_rr_ANT = 1e-6;   % km/s (1 mm/s)

    Rinv_DSN = diag([1/sigma_r_DSN^2,  1/sigma_rr_DSN^2]);
    Rinv_ANT = diag([1/sigma_r_ANT^2,  1/sigma_rr_ANT^2]);

    %% 5. ODE options
    ode_opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

    %% 6. Batch iterations over full arc
    max_iter   = 9;
    tol_update = 1e-6;
    n_state    = 10;

    iteration_results = struct([]);
    fprintf('\nBeginning Final Batch iterations...\n');

    for iter = 1:max_iter
        fprintf('\n=== FINAL BATCH ITER %d ===\n', iter);

        X0_full = X0_est;
        Phi0    = eye(n_state);

        % 6.1 Single propagation for all measurement times (state + STM)
        % State dim: 10; STM dim: 10×10 = 100
        X_aug0 = [X0_full; Phi0(:)];      % 10 + 100 = 110×1

        % Handle measurements at t=0 separately
        is_t0    = abs(t_sec) < 1e-9;
        t_pos    = t_sec(~is_t0);
        t_et_pos = t0_et + t_pos;
        t_et_u   = unique(t_et_pos);

        X_prop = [];
        if ~isempty(t_et_u)
            tspan = [t0_et; t_et_u(:)];
            [~, X_prop] = ode45(@(t,X) lib_dynamics(t,X,C), ...
                                 tspan, X_aug0, ode_opts);
        end

        % Allocate containers
        Xk_all  = zeros(n_state, n_meas);
        Phi_all = zeros(n_state, n_state, n_meas);

        % Measurements at t = 0
        if any(is_t0)
            idx0 = find(is_t0);
            for j = 1:numel(idx0)
                k0 = idx0(j);
                Xk_all(:,k0)    = X0_full;
                Phi_all(:,:,k0) = Phi0;
            end
        end

        % Measurements at t > 0
        if ~isempty(t_et_u)
            [~, idx_u_for_meas] = ismember(t_et_pos, t_et_u);
            idx_meas_pos        = find(~is_t0);

            for j = 1:numel(idx_meas_pos)
                k_meas = idx_meas_pos(j);
                iu     = idx_u_for_meas(j);

                row     = iu + 1;   % +1 for initial t0 row in tspan
                X_aug_k = X_prop(row,:).';

                Xk_all(:,k_meas)    = X_aug_k(1:n_state);
                Phi_all(:,:,k_meas) = reshape(X_aug_k(n_state+1:end), n_state, n_state);
            end
        end

        % 6.2 Normal equations with prior
        Delta   = P0_inv;
        N_vec   = P0_inv * (X0_prior - X0_full);
        prefit_rr  = zeros(n_meas,1);

        for k = 1:n_meas
            t_k = t_et_all(k);
            st  = st_id(k);

            X_k   = Xk_all(:,k);
            Phi_k = Phi_all(:,:,k);

            [Y_comp, H_tilde] = lib_measurements(t_k, X_k, st, C);   % [rho; rho_dot]
            H_full = H_tilde * Phi_k;                                % 2×10

            rho_obs = rng_km(k);
            rr_obs  = rr_kmps(k);
            has_r   = ~isnan(rho_obs);

            if has_r
                dy = [rho_obs - Y_comp(1);
                      rr_obs  - Y_comp(2)];
                Hk = H_full;

                if st == 4
                    Rinv = Rinv_ANT;
                else
                    Rinv = Rinv_DSN;
                end

                prefit_rr(k) = dy(2);

                Delta = Delta + Hk' * (Rinv * Hk);
                N_vec = N_vec + Hk' * (Rinv * dy);

            else
                dy = rr_obs - Y_comp(2);
                Hk = H_full(2,:);

                if st == 4
                    Rinv_rr = 1 / sigma_rr_ANT^2;
                else
                    Rinv_rr = 1 / sigma_rr_DSN^2;
                end

                prefit_rr(k) = dy;

                Delta = Delta + (Hk' * Hk) * Rinv_rr;
                N_vec = N_vec + (Hk' * dy) * Rinv_rr;
            end
        end

        % 6.3 Solve for ΔX0 (normal equations)
        delta_X0 = Delta \ N_vec;
        d_norm   = norm(delta_X0);
        fprintf('  ||ΔX0||_2 = %.3e\n', d_norm);

        X0_est = X0_est + delta_X0;

        % Posterior covariance via Cholesky of Delta
        [L,p] = chol(Delta,'lower');
        if p ~= 0
            warning('Delta not PD; using inv(Delta) for P_post.');
            P_post = inv(Delta); %#ok<MINV>
        else
            Li    = L \ eye(size(L));
            P_post = Li' * Li;
        end

        % 6.4 Postfit Doppler (linearized about updated X0_est)
        postfit_rr = zeros(n_meas,1);
        for k = 1:n_meas
            t_k = t_et_all(k);
            st  = st_id(k);

            X_k   = Xk_all(:,k);
            Phi_k = Phi_all(:,:,k);

            X_k_upd = X_k + Phi_k * delta_X0;
            [Y_upd, ~] = lib_measurements(t_k, X_k_upd, st, C);

            postfit_rr(k) = rr_kmps(k) - Y_upd(2);
        end

        % 6.5 Diagnostics
        rr_pre_rms  = sqrt(mean(prefit_rr.^2));
        rr_post_rms = sqrt(mean(postfit_rr.^2));
        fprintf('  RR RMS pre/post: %.3e / %.3e km/s (%.2f / %.2f mm/s)\n', ...
            rr_pre_rms, rr_post_rms, rr_pre_rms*1e6, rr_post_rms*1e6);

        % Store iteration data
        iteration_results(iter).X0_est     = X0_est;
        iteration_results(iter).P_post     = P_post;
        iteration_results(iter).time_days  = time_days;
        iteration_results(iter).prefit_rr  = prefit_rr;
        iteration_results(iter).postfit_rr = postfit_rr;

        % 6.6 Convergence check
        if d_norm < tol_update
            fprintf('  Converged after %d iterations.\n', iter);
            break;
        end
    end

    % Final estimate from last iteration
    X0_est_final  = X0_est;
    P0_post_final = iteration_results(end).P_post;
    prefit_final  = iteration_results(end).prefit_rr;
    postfit_final = iteration_results(end).postfit_rr;

    %% 7. Summary
    fprintf('\n=== FINAL 10-STATE BATCH RESULT (t0) ===\n');
    fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', X0_est_final(1:3));
    fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', X0_est_final(4:6));
    fprintf('k_SRP:        %.6f\n',  X0_est_final(7));
    fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', ...
            X0_est_final(8), X0_est_final(8)*1e6);
    fprintf('lat_4 (deg):  %.6f\n', X0_est_final(9)*C.rad2deg);
    fprintf('lon_4 (deg):  %.6f\n', X0_est_final(10)*C.rad2deg);

    %% 8. Save
    final_batch = struct();
    final_batch.X0_prior         = X0_prior;
    final_batch.P0_prior         = P0_prior;
    final_batch.X0_batch         = X0_est_final;
    final_batch.P0_batch         = P0_post_final;
    final_batch.time_days        = time_days;
    final_batch.st_id            = st_id;
    final_batch.prefit_rr        = prefit_final;
    final_batch.postfit_rr       = postfit_final;
    final_batch.iteration_results = iteration_results;

    save('ASTE583_FinalBatch_Results.mat', 'final_batch', 'meas_full');

    fprintf('\nSaved final batch to ASTE583_FinalBatch_Results.mat\n');
    toc;
end
