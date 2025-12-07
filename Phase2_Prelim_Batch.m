% Phase2_Prelim_Batch.m
% Preliminary 10-state Batch LS OD using radiometric data 0–6 days.
% X0 = [r0(3); v0(3); k_SRP; rho_dot_bias; lat_4; lon_4]

clear; clc; close all;
tic
%% 1. Setup
init_project();
C     = lib_constants();
t0_et = C.t_detect_et;

DCO_days = 6.0;
t_DCO    = DCO_days * C.day2sec;

%% 2. Measurements (0–6 days only for batch, full set saved later)
meas_full = load_project_measurements();    % 0–14 days

use      = meas_full.time_sec <= t_DCO;
t_sec    = meas_full.time_sec(use);
st_id    = meas_full.station_id(use);
rng_km   = meas_full.range_km(use);
rr_kmps  = meas_full.rr_kmps(use);
n_meas   = numel(t_sec);

fprintf('Prelim Batch OD: %d meas, 0–%.3f days.\n', n_meas, DCO_days);

time_days = t_sec / C.day2sec;
t_et_all  = t0_et + t_sec;

% For later plotting / station-wise residuals
has_range = ~isnan(rng_km);

%% 3. A priori state & covariance (10-state)
st4 = C.stations(4);

X0_est = [ C.X0_ref; ...     % 1–6: r0, v0
           C.k_SRP_0; ...    % 7:   SRP scale
           0.0;       ...    % 8:   Doppler bias (km/s)
           st4.lat;   ...    % 9:   lat_4 (rad)
           st4.lon ];        % 10:  lon_4 (rad)

P0_est = C.P0_aug;           % 10×10

% P0_inv via Cholesky
[L0,p0] = chol(P0_est,'lower');
if p0 ~= 0, error('P0_est not PD.'); end
L0i    = L0 \ eye(size(L0));
P0_inv = L0i' * L0i;

X0_prior = X0_est;           % prior mean anchor at detection

%% 4. Measurement noise
sigma_r_DSN  = 1e-3;   % km   (1 m)
sigma_rr_DSN = 1e-7;   % km/s (0.1 mm/s)
sigma_r_ANT  = 1e-2;   % km   (10 m)
sigma_rr_ANT = 1e-6;   % km/s (1 mm/s)

Rinv_DSN = diag([1/sigma_r_DSN^2,  1/sigma_rr_DSN^2]);
Rinv_ANT = diag([1/sigma_r_ANT^2,  1/sigma_rr_ANT^2]);

%% 5. ODE options
ode_opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

%% 6. Batch iterations
max_iter   = 6;
tol_update = 1e-6;

results = struct([]);
n_state = 10;

for iter = 1:max_iter
    fprintf('\n=== PRELIM BATCH ITER %d ===\n', iter);

    X0_full = X0_est;

    % 6.1 Single propagation for all measurement times (state + STM)
    Phi0   = eye(n_state);
    X_aug0 = [X0_full; Phi0(:)];      % 10 + 100 = 110×1

    is_t0    = abs(t_sec) < 1e-9;
    t_pos    = t_sec(~is_t0);
    t_et_pos = t0_et + t_pos;
    t_et_u   = unique(t_et_pos);

    X_prop = [];
    if ~isempty(t_et_u)
        tspan = [t0_et; t_et_u(:)];
        [~, X_prop] = ode45(@(t,X) lib_dynamics(t,X,C), tspan, X_aug0, ode_opts);
    end

    Xk_all  = zeros(n_state, n_meas);
    Phi_all = zeros(n_state, n_state, n_meas);

    % t = 0 measurements
    if any(is_t0)
        idx0 = find(is_t0);
        for j = 1:numel(idx0)
            k0 = idx0(j);
            Xk_all(:,k0)    = X0_full;
            Phi_all(:,:,k0) = Phi0;
        end
    end

    % t > 0 measurements
    if ~isempty(t_et_u)
        [~, idx_u_for_meas] = ismember(t_et_pos, t_et_u);
        idx_meas_pos = find(~is_t0);

        for j = 1:numel(idx_meas_pos)
            k_meas = idx_meas_pos(j);
            iu     = idx_u_for_meas(j);

            row     = iu + 1;               % offset for t0 row
            X_aug_k = X_prop(row,:).';

            Xk_all(:,k_meas)    = X_aug_k(1:n_state);
            Phi_all(:,:,k_meas) = reshape(X_aug_k(n_state+1:end), n_state, n_state);
        end
    end

    % 6.2 Normal equations with prior
    Delta = P0_inv;
    N_vec = P0_inv * (X0_prior - X0_full);

    prefit_rr = zeros(n_meas,1);

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

    % 6.3 Solve for ΔX0
    delta_X0 = Delta \ N_vec;
    d_norm   = norm(delta_X0);
    fprintf('  ||ΔX0||_2 = %.3e\n', d_norm);

    X0_est = X0_est + delta_X0;

    % Posterior covariance via Cholesky of Delta
    [L,p] = chol(Delta,'lower');
    if p ~= 0
        warning('Delta not PD; using inv(Delta).');
        P_post = inv(Delta); %#ok<MINV>
    else
        Li     = L \ eye(size(L));
        P_post = Li' * Li;
    end

    % 6.4 Postfit Doppler (linearized)
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

    % 6.5 Store and basic diagnostics
    results(iter).X0_est     = X0_est;
    results(iter).P_post     = P_post;
    results(iter).time_days  = time_days;
    results(iter).prefit_rr  = prefit_rr;
    results(iter).postfit_rr = postfit_rr;

    rr_pre_rms  = sqrt(mean(prefit_rr.^2));
    rr_post_rms = sqrt(mean(postfit_rr.^2));
    fprintf('  RR RMS pre/post: %.3e / %.3e km/s (%.2f / %.2f mm/s)\n', ...
        rr_pre_rms, rr_post_rms, rr_pre_rms*1e6, rr_post_rms*1e6);

    % 6.6 Convergence
    if d_norm < tol_update
        fprintf('  Converged after %d iterations.\n', iter);
        break;
    end
end

X0_est_final  = X0_est;
P0_post_final = results(end).P_post;

fprintf('\n=== PRELIM 10-STATE BATCH RESULT ===\n');
fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', X0_est_final(1:3));
fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', X0_est_final(4:6));
fprintf('k_SRP:        %.6f\n', X0_est_final(7));
fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', ...
        X0_est_final(8), X0_est_final(8)*1e6);
fprintf('lat_4 (deg):  %.6f\n', X0_est_final(9)*C.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', X0_est_final(10)*C.rad2deg);

%% 7. Save
prelim_results = struct();
prelim_results.X0_prior          = X0_prior;
prelim_results.P0_prior          = P0_est;
prelim_results.X0_batch          = X0_est_final;
prelim_results.P0_batch          = P0_post_final;
prelim_results.const             = C;
prelim_results.iteration_results = results;

% Measurement / indexing info for by-station plots
prelim_results.time_days = time_days;   % days since detection
prelim_results.t_sec     = t_sec;       % seconds since detection
prelim_results.t_et      = t_et_all;    % ET (s)
prelim_results.st_id     = st_id;       % station ID per measurement
prelim_results.has_range = has_range;   % logical (range present)

save('ASTE583_PrelimBatch_Results.mat', 'prelim_results', 'meas_full');
fprintf('\nSaved prelim batch to ASTE583_PrelimBatch_Results.mat\n');
toc
