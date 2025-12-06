% Phase2_Prelim_Batch.m
% Preliminary OD – 10-state Batch Least Squares using radiometric data
% up to 6 days after detection.
%
% Estimated state at detection (10×1):
%   X0_est = [ r0(3); v0(3); k_SRP; rho_dot_bias; lat_4; lon_4 ]
%
% Uses:
%   - init_project.m
%   - lib_constants.m
%   - lib_dynamics.m  (10-state + 10x10 STM)
%   - lib_measurements.m
%   - load_project_measurements.m

clear; clc; close all;

%% 1. Setup & constants
init_project();
const = lib_constants();

t0_et = const.t_detect_et;             % ET at detection epoch

% Preliminary DCO at 6 days after detection
DCO_days_prelim = 6.0;
t_DCO_sec       = DCO_days_prelim * const.day2sec;

%% 2. Load measurements & restrict to prelim window
% Keep full measurement set for later phases (final batch + Kalman)
meas_full = load_project_measurements();   % full 0–14 days

% Prelim uses only 0–6 days
use_idx = (meas_full.time_sec <= t_DCO_sec);
t_sec   = meas_full.time_sec(use_idx);      % seconds since detection
st_id   = meas_full.station_id(use_idx);
rng_km  = meas_full.range_km(use_idx);
rr_kmps = meas_full.rr_kmps(use_idx);

n_meas = numel(t_sec);
fprintf('Prelim 10-state Batch OD: using %d measurements up to %.3f days.\n', ...
        n_meas, t_DCO_sec/const.day2sec);

time_days = t_sec / const.day2sec;
t_et_all  = t0_et + t_sec;                 % ET for each measurement

%% 3. Estimated state and a priori covariance (10-state)

% State order must match lib_dynamics/lib_measurements:
% [ r(3); v(3); k_SRP; bias; lat_4; lon_4 ]

stat4 = const.stations(4);

X0_est = [ const.X0_ref;    ... % 1–6: r0, v0  (Sun-centered EMO2000)
           const.k_SRP_0;   ... % 7:   SRP scale factor
           0.0;             ... % 8:   Doppler bias (km/s) during first 6 days
           stat4.lat;       ... % 9:   Antarctica latitude  (rad)
           stat4.lon ];         % 10:  Antarctica longitude (rad)

n_est = numel(X0_est);          % = 10

% Full augmented a priori covariance from lib_constants
P0_est = const.P0_aug;          % 10×10

% A priori information matrix via Cholesky: P0_inv = inv(P0_est)
[L0,p0] = chol(P0_est,'lower');
if p0 ~= 0
    error('P0_est is not positive definite.');
end
L0i    = L0 \ eye(size(L0));    % L0 * L0' = P0  => inv(P0) = (L0^{-1})' * L0^{-1}
P0_inv = L0i' * L0i;            % 10×10

% Fixed a priori mean at detection (anchor)
X0_prior = X0_est;              % 10×1

%% 4. Measurement noise (project handout)
% DSN:    σ_range = 1 m,  σ_rdot = 0.1 mm/s
sigma_r_DSN   = 1e-3;           % km
sigma_rr_DSN  = 1e-7;           % km/s

% Antarctica: σ_range = 10 m, σ_rdot = 1 mm/s
sigma_r_ANT   = 1e-2;           % km
sigma_rr_ANT  = 1e-6;           % km/s

%% 5. ODE options
ode_opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

%% 6. Batch LS iterations
max_iter   = 6;
tol_update = 1e-6;              % ||ΔX0||2 convergence tolerance

results = struct([]);

for iter = 1:max_iter
    fprintf('\n=== PRELIM 10-STATE BATCH ITERATION %d ===\n', iter);

    % Nominal 10-state at detection for this linearization
    X0_full = X0_est;           % 10×1

    % ---------- 6.1 ONE PROPAGATION FOR ALL MEASUREMENTS ----------
    Phi0   = eye(10);
    X_aug0 = [ X0_full;
               Phi0(:) ];       % 10 + 100 = 110×1

    % Separate t = 0 and t > 0 measurement times
    is_t0        = (abs(t_sec) < 1e-9);
    has_t0_meas  = any(is_t0);
    t_sec_pos    = t_sec(~is_t0);
    t_et_pos     = t0_et + t_sec_pos;
    t_et_pos_unique = unique(t_et_pos);

    X_prop = [];
    if ~isempty(t_et_pos_unique)
        % tspan: from t0 to last positive time; solver returns at given times
        tspan = [t0_et; t_et_pos_unique(:)];
        [~, X_prop] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                            tspan, X_aug0, ode_opts);
    end

    % Pre-allocate state & STM at each measurement time
    Xk_all   = zeros(10, n_meas);
    Phi_all  = zeros(10, 10, n_meas);

    % Fill for t = 0 measurements: state is X0_full, Φ = I
    if has_t0_meas
        idx0 = find(is_t0);
        for j = 1:numel(idx0)
            k0 = idx0(j);
            Xk_all(:,k0)    = X0_full;
            Phi_all(:,:,k0) = Phi0;
        end
    end

    % Fill for positive times using the single propagation
    if ~isempty(t_et_pos_unique)
        % Map each positive measurement time to its row in t_et_pos_unique
        [~, idx_unique_for_meas] = ismember(t_et_pos, t_et_pos_unique);
        idx_meas_all = find(~is_t0);   % indices in 1..n_meas that are >0

        for j = 1:numel(idx_meas_all)
            k_meas = idx_meas_all(j);          % measurement index
            iu     = idx_unique_for_meas(j);   % index into unique times (1..Nu)

            row = iu + 1;                      % row in X_prop (offset for t0)
            X_aug_k = X_prop(row,:).';

            Xk_all(:,k_meas)    = X_aug_k(1:10);
            Phi_all(:,:,k_meas) = reshape(X_aug_k(11:end), 10, 10);
        end
    end

    % ---------- 6.2 BUILD NORMAL EQUATIONS (WITH PRIOR ANCHOR) ----------
    % Lambda = Σ H^T R^{-1} H + P0_inv
    % N      = Σ H^T R^{-1} dy + P0_inv (X0_prior - X0_full)
    Delta = P0_inv;                              % 10×10
    N_vec = P0_inv * (X0_prior - X0_full);       % 10×1

    prefit_rr = zeros(n_meas,1);                 % Doppler prefit residuals

    for k = 1:n_meas

        t_k_et = t_et_all(k);
        st     = st_id(k);

        X_k   = Xk_all(:,k);
        Phi_k = Phi_all(:,:,k);

        % Measurement model at nominal X_k
        [Y_comp, H_tilde] = lib_measurements(t_k_et, X_k, st, const);
        % Y_comp = [rho; rho_dot], H_tilde = 2×10

        % Map H to initial time
        H_full = H_tilde * Phi_k;      % 2×10

        % Observations
        rho_obs   = rng_km(k);         % may be NaN
        rr_obs    = rr_kmps(k);        % always present
        has_range = ~isnan(rho_obs);

        if has_range
            % 2×1 residual: [range; range-rate]
            dy = [rho_obs - Y_comp(1);
                  rr_obs  - Y_comp(2)];
            Hk = H_full;               % 2×10

            if st == 4
                Rk = diag([sigma_r_ANT^2; sigma_rr_ANT^2]);  % 2×2
            else
                Rk = diag([sigma_r_DSN^2; sigma_rr_DSN^2]);
            end

            prefit_rr(k) = dy(2);

            % 2×2 accumulation
            Delta = Delta + Hk' * (Rk \ Hk);
            N_vec = N_vec + Hk' * (Rk \ dy);

        else
            % Doppler-only (1×1 residual)
            dy = rr_obs - Y_comp(2);   % scalar
            Hk = H_full(2,:);          % 1×10

            if st == 4
                sigma_rr = sigma_rr_ANT;
            else
                sigma_rr = sigma_rr_DSN;
            end

            Rk = sigma_rr^2;           % scalar
            prefit_rr(k) = dy;

            % Scalar accumulation
            Delta = Delta + (Hk' * Hk) / Rk;
            N_vec = N_vec + (Hk' * dy) / Rk;
        end
    end

    % ---------- 6.3 SOLVE FOR ΔX0 ----------
    % Solve Delta * ΔX0 = N_vec via backslash (no explicit inv)
    delta_X0 = Delta \ N_vec;
    fprintf('||ΔX0||_2 = %.3e\n', norm(delta_X0));

    % Update estimate
    X0_est = X0_est + delta_X0;

    % Posterior covariance via Cholesky: P_post = inv(Delta)
    [L,p] = chol(Delta,'lower');
    if p ~= 0
        warning('Delta not positive definite; computing naive inverse for P_post.');
        P_post = inv(Delta);   %#ok<MINV>
    else
        Li    = L \ eye(size(L));
        P_post = Li' * Li;
    end

    % ---------- 6.4 POSTFIT DOPPLER (LINEARIZED) ----------
    postfit_rr = zeros(n_meas,1);

    for k = 1:n_meas
        t_k_et = t_et_all(k);
        st     = st_id(k);

        X_k   = Xk_all(:,k);
        Phi_k = Phi_all(:,:,k);

        % Updated state at t_k via linear correction
        X_k_upd = X_k + Phi_k * delta_X0;

        [Y_comp_upd, ~] = lib_measurements(t_k_et, X_k_upd, st, const);

        rr_obs = rr_kmps(k);
        postfit_rr(k) = rr_obs - Y_comp_upd(2);
    end

    % ---------- 6.5 STORE & PLOT ----------
    results(iter).X0_est     = X0_est;
    results(iter).P_post     = P_post;
    results(iter).time_days  = time_days;
    results(iter).prefit_rr  = prefit_rr;
    results(iter).postfit_rr = postfit_rr;

    figure('Name', sprintf('Prelim 10-state Batch – Iteration %d', iter), ...
           'Color', 'w');

    % km/s
    subplot(2,1,1); hold on; grid on;
    plot(time_days, prefit_rr,  '.', 'MarkerSize', 4);
    plot(time_days, postfit_rr, '.', 'MarkerSize', 4);
    xlabel('Time since detection (days)');
    ylabel('\rhȯ residual (km/s)');
    title(sprintf('Iteration %d: Range-Rate Residuals (km/s)', iter));
    legend('Prefit','Postfit','Location','best');

    % mm/s
    subplot(2,1,2); hold on; grid on;
    plot(time_days, prefit_rr*1e6,  '.', 'MarkerSize', 4);
    plot(time_days, postfit_rr*1e6, '.', 'MarkerSize', 4);
    xlabel('Time since detection (days)');
    ylabel('\rhȯ residual (mm/s)');
    title(sprintf('Iteration %d: Range-Rate Residuals (mm/s)', iter));
    legend('Prefit','Postfit','Location','best');

    drawnow;

    % ---------- 6.6 CONVERGENCE CHECK ----------
    if norm(delta_X0) < tol_update
        fprintf('Converged after %d iterations.\n', iter);
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
fprintf('lat_4 (deg):  %.6f\n', X0_est_final(9)*const.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', X0_est_final(10)*const.rad2deg);

%% 7. SAVE RESULTS FOR FINAL BATCH OD & KALMAN

prelim_results = struct();
prelim_results.X0_prior   = X0_prior;      % original a priori mean at detection
prelim_results.P0_prior   = P0_est;        % original a priori covariance
prelim_results.X0_batch   = X0_est_final;  % prelim batch posterior mean at detection
prelim_results.P0_batch   = P0_post_final; % prelim batch posterior covariance
prelim_results.const      = const;         % snapshot of constants (for convenience)
prelim_results.iteration_results = results; % residuals, etc. (optional but nice)

% Save full measurement set as well for later phases (0–14 days)
% meas_full already contains time_sec, station_id, range_km, rr_kmps, etc.

save('ASTE583_PrelimBatch_Results.mat', ...
     'prelim_results', 'meas_full');

fprintf('\nSaved prelim batch & measurement data to ASTE583_PrelimBatch_Results.mat\n');
