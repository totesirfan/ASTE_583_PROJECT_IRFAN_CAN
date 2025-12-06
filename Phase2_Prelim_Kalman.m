% Phase2_Prelim_Kalman_.m
% Preliminary Sequential Kalman Filter (CKF-style) for 0–6 days of data.
% Includes DEBUG printing so you can monitor progress.

clear; clc; close all;

%% 1. Setup
init_project();
const = lib_constants();

t0_et = const.t_detect_et;
DCO_days_prelim = 6.0;
t_DCO_sec = DCO_days_prelim * const.day2sec;

%% 2. Load + restrict measurements (0–6 days)
meas_full = load_project_measurements();

use_idx = (meas_full.time_sec <= t_DCO_sec);
t_sec = meas_full.time_sec(use_idx);
st_id = meas_full.station_id(use_idx);
rng_km = meas_full.range_km(use_idx);
rr_kmps = meas_full.rr_kmps(use_idx);

% Sort by time
[t_sec, sort_idx] = sort(t_sec);
st_id  = st_id(sort_idx);
rng_km = rng_km(sort_idx);
rr_kmps = rr_kmps(sort_idx);

n_meas = numel(t_sec);
fprintf('Prelim KF: using %d measurements up to %.3f days.\n', ...
        n_meas, t_DCO_sec/const.day2sec);

time_days = t_sec/const.day2sec;
t_et_all = const.t_detect_et + t_sec;

%% 3. Reference initial state (10-state augmented)
stat4 = const.stations(4);
X0_ref = [ const.X0_ref;
           const.k_SRP_0;
           0.0;                % Doppler bias at detection
           stat4.lat;
           stat4.lon ];

n_est = 10;

% A priori covariance
P0 = const.P0_aug;

deltaX_hat_prev = zeros(n_est,1);
P_prev = P0;

%% 4. Measurement noise (handout)
sigma_r_DSN  = 1e-3;
sigma_rr_DSN = 1e-7;

sigma_r_ANT  = 1e-2;
sigma_rr_ANT = 1e-6;

%% 5. Propagate nominal trajectory + STM to all measurement times
fprintf('Propagating nominal trajectory + STM...\n');

ode_opts = odeset('RelTol',1e-11,'AbsTol',1e-11);

t_et_unique = unique(t_et_all);
tspan = [t0_et; t_et_unique(:)];

Phi0_init = eye(n_est);
X_aug0 = [X0_ref; Phi0_init(:)];

[T_prop, X_prop] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                         tspan, X_aug0, ode_opts);

[~, idx_unique_map] = ismember(t_et_unique, T_prop);

X_nom = zeros(n_est,n_meas);
Phi0  = zeros(n_est,n_est,n_meas);

for k = 1:n_meas
    iu = find(t_et_unique == t_et_all(k),1,'first');
    row = idx_unique_map(iu);
    X_aug_k = X_prop(row,:)';

    X_nom(:,k) = X_aug_k(1:n_est);
    Phi0(:,:,k) = reshape(X_aug_k(n_est+1:end),n_est,n_est);
end

fprintf('Propagation complete. Beginning KF updates...\n\n');

%% 6. Kalman updates measurement-by-measurement
prefit_rr  = zeros(n_meas,1);
postfit_rr = zeros(n_meas,1);

for k = 1:n_meas
    t_k_et = t_et_all(k);
    st     = st_id(k);

    %% --- TIME UPDATE ---
    Phi_k0 = Phi0(:,:,k);

    if k == 1
        Phi_prev0 = eye(n_est);
    else
        Phi_prev0 = Phi0(:,:,k-1);
    end

    Phi_k_prev = Phi_k0 / Phi_prev0;

    deltaX_bar = Phi_k_prev * deltaX_hat_prev;
    P_bar      = Phi_k_prev * P_prev * Phi_k_prev';

    %% --- MEASUREMENT MODEL ---
    X_nom_k = X_nom(:,k);

    [Y_nom,H_tilde] = lib_measurements(t_k_et,X_nom_k,st,const);

    rho_obs = rng_km(k);
    rr_obs  = rr_kmps(k);
    has_range = ~isnan(rho_obs);

    prefit_rr(k) = rr_obs - Y_nom(2);

    if has_range
        dy = [rho_obs - Y_nom(1);
              rr_obs  - Y_nom(2)];
        Hk = H_tilde;

        if st == 4
            Rk = diag([sigma_r_ANT^2, sigma_rr_ANT^2]);
        else
            Rk = diag([sigma_r_DSN^2, sigma_rr_DSN^2]);
        end
    else
        dy = rr_obs - Y_nom(2);
        Hk = H_tilde(2,:);

        if st == 4
            Rk = sigma_rr_ANT^2;
        else
            Rk = sigma_rr_DSN^2;
        end
    end

    %% --- MEASUREMENT UPDATE ---
    if has_range
        S = Hk*P_bar*Hk' + Rk;
        K = P_bar*Hk'/S;
        innov = dy - Hk*deltaX_bar;
        deltaX_hat_k = deltaX_bar + K*innov;
        P_k = (eye(n_est)-K*Hk)*P_bar;
    else
        S = Hk*P_bar*Hk' + Rk;
        K = (P_bar*Hk')/S;
        innov = dy - Hk*deltaX_bar;
        deltaX_hat_k = deltaX_bar + K*innov;
        P_k = (eye(n_est)-K*Hk)*P_bar;
    end

    %% --- POSTFIT ---
    X_est_k = X_nom_k + deltaX_hat_k;
    [Y_post,~] = lib_measurements(t_k_et,X_est_k,st,const);
    postfit_rr(k) = rr_obs - Y_post(2);

    %% --- STORE ---
    deltaX_hat_prev = deltaX_hat_k;
    P_prev = P_k;

    %% --- DEBUG PRINT ---
    if mod(k,5000) == 0
        fprintf('[%5d/%5d] t = %.3f d, st=%d  prefit= %+9.3e  postfit= %+9.3e\n', ...
                k, n_meas, time_days(k), st, prefit_rr(k), postfit_rr(k));
    end
end

%% 7. Map back to detection epoch
ell = n_meas;

Phi_ell0 = Phi0(:,:,ell);
Phi_0ell = eye(n_est) / Phi_ell0;

deltaX0_hat = Phi_0ell * deltaX_hat_k;
P0_KF       = Phi_0ell * P_k * Phi_0ell';

X0_KF = X0_ref + deltaX0_hat;

%% 8. Final summary
fprintf('\n=== PRELIM 10-STATE KALMAN RESULT (t0) ===\n');
fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', X0_KF(1:3));
fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', X0_KF(4:6));
fprintf('k_SRP:        %.6f\n', X0_KF(7));
fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', ...
        X0_KF(8), X0_KF(8)*1e6);
fprintf('lat_4 (deg):  %.6f\n', X0_KF(9)*const.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', X0_KF(10)*const.rad2deg);

%% 9. Save
kf_results = struct();
kf_results.X0_KF = X0_KF;
kf_results.P0_KF = P0_KF;
kf_results.prefit_rr  = prefit_rr;
kf_results.postfit_rr = postfit_rr;
kf_results.time_days  = time_days;

save('ASTE583_PrelimKalman_Results.mat','kf_results');

fprintf('\nSaved to ASTE583_PrelimKalman_Debug.mat\n');
