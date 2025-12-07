% Phase2_Prelim_Kalman.m
% Preliminary 10-state CKF over 0–6 days.
% State: [r0(3); v0(3); k_SRP; rho_dot_bias; lat_4; lon_4]

clear; clc; close all;

%% 1. Setup
init_project();
C      = lib_constants();
t0_et  = C.t_detect_et;
DCO_d  = 6.0;
t_DCO  = DCO_d * C.day2sec;

%% 2. Measurements (0–6 days)
meas   = load_project_measurements();

use    = meas.time_sec <= t_DCO;
t_sec  = meas.time_sec(use);
st_id  = meas.station_id(use);
rng_km = meas.range_km(use);
rr_kmps= meas.rr_kmps(use);

% Sort by time
[t_sec, ksort] = sort(t_sec);
st_id   = st_id(ksort);
rng_km  = rng_km(ksort);
rr_kmps = rr_kmps(ksort);

n_meas   = numel(t_sec);
time_days= t_sec / C.day2sec;
t_et_all = t0_et + t_sec;

fprintf('Prelim KF: using %d measurements up to %.3f days.\n', ...
        n_meas, DCO_d);

%% 3. A priori state (10-state)
st4   = C.stations(4);
n_x   = 10;

X0_ref = [ C.X0_ref; ...
           C.k_SRP_0; ...
           0.0; ...       % Doppler bias
           st4.lat; ...
           st4.lon ];

P0           = C.P0_aug;
deltaX_prev  = zeros(n_x,1);
P_prev       = P0;

%% 4. Measurement noise
sigma_r_DSN  = 1e-3;   % km
sigma_rr_DSN = 1e-7;   % km/s
sigma_r_ANT  = 1e-2;   % km
sigma_rr_ANT = 1e-6;   % km/s

%% 5. Nominal trajectory + STM to all measurement times
fprintf('Propagating nominal trajectory + STM...\n');
ode_opts = odeset('RelTol',1e-11,'AbsTol',1e-11);

t_unique = unique(t_et_all);
tspan    = [t0_et; t_unique(:)];

Phi0_init = eye(n_x);
X_aug0    = [X0_ref; Phi0_init(:)];

[T_prop, X_prop] = ode45(@(t,X) lib_dynamics(t,X,C), ...
                         tspan, X_aug0, ode_opts);

% Map each measurement time directly to a row in T_prop
[~, idx_meas_rows] = ismember(t_et_all, T_prop);

X_nom = zeros(n_x,  n_meas);
Phi0  = zeros(n_x, n_x, n_meas);

for k = 1:n_meas
    row     = idx_meas_rows(k);
    X_aug_k = X_prop(row,:).';
    X_nom(:,k)    = X_aug_k(1:n_x);
    Phi0(:,:,k)   = reshape(X_aug_k(n_x+1:end), n_x, n_x);
end

fprintf('Propagation complete. Beginning KF updates...\n\n');

%% 6. Sequential CKF updates
I     = eye(n_x);
has_r = ~isnan(rng_km);   % <- also saved later for plotting

prefit_rr  = zeros(n_meas,1);
postfit_rr = zeros(n_meas,1);

for k = 1:n_meas
    t_k = t_et_all(k);
    st  = st_id(k);

    % --- TIME UPDATE ---
    Phi_k0 = Phi0(:,:,k);
    if k == 1
        Phi_prev0 = I;
    else
        Phi_prev0 = Phi0(:,:,k-1);
    end
    Phi_k_prev = Phi_k0 / Phi_prev0;        % Φ(k,prev)

    deltaX_bar = Phi_k_prev * deltaX_prev;
    P_bar      = Phi_k_prev * P_prev * Phi_k_prev';

    % --- MEAS MODEL @ nominal ---
    X_nom_k = X_nom(:,k);
    [Y_nom,H_tilde] = lib_measurements(t_k, X_nom_k, st, C);

    rho_obs = rng_km(k);
    rr_obs  = rr_kmps(k);

    prefit_rr(k) = rr_obs - Y_nom(2);       % km/s

    if has_r(k)
        % 2x1 measurement
        dy = [rho_obs - Y_nom(1);
              rr_obs  - Y_nom(2)];
        Hk = H_tilde;

        if st == 4
            Rk = diag([sigma_r_ANT^2, sigma_rr_ANT^2]);
        else
            Rk = diag([sigma_r_DSN^2, sigma_rr_DSN^2]);
        end
    else
        % 1x1 Doppler-only
        dy = rr_obs - Y_nom(2);
        Hk = H_tilde(2,:);

        if st == 4
            Rk = sigma_rr_ANT^2;
        else
            Rk = sigma_rr_DSN^2;
        end
    end

    % --- MEAS UPDATE ---
    S     = Hk * P_bar * Hk' + Rk;
    K     = (P_bar * Hk') / S;
    innov = dy - Hk * deltaX_bar;

    deltaX_k = deltaX_bar + K * innov;
    P_k      = (I - K * Hk) * P_bar;

    % --- POSTFIT ---
    X_est_k = X_nom_k + deltaX_k;
    [Y_post,~] = lib_measurements(t_k, X_est_k, st, C);
    postfit_rr(k) = rr_obs - Y_post(2);

    % --- STORE ---
    deltaX_prev = deltaX_k;
    P_prev      = P_k;

    if mod(k,10000) == 0
        fprintf('[%5d/%5d] t=%.3f d, st=%d  pre=%+9.3e  post=%+9.3e\n', ...
            k, n_meas, time_days(k), st, prefit_rr(k), postfit_rr(k));
    end
end

% RR RMS debug
rr_pre_rms  = sqrt(mean(prefit_rr.^2));
rr_post_rms = sqrt(mean(postfit_rr.^2));
fprintf('RR RMS pre/post: %.3e / %.3e km/s (%.2f / %.2f mm/s)\n', ...
    rr_pre_rms, rr_post_rms, rr_pre_rms*1e6, rr_post_rms*1e6);

%% 7. Map estimate back to detection epoch
ell     = n_meas;
Phi_ell0= Phi0(:,:,ell);
Phi_0ell= I / Phi_ell0;

deltaX0 = Phi_0ell * deltaX_k;
P0_KF   = Phi_0ell * P_k * Phi_0ell';
X0_KF   = X0_ref + deltaX0;

%% 8. Final summary
fprintf('\n=== PRELIM 10-STATE KALMAN RESULT (t0) ===\n');
fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', X0_KF(1:3));
fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', X0_KF(4:6));
fprintf('k_SRP:        %.6f\n', X0_KF(7));
fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', ...
        X0_KF(8), X0_KF(8)*1e6);
fprintf('lat_4 (deg):  %.6f\n', X0_KF(9)*C.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', X0_KF(10)*C.rad2deg);

%% 9. Save
kf_results = struct();
kf_results.X0_KF      = X0_KF;
kf_results.P0_KF      = P0_KF;

% Global residuals (all stations)
kf_results.prefit_rr  = prefit_rr;          % km/s
kf_results.postfit_rr = postfit_rr;         % km/s

% Timing + station info for by-station plots later
kf_results.time_days  = time_days;          % days since detection
kf_results.t_sec      = t_sec;              % seconds since detection
kf_results.t_et       = t_et_all;           % ET (s)
kf_results.st_id      = st_id;              % station ID per measurement
kf_results.has_range  = has_r;              % logical (range present)

save('ASTE583_PrelimKalman_Results.mat','kf_results');
fprintf('\nSaved to ASTE583_PrelimKalman_Results.mat\n');
