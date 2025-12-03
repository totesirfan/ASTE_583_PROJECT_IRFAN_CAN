% Phase2_Kalman_7state.m
% PRELIMINARY OD (Days 0–6 after detection)
% 7-state Extended Kalman Filter on [r0(3); v0(3); Doppler bias]
% k_SRP and Station 4 location held fixed.

clear; clc; close all;

%% 1. Initialization
init_project();
const = lib_constants();

t_detection = cspice_str2et('2025 DEC 01 00:00:00.00');  % detection epoch t0
t_dco       = t_detection + 6*86400;                    % 6 days after detection

% Load measurements and keep only 0–6 days
all_obs = Phase2_Load_Data();
valid   = ([all_obs.t] <= t_dco);
obs     = all_obs(valid);

fprintf('KALMAN PRELIM OD (7-state): %d measurements (Days 0–6)\n', length(obs));

% Ensure measurements are time-ordered
[~, idx_sort] = sort([obs.t]);
obs = obs(idx_sort);

%% 2. State definition and a priori (7 estimated params)
% Full 10-state used by dynamics: [r; v; k_SRP; bias; lat4; lon4]
X_full_0 = [const.X0_ref; ...
            1.0; ...                        % k_SRP (held fixed)
            0.0; ...                        % Doppler bias (to be estimated)
            const.stations(4).lat; ...      % Station 4 latitude (fixed)
            const.stations(4).lon];         % Station 4 longitude (fixed)

% 7-state index mapping into 10-state
est_idx = [1:6, 8];     % r(1:3), v(4:6), bias(8)
n_est   = numel(est_idx);

% A priori sigmas for estimated 7-state at t0
sig_r = 100;     % km
sig_v = 1e-3;    % km/s
sig_b = 1e-3;    % km/s

sig_vec = [sig_r*ones(3,1); sig_v*ones(3,1); sig_b];
P0      = diag(sig_vec.^2);     % 7x7 prior covariance

% Initial 7-state at t0 from nominal full state
x0 = X_full_0(est_idx);         % [r; v; bias] at detection

%% 3. Warm start for bias: prefit Doppler residuals with nominal orbit
fprintf('\n--- Warm Start: Doppler Bias ---\n');
Phi0     = eye(10);
ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Propagate nominal full 10-state + STM from t0 across all obs
t_span = unique([t_detection; [obs.t]']);
X_aug0 = [X_full_0; Phi0(:)];
[T_nom, X_nom] = ode45(@(t,x) lib_dynamics(t,x,const), t_span, X_aug0, ode_opts);

resids_prefit = zeros(length(obs),1);
times_prefit  = zeros(length(obs),1);

for k = 1:length(obs)
    Xk = interp1(T_nom, X_nom, obs(k).t)';
    [Yc, ~] = lib_measurements(obs(k).t, Xk(1:10), obs(k).ID, const);
    if strcmp(obs(k).type,'Doppler')
        resids_prefit(k) = obs(k).value - Yc(2);
    else
        resids_prefit(k) = NaN;
    end
    times_prefit(k) = (obs(k).t - t_detection)/86400;
end

bias_guess = mean(resids_prefit(~isnan(resids_prefit)));
fprintf('Initial Doppler Mean Residual (bias guess): %.6f km/s\n', bias_guess);

% Update initial full state and 7-state with warm-start bias
X_full_0(8) = bias_guess;
x0(end)     = bias_guess;

%% 4. Kalman filter configuration
x_hat = x0;    % 7x1, a priori mean
P_hat = P0;    % 7x7, a priori covariance

X_full = X_full_0; 

Q_proc = zeros(n_est);        % no process noise for prelim OD
sigma_floor = 1.0e-2;         % 0.01 km/s floor for Doppler

resids_kf_prefit  = zeros(length(obs),1);
resids_kf_postfit = zeros(length(obs),1);
times_kf          = zeros(length(obs),1);

t_prev = t_detection;

fprintf('\n--- Running 7-State Extended Kalman Filter ---\n');

%% 5. Sequential EKF over all measurements
for k = 1:length(obs)
    tk = obs(k).t;
    dt = tk - t_prev;
    
    % 5a. TIME UPDATE: handle zero or positive dt
    if dt > 0
        % Propagate from t_prev to tk
        t_prop = [t_prev, tk];
        X_aug  = [X_full; Phi0(:)];  % full state + STM reset to I at t_prev
        [Tprop, Xprop] = ode45(@(t,x) lib_dynamics(t,x,const), t_prop, X_aug, ode_opts);
        
        X_end   = Xprop(end,:)';
        X_full  = X_end(1:10);                    % propagated full 10-state
        Phi_10  = reshape(X_end(11:end),10,10);   % STM(tk, t_prev)
    else
        % No propagation if tk == t_prev (multiple meas at same epoch)
        Phi_10 = eye(10);
        % X_full already represents the state at this epoch
    end
    
    % Map to 7-state
    F_7  = Phi_10(est_idx, est_idx);       % 7x7 STM for [r; v; bias]
    x_pr = X_full(est_idx);               % predicted 7-state at tk
    P_pr = F_7 * P_hat * F_7.' + Q_proc;  % predicted covariance
    
    % 5b. MEASUREMENT UPDATE
    [Yc, H_10] = lib_measurements(tk, X_full, obs(k).ID, const);
    
    if strcmp(obs(k).type,'Range')
        z_k  = obs(k).value;
        yhat = Yc(1);
        Hrow = H_10(1,:);
    else
        z_k  = obs(k).value;
        yhat = Yc(2);
        Hrow = H_10(2,:);
    end
    
    innov = z_k - yhat;                       % scalar residual (prefit)
    resids_kf_prefit(k) = innov;
    times_kf(k)         = (tk - t_detection)/86400;
    
    % Measurement sensitivity wrt 7-state
    H_k = Hrow(est_idx);                      % 1x7
    
    % Measurement noise (with inflation)
    sig_meas = max(obs(k).sigma, sigma_floor);
    R_k      = sig_meas^2;                    % scalar
    
    % Kalman gain and covariance update
    S_k = H_k * P_pr * H_k.' + R_k;           % scalar innovation variance
    K_k = (P_pr * H_k.') / S_k;               % 7x1
    
    x_hat_new = x_pr + K_k * innov;
    I7        = eye(n_est);
    P_hat_new = (I7 - K_k*H_k) * P_pr * (I7 - K_k*H_k).' + K_k*R_k*K_k.';
    
    % Postfit residual (innovation after update)
    resids_kf_postfit(k) = innov - H_k * (x_hat_new - x_pr);
    
    % Commit update
    x_hat = x_hat_new;
    P_hat = P_hat_new;
    X_full(est_idx) = x_hat;
    
    t_prev = tk;
    
    fprintf('Meas %5d: t = %.3f days, prefit = %+8.3e km/s, |K|max = %.3e\n', ...
        k, times_kf(k), innov, max(abs(K_k)));
end

%% 6. Results at detection epoch and DCO

fprintf('\n=== A PRIORI AT DETECTION EPOCH (t0) ===\n');
fprintf('r0_prior (km):   [%12.4f, %12.4f, %12.4f]\n', X_full_0(1:3));
fprintf('                 +/- [%.3f, %.3f, %.3f]  (1-sigma)\n', sig_vec(1:3));
fprintf('v0_prior (km/s): [%12.7f, %12.7f, %12.7f]\n', X_full_0(4:6));
fprintf('                 +/- [%.1e, %.1e, %.1e]  (1-sigma)\n', sig_vec(4:6));
fprintf('bias_prior (km/s):  %12.8f  +/- %.1e\n', X_full_0(8), sig_vec(7));

fprintf('\n=== FILTERED 7-STATE AT LAST MEASUREMENT (t_last ≈ 6 days) ===\n');
fprintf('r_last (km):   [%12.4f, %12.4f, %12.4f]\n', x_hat(1:3));
fprintf('               +/- [%.3f, %.3f, %.3f]  (1-sigma)\n', sqrt(diag(P_hat(1:3,1:3))));
fprintf('v_last (km/s): [%12.7f, %12.7f, %12.7f]\n', x_hat(4:6));
fprintf('               +/- [%.1e, %.1e, %.1e]  (1-sigma)\n', sqrt(diag(P_hat(4:6,4:6))));
fprintf('Bias_last (km/s):  %12.8f  +/- %.1e\n', x_hat(7), sqrt(P_hat(7,7)));

%% 7. Propagate filtered state to DCO and map covariance
fprintf('\n--- Propagating Filtered State to DCO ---\n');
t_last = obs(end).t;

X_aug_last = [X_full; Phi0(:)];
[T_d, X_d] = ode45(@(t,x) lib_dynamics(t,x,const), [t_last; t_dco], X_aug_last, ode_opts);
X_dco_aug  = X_d(end,:)';
state_dco  = X_dco_aug(1:10);
Phi_10_dco = reshape(X_dco_aug(11:end),10,10);

F_dco    = Phi_10_dco(est_idx, est_idx);
P_dco_7  = F_dco * P_hat * F_dco.';
P_r_dco  = P_dco_7(1:3,1:3);
sig_r_dc = sqrt(diag(P_r_dco));

fprintf('\n=== FILTERED STATE AT DCO (t0 + 6 days) ===\n');
fprintf('r_DCO (km):   [%12.4f, %12.4f, %12.4f]\n', state_dco(1:3));
fprintf('              +/- [%.3f, %.3f, %.3f]  (1-sigma)\n', sig_r_dc(1), sig_r_dc(2), sig_r_dc(3));
fprintf('v_DCO (km/s): [%12.7f, %12.7f, %12.7f]\n', state_dco(4:6));

%% 8. Doppler residual plots (prefit / postfit)
idx_dopp = ~isnan(resids_prefit);

figure('Name','Kalman OD Doppler Residuals','Color','w');
subplot(2,1,1);
plot(times_prefit(idx_dopp), resids_prefit(idx_dopp)*1e6, 'k.', 'MarkerSize', 4);
xlabel('Time since detection (days)');
ylabel('Prefit Doppler (mm/s)');
title('Kalman Filter: Prefit Doppler Residuals');
grid on; yline(0,'k-');

subplot(2,1,2);
plot(times_kf(idx_dopp), resids_kf_postfit(idx_dopp)*1e6, 'b.', 'MarkerSize', 4);
xlabel('Time since detection (days)');
ylabel('Postfit Doppler (mm/s)');
title('Kalman Filter: Postfit Doppler Residuals');
grid on; yline(0,'k-');

%% 9. 3-sigma XY covariance ellipse at DCO
C_xy   = P_r_dco(1:2,1:2);
[vec, val] = eig(C_xy);
sig_xy = sqrt(diag(val));

theta        = linspace(0,2*pi,200);
ellipse_unit = [cos(theta); sin(theta)];
ellipse_xy   = state_dco(1:2) + 3 * (vec * (sig_xy .* ellipse_unit));

figure('Name','Kalman 3-sigma XY Covariance at DCO','Color','w');
plot(state_dco(1), state_dco(2), 'rx', 'MarkerSize',8, 'LineWidth',2); hold on;
plot(ellipse_xy(1,:), ellipse_xy(2,:), 'b-');
axis equal;
grid on;
xlabel('X (km)'); ylabel('Y (km)');
title('Kalman Filter: 3\sigma XY Position Covariance at DCO');
