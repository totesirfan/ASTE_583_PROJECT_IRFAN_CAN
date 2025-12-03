% Phase2_Batch_Filter_7state.m
% PRELIMINARY OD (Days 0–6 after detection)
% 7-state batch LSQ on [r0(3); v0(3); Doppler bias]
% k_SRP and Station 4 location held fixed.
%
% Outputs for handout requirements:
%   - State and 1-sigma at detection epoch and at DCO (t0 + 6 days)
%   - Converged prefit/postfit Doppler residuals (mm/s vs time)
%   - 3-sigma XY position covariance ellipse at DCO

clear; clc; close all;

%% 1. Initialization
init_project();
const = lib_constants();

t_detection = cspice_str2et('2025 DEC 01 00:00:00.00');   % detection epoch
t_dco       = t_detection + 6*86400;                      % 6 days after detection

% Load measurements and keep only 0–6 days
all_obs   = Phase2_Load_Data();
valid     = ([all_obs.t] <= t_dco);
obs       = all_obs(valid);

fprintf('PRELIMINARY OD (7-state): %d measurements (Days 0–6)\n', length(obs));

% Unique times for propagation (include t0 explicitly)
t_span = unique([t_detection; [obs.t]']);

%% 2. A priori state and covariance (7 estimated params)
% Full 10-state used by dynamics: [r; v; k_SRP; bias; lat4; lon4]
X_full_nom = [const.X0_ref; ...
              1.0; ...                        % k_SRP (held fixed)
              0.0; ...                        % Doppler bias (to be estimated)
              const.stations(4).lat; ...      % Station 4 latitude (fixed here)
              const.stations(4).lon];         % Station 4 longitude (fixed here)

% Indices of estimated parameters inside 10-state
est_idx = [1:6, 8];    % r(1:3), v(4:6), bias(8)
n_est   = numel(est_idx);

% A priori sigmas for estimated sub-state
sig_r   = 100;     % km
sig_v   = 1e-3;    % km/s
sig_b   = 1e-3;    % km/s (after warm-start we expect ~1e-1, keep prior tight-ish)

sig_vec = [sig_r*ones(3,1); sig_v*ones(3,1); sig_b];
P_bar_inv = diag(1 ./ sig_vec.^2);      % 7x7 information matrix

%% 3. Warm start for Doppler bias (pre-fit residuals with nominal orbit)
fprintf('\n--- Warm Start: Doppler Bias ---\n');
Phi0        = eye(10);
X_aug_0     = [X_full_nom; Phi0(:)];
ode_opts    = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T0, X0]    = ode45(@(t,x) lib_dynamics(t,x,const), t_span, X_aug_0, ode_opts);

% Pre-fit residuals with nominal full state (bias = 0 initially)
resids_prefit = zeros(length(obs),1);
times_prefit  = zeros(length(obs),1);

for k = 1:length(obs)
    Xk = interp1(T0, X0, obs(k).t)';
    [Yc, ~] = lib_measurements(obs(k).t, Xk(1:10), obs(k).ID, const);
    if strcmp(obs(k).type, 'Doppler')
        resids_prefit(k) = obs(k).value - Yc(2);
    else
        resids_prefit(k) = NaN;  % no Doppler in this slot for days 0–6, but keep size
    end
    times_prefit(k) = (obs(k).t - t_detection)/86400;
end

% Warm-start bias estimate = mean Doppler residual (ignore NaNs)
bias_guess = mean(resids_prefit(~isnan(resids_prefit)));
fprintf('Initial Doppler Mean Residual (bias guess): %.6f km/s\n', bias_guess);

% Update full nominal state with warm-start bias
X_full_nom(8) = bias_guess;
X_bar_full    = X_full_nom;   % full-state a priori "anchor"

% Extract 7-state a priori vector
x_nom = X_full_nom(est_idx);  % [r; v; bias]

%% 4. Batch filter settings
max_iter    = 20;
conv_tol    = 1e-4;
sigma_floor = 1.0e-2;   % inflated Doppler noise (0.01 km/s) for stability

fprintf('\n--- Starting 7-Parameter Batch Filter (r, v, bias) ---\n');

for iter = 1:max_iter
    % Normal equations initialised with a priori information
    Lambda = P_bar_inv;                         % 7x7
    dx_bar = X_bar_full(est_idx) - x_nom;       % 7x1
    N_vec  = P_bar_inv * dx_bar;               % 7x1
    
    % Propagate from t_detection with current full-state nominal
    X_aug_iter = [zeros(10,1); zeros(100,1)];   % allocate
    X_aug_iter(1:10) = X_full_nom;
    X_aug_iter(11:end) = Phi0(:);               % STM(t0) = I
    [T_iter, X_iter] = ode45(@(t,x) lib_dynamics(t,x,const), t_span, X_aug_iter, ode_opts);
    
    rms_resid = 0;
    n_dopp    = 0;
    
    for k = 1:length(obs)
        tk   = obs(k).t;
        Xval = interp1(T_iter, X_iter, tk)';
        
        state_k = Xval(1:10);
        Phi_k   = reshape(Xval(11:end),10,10);  % STM(tk,t0)
        
        [Y_comp, H_10] = lib_measurements(tk, state_k, obs(k).ID, const);
        
        if strcmp(obs(k).type,'Range')
            y_c = Y_comp(1); h_t = H_10(1,:);
        else
            y_c = Y_comp(2); h_t = H_10(2,:);
        end
        
        resid = obs(k).value - y_c;
        if strcmp(obs(k).type,'Doppler')
            rms_resid = rms_resid + resid^2;
            n_dopp    = n_dopp + 1;
        end
        
        % Map measurement sensitivity back to epoch and extract 7 estimated cols
        H_k_full = h_t * Phi_k;          % 1x10
        H_k      = H_k_full(est_idx);    % 1x7
        
        % Measurement weight (inflated)
        sig_meas = max(obs(k).sigma, sigma_floor);
        w        = 1 / (sig_meas^2);
        
        Lambda = Lambda + (H_k' * w * H_k);
        N_vec  = N_vec  + (H_k' * w * resid);
    end
    
    % Solve normal equations with small ridge if needed
    if rcond(Lambda) < 1e-14
        Lambda = Lambda + 1e-6 * diag(diag(Lambda));
    end
    delta_x = Lambda \ N_vec;
    
    % Step limiting: keep updates physically small
    scl = 1.0;
    if norm(delta_x(1:3)) > 1000,  scl = min(scl, 1000/norm(delta_x(1:3))); end   % km
    if norm(delta_x(4:6)) > 0.02,  scl = min(scl, 0.02/norm(delta_x(4:6))); end   % km/s
    if abs(delta_x(7))    > 0.10,  scl = min(scl, 0.10/abs(delta_x(7)));   end    % km/s bias
    
    delta_x = scl * delta_x;
    
    % Update 7-state and embed back into full state vector
    x_nom = x_nom + delta_x;
    X_full_nom(est_idx) = x_nom;
    
    % Report iteration metrics (RMS only over Doppler)
    rms_kmps = sqrt(rms_resid / max(n_dopp,1));
    fprintf('Iteration %2d: RMS = %.6f km/s, Update Max = %.3e (scale = %.2f)\n', ...
            iter, rms_kmps, max(abs(delta_x)), scl);
    
    if max(abs(delta_x)) < conv_tol
        break;
    end
end

% Posterior covariance for 7-state
P_post_7 = inv(Lambda);
sig_7    = sqrt(diag(P_post_7));

%% 5. Final state at detection epoch and at DCO
fprintf('\n=== FINAL 7-STATE ESTIMATE (Epoch: %s) ===\n', const.epoch_utc_str);

% Position
fprintf('r0 (km):   [%12.4f, %12.4f, %12.4f]\n', x_nom(1:3));
fprintf('            +/- [%.3f, %.3f, %.3f]  (1-sigma)\n', sig_7(1:3));

% Velocity
fprintf('v0 (km/s): [%12.7f, %12.7f, %12.7f]\n', x_nom(4:6));
fprintf('            +/- [%.1e, %.1e, %.1e]  (1-sigma)\n', sig_7(4:6));

% Doppler bias
fprintf('Bias (km/s):  %12.8f  +/- %.1e\n', x_nom(7), sig_7(7));

% Fixed parameters (not estimated in this 7-state prelim OD)
fprintf('k_SRP:            1.0000  (held fixed in prelim OD)\n');
fprintf('Stn 4 Lat:      -80.0000 deg (held fixed)\n');
fprintf('Stn 4 Lon:        0.0000 deg (held fixed)\n');

% Propagate final estimate to DCO and map covariance
X_aug_final0 = [X_full_nom; Phi0(:)];
[Tf, Xf]     = ode45(@(t,x) lib_dynamics(t,x,const), [t_detection; t_dco], X_aug_final0, ode_opts);
X_dco_aug    = Xf(end,:)';
state_dco    = X_dco_aug(1:10);
Phi_dco      = reshape(X_dco_aug(11:end),10,10);

Phi_red      = Phi_dco(1:3, est_idx);      % position rows, estimated cols
P_r_dco      = Phi_red * P_post_7 * Phi_red.';   % 3x3 position cov at DCO

sig_r_dco = sqrt(diag(P_r_dco));
fprintf('\n--- Mapped State at DCO (6 days after detection) ---\n');
fprintf('r_DCO (km):   [%12.4f, %12.4f, %12.4f] +/- [%.3f, %.3f, %.3f]\n', ...
    state_dco(1:3), sig_r_dco(1), sig_r_dco(2), sig_r_dco(3));
fprintf('v_DCO (km/s): [%12.7f, %12.7f, %12.7f]\n', state_dco(4:6));

%% 6. Post-fit residuals with final solution
resids_post  = zeros(length(obs),1);
times_post   = zeros(length(obs),1);

% reuse Tf,Xf propagation (covers t_detection..t_dco densely)
for k = 1:length(obs)
    Xk = interp1(Tf, Xf, obs(k).t)';
    [Yc, ~] = lib_measurements(obs(k).t, Xk(1:10), obs(k).ID, const);
    if strcmp(obs(k).type,'Doppler')
        resids_post(k) = obs(k).value - Yc(2);
    else
        resids_post(k) = NaN;
    end
    times_post(k) = (obs(k).t - t_detection)/86400;
end

%% 7. Plots: prefit / postfit Doppler residuals & 3-sigma XY ellipse
% Doppler residuals (mm/s)
idx_dopp = ~isnan(resids_prefit);

figure('Name','Prelim OD Doppler Residuals','Color','w');
subplot(2,1,1);
plot(times_prefit(idx_dopp), resids_prefit(idx_dopp)*1e6, 'k.', 'MarkerSize', 4);
xlabel('Time since detection (days)');
ylabel('Prefit Doppler (mm/s)');
title('Batch Filter: Prefit Doppler Residuals');
grid on; yline(0,'k-');

subplot(2,1,2);
plot(times_post(idx_dopp),  resids_post(idx_dopp)*1e6, 'b.', 'MarkerSize', 4);
xlabel('Time since detection (days)');
ylabel('Postfit Doppler (mm/s)');
title('Batch Filter: Postfit Doppler Residuals');
grid on; yline(0,'k-');

% 3-sigma XY covariance ellipse at DCO
C_xy   = P_r_dco(1:2,1:2);    % 2x2 covariance in X-Y
[vec, val] = eig(C_xy);
sig_xy = sqrt(diag(val));     % 1-sigma along principal axes
theta  = linspace(0,2*pi,200);
ellipse_unit = [cos(theta); sin(theta)];
ellipse_xy   = state_dco(1:2) + 3 * (vec * (sig_xy .* ellipse_unit));

figure('Name','3-sigma XY Covariance at DCO','Color','w');
plot(state_dco(1), state_dco(2), 'rx', 'MarkerSize',8, 'LineWidth',2); hold on;
plot(ellipse_xy(1,:), ellipse_xy(2,:), 'b-');
axis equal;
grid on;
xlabel('X (km)'); ylabel('Y (km)');
title('Batch Filter: 3\sigma XY Position Covariance at DCO');
