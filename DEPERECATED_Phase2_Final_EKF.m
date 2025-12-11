function Phase2_Final_EKF(prior_source)
% FINAL 10-STATE EKF OD (0–14 DAYS, CURRENT-STATE FILTER)
% prior_source: 'kalman' | 'batch' | 'handout' (default)

tic
if nargin < 1, prior_source = 'kalman'; end
prior_source = lower(strtrim(prior_source));

clearvars -except prior_source; clc; close all;
fprintf('=== PHASE 2: FINAL 10-STATE EKF (0–14 DAYS) ===\n');

%% 1. Setup & a priori
try init_project(); catch, end
C       = lib_constants();
t0_et   = C.t_detect_et;
N       = 10;
I10     = eye(N);

% Handout a priori
X0_ref = [C.X0_ref; C.k_SRP_0; 0.0; ...
          C.stations(4).lat; C.stations(4).lon];
P0_ref = C.P0_aug;
used_prior = 'Handout a priori';

switch prior_source
    case 'kalman'
        if exist('ASTE583_PrelimKalman_Results.mat','file')
            S = load('ASTE583_PrelimKalman_Results.mat');
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
            used_prior = 'Prelim Kalman solution';
        else
            warning('Prelim Kalman file not found, using handout.');
        end

    case 'batch'
        if exist('ASTE583_PrelimBatch_Results.mat','file')
            S = load('ASTE583_PrelimBatch_Results.mat');
            if isfield(S,'prelim_results')
                pr = S.prelim_results;
                if isfield(pr,'X0_batch'), X0_ref = pr.X0_batch(:); end
                if isfield(pr,'P0_batch'), P0_ref = pr.P0_batch;    end
            else
                if isfield(S,'X0_batch'), X0_ref = S.X0_batch(:); end
                if isfield(S,'P0_batch'), P0_ref = S.P0_batch;    end
            end
            used_prior = 'Prelim Batch solution';
        else
            warning('Prelim Batch file not found, using handout.');
        end

    case 'handout'
        used_prior = 'Handout a priori';

    otherwise
        warning('Unknown prior_source "%s". Using handout.', prior_source);
end

if numel(X0_ref) ~= N
    error('A priori state must be 10x1; got %d elements.', numel(X0_ref));
end
fprintf('A priori: %s\n', used_prior);

%% 2. Measurements (0–14 days)
fprintf('Loading measurements...\n');
meas      = load_project_measurements();
t_sec     = meas.time_sec(:);
st_id     = meas.station_id(:);
rng_km    = meas.range_km(:);
rr_kmps   = meas.rr_kmps(:);
has_range = ~isnan(rng_km);

t_days  = t_sec / C.day2sec;
t_et    = t0_et + t_sec;
use_idx = (t_days >= 0.0) & (t_days <= 14.0);

t_et      = t_et(use_idx);
st_id     = st_id(use_idx);
rng_km    = rng_km(use_idx);
rr_kmps   = rr_kmps(use_idx);
has_range = has_range(use_idx);

[t_et, idx] = sort(t_et);
st_id       = st_id(idx);
rng_km      = rng_km(idx);
rr_kmps     = rr_kmps(idx);
has_range   = has_range(idx);

n_meas   = numel(t_et);
t_last   = t_et(end);
t_last_d = (t_last - t0_et)/C.day2sec;
fprintf('EKF: %d measurements in [0,%.3f] days.\n', n_meas, t_last_d);

%% 3. Measurement noise (km, km/s)
sig_r_DSN  = 1e-3;
sig_rr_DSN = 1e-7;
sig_r_ANT  = 1e-2;
sig_rr_ANT = 1e-6;

%% 4. EKF settings
opt = odeset('RelTol',1e-10,'AbsTol',1e-9);

rho_prefit  = NaN(n_meas,1);
rho_postfit = NaN(n_meas,1);
rr_prefit   = NaN(n_meas,1);
rr_postfit  = NaN(n_meas,1);

%% 5. EKF forward pass (state at measurement epoch)
Xk   = X0_ref(:);
Pk   = P0_ref;
t_prev = t0_et;

for k = 1:n_meas
    tk   = t_et(k);
    sid  = st_id(k);
    useR = has_range(k);
    rhoo = rng_km(k);
    rro  = rr_kmps(k);

    % 5.1 Time update: propagate [X; Phi] from t_prev to tk
    dt = tk - t_prev;
    if dt > 0
        Z0 = [Xk; I10(:)];
        [~,Z] = ode45(@(t,Z) lib_dynamics(t,Z,C), [t_prev tk], Z0, opt);
        Z_end = Z(end,:).';
        Xk    = Z_end(1:N);
        Phi   = reshape(Z_end(N+1:end), N, N);
        Pk    = Phi * Pk * Phi.';   % Q = 0
        t_prev = tk;
    end

    % 5.2 Measurement model at current state
    [Ynom,Htilde] = lib_measurements(tk, Xk, sid, C);

    if useR
        % 2x1 measurement (rho, rr)
        y_obs = [rhoo; rro];
        dy    = y_obs - Ynom;      % prefit
        Hk    = Htilde;            % 2x10

        if sid == 4
            Rk = diag([sig_r_ANT^2; sig_rr_ANT^2]);
        else
            Rk = diag([sig_r_DSN^2; sig_rr_DSN^2]);
        end

        rho_prefit(k) = dy(1);
        rr_prefit(k)  = dy(2);

        S_k = Hk * Pk * Hk.' + Rk;
        K   = (Pk * Hk.') / S_k;

        Xk = Xk + K * dy;
        Pk = (I10 - K*Hk) * Pk;

        % Postfit: recompute h(x) at updated state
        Y_post = lib_measurements(tk, Xk, sid, C);
        dy_post = y_obs - Y_post;
        rho_postfit(k) = dy_post(1);
        rr_postfit(k)  = dy_post(2);

    else
        % Range-rate only (1x1)
        dy  = rro - Ynom(2);
        Hk  = Htilde(2,:);   % 1x10

        if sid == 4
            Rk = sig_rr_ANT^2;
        else
            Rk = sig_rr_DSN^2;
        end

        rr_prefit(k) = dy;

        S_k = Hk * (Pk * Hk.') + Rk;   % scalar
        K   = (Pk * Hk.') / S_k;

        Xk = Xk + K * dy;
        Pk = (I10 - K*Hk) * Pk;

        Y_post = lib_measurements(tk, Xk, sid, C);
        dy_post = rro - Y_post(2);
        rr_postfit(k) = dy_post;
    end
end

%% 6. Back-propagate solution to t0
fprintf('Back-propagating EKF solution to detection epoch...\n');
Z0_back = [Xk; I10(:)];
[~,Zb]  = ode45(@(t,Z) lib_dynamics(t,Z,C), [t_prev t0_et], Z0_back, opt);
Zb_end  = Zb(end,:).';
X0_EKF  = Zb_end(1:N);
Phi_b   = reshape(Zb_end(N+1:end), N, N);
P0_EKF  = Phi_b * Pk * Phi_b.';

%% 7. Report
r0 = X0_EKF(1:3);
v0 = X0_EKF(4:6);
k_SRP = X0_EKF(7);
bias  = X0_EKF(8);
lat4  = X0_EKF(9);
lon4  = X0_EKF(10);

rr_rms_mm = sqrt(nanmean(rr_prefit.^2)) * 1e6;

fprintf('\nFINAL EKF (mapped to t0):\n');
fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', r0);
fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', v0);
fprintf('k_SRP:        %.6f\n', k_SRP);
fprintf('rho_dot_bias: %.3f mm/s\n', bias*1e6);
fprintf('lat_4 (deg):  %.6f\n', lat4 * C.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', lon4 * C.rad2deg);
fprintf('RR prefit RMS: %.2f mm/s\n', rr_rms_mm);

%% 8. Save
X0_final_EKF = X0_EKF;
P0_final_EKF = P0_EKF;

save('ASTE583_FinalEKF_Results.mat', ...
     'X0_final_EKF','P0_final_EKF', ...
     't_et','st_id','has_range', ...
     'rho_prefit','rho_postfit', ...
     'rr_prefit','rr_postfit');

fprintf('Saved EKF results to ASTE583_FinalEKF_Results.mat\n');
toc
end
