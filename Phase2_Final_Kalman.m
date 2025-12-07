function Phase2_Final_Kalman(prior_source)
% FINAL 10-STATE CKF OD (0–14 DAYS, ITERATED AT t0)
% prior_source: 'kalman' | 'batch' | 'handout' (default)

tic
if nargin < 1, prior_source = 'handout'; end
prior_source = lower(strtrim(prior_source));

clearvars -except prior_source; clc; close all;
fprintf('=== PHASE 2: FINAL 10-STATE CKF (0–14 DAYS) ===\n');

%% 1. Setup & a priori
try init_project(); catch, end
C       = lib_constants();
t0_et   = C.t_detect_et;
N       = 10;

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

t_et    = t_et(use_idx);
st_id   = st_id(use_idx);
rng_km  = rng_km(use_idx);
rr_kmps = rr_kmps(use_idx);
has_range = has_range(use_idx);

[t_et, idx] = sort(t_et);
st_id       = st_id(idx);
rng_km      = rng_km(idx);
rr_kmps     = rr_kmps(idx);
has_range   = has_range(idx);

n_meas  = numel(t_et);
t_last  = t_et(end);
t_last_d = (t_last - t0_et)/C.day2sec;
fprintf('Using %d measurements up to %.3f days.\n', n_meas, t_last_d);

%% 3. Measurement noise (km, km/s)
sig_r_DSN  = 1e-3;
sig_rr_DSN = 1e-7;
sig_r_ANT  = 1e-2;
sig_rr_ANT = 1e-6;

%% 4. CKF settings
opt   = odeset('RelTol',1e-10,'AbsTol',1e-9);
max_it = 9;
tol_dx = 1e-3;

rho_prefit  = NaN(n_meas,1);
rho_postfit = NaN(n_meas,1);
rr_prefit   = NaN(n_meas,1);
rr_postfit  = NaN(n_meas,1);

%% 5. Iterated CKF at t0
X0_curr = X0_ref(:);
P0_curr = P0_ref;
I10     = eye(N);

for it = 1:max_it
    fprintf('\n=== CKF ITERATION %d ===\n', it);

    % 5.1 Propagate nominal + STM t0 -> t_last
    Phi0 = I10;
    Z0   = [X0_curr; Phi0(:)];
    fprintf('  Propagating to %.3f days...\n', t_last_d);
    [Tprop,Zprop] = ode45(@(t,Z) lib_dynamics(t,Z,C), [t0_et t_last], Z0, opt);
    X_hist   = Zprop(:,1:N);
    Phi_hist = Zprop(:,N+1:end);
    fprintf('  Prop done, %d steps. Interpolating...\n', numel(Tprop));

    X_meas   = interp1(Tprop, X_hist,   t_et, 'spline');
    Phi_meas = interp1(Tprop, Phi_hist, t_et, 'spline');

    % 5.2 Sequential CKF at t0
    deltaX0 = zeros(N,1);
    P0      = P0_curr;

    rho_pre_it  = NaN(n_meas,1);
    rho_post_it = NaN(n_meas,1);
    rr_pre_it   = zeros(n_meas,1);
    rr_post_it  = zeros(n_meas,1);

    for k = 1:n_meas
        tk   = t_et(k);
        sid  = st_id(k);
        useR = has_range(k);
        rhoo = rng_km(k);
        rro  = rr_kmps(k);

        Xk    = X_meas(k,:).';
        Phik  = reshape(Phi_meas(k,:).', N, N);
        [Ynom,Htilde] = lib_measurements(tk, Xk, sid, C);
        Hfull = Htilde * Phik;  % 2x10

        if useR
            y_obs = [rhoo; rro];
            dy    = y_obs - Ynom;
            Hk    = Hfull;

            if sid == 4
                Rk = diag([sig_r_ANT^2; sig_rr_ANT^2]);
            else
                Rk = diag([sig_r_DSN^2; sig_rr_DSN^2]);
            end

            rho_pre_it(k) = dy(1);
            rr_pre_it(k)  = dy(2);

            innov = dy - Hk*deltaX0;
            S_k   = Hk * P0 * Hk' + Rk;
            K0    = (P0 * Hk') / S_k;

            deltaX0 = deltaX0 + K0 * innov;
            P0      = (I10 - K0*Hk) * P0;

            dy_post = dy - Hk*deltaX0;
            rho_post_it(k) = dy_post(1);
            rr_post_it(k)  = dy_post(2);

        else
            dy = rro - Ynom(2);
            Hk = Hfull(2,:);    % 1x10

            if sid == 4
                Rk = sig_rr_ANT^2;
            else
                Rk = sig_rr_DSN^2;
            end

            rr_pre_it(k) = dy;

            innov = dy - Hk*deltaX0;
            S_k   = Hk * (P0 * Hk') + Rk;
            K0    = (P0 * Hk') / S_k;

            deltaX0 = deltaX0 + K0 * innov;
            P0      = (I10 - K0*Hk) * P0;

            dy_post = dy - Hk*deltaX0;
            rr_post_it(k) = dy_post;
        end
    end

    % 5.3 Update reference at t0 and check convergence
    dX_norm = norm(deltaX0);
    idx_rho = ~isnan(rho_pre_it);

    rr_pre_rms_mm  = sqrt(mean(rr_pre_it.^2))   * 1e6;
    rr_post_rms_mm = sqrt(mean(rr_post_it.^2))  * 1e6;
    rho_pre_rms_km = sqrt(mean(rho_pre_it(idx_rho).^2));
    rho_post_rms_km= sqrt(mean(rho_post_it(idx_rho).^2));

    fprintf('  ||ΔX0|| = %.3e\n', dX_norm);
    fprintf('  RR RMS  pre/post ≈ %.2f / %.2f mm/s\n', ...
            rr_pre_rms_mm, rr_post_rms_mm);
    fprintf('  R  RMS  pre/post ≈ %.3f / %.3f km\n', ...
            rho_pre_rms_km, rho_post_rms_km);

    X0_curr = X0_curr + deltaX0;
    P0_curr = P0;

    rho_prefit  = rho_pre_it;
    rho_postfit = rho_post_it;
    rr_prefit   = rr_pre_it;
    rr_postfit  = rr_post_it;

    if dX_norm < tol_dx
        fprintf('  Converged: ||ΔX0|| < %.1e, stopping.\n', tol_dx);
        break;
    end
end

%% 6. Final solution at t0
X0_final = X0_curr;
P0_final = P0_curr;

r0 = X0_final(1:3);
v0 = X0_final(4:6);
k_SRP_final = X0_final(7);
bias_final  = X0_final(8);
lat4_final  = X0_final(9);
lon4_final  = X0_final(10);

fprintf('\n=== FINAL 10-STATE CKF RESULT (t0) ===\n');
fprintf('r0 (km):      [%.6f  %.6f  %.6f]\n', r0);
fprintf('v0 (km/s):    [%.6f  %.6f  %.6f]\n', v0);
fprintf('k_SRP:        %.6f\n', k_SRP_final);
fprintf('rho_dot_bias: %.6e km/s (%.3f mm/s)\n', bias_final, bias_final*1e6);
fprintf('lat_4 (deg):  %.6f\n', lat4_final * C.rad2deg);
fprintf('lon_4 (deg):  %.6f\n', lon4_final * C.rad2deg);

%% 7. Save
X0_final_kalman = X0_final;
P0_final_kalman = P0_final;

save('ASTE583_FinalKalman_Results.mat', ...
     'X0_final_kalman','P0_final_kalman', ...
     't_et','st_id','has_range', ...
     'rho_prefit','rho_postfit', ...
     'rr_prefit','rr_postfit');

fprintf('Saved to ASTE583_FinalKalman_Results.mat\n');
toc
end
