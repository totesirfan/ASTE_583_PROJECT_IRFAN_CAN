function Phase2_Prelim_Plots()
% PHASE2_PRELIM_PLOTS
% Post-processing for PRELIM OD (0–6 days):
%   - Converged prefit/postfit range-rate residuals (Batch & Kalman)
%   - 3-sigma covariance ellipses in RTN in three separate planes:
%       (R,T), (R,N), (T,N)
%
% Requires:
%   - Phase2_Prelim_Batch*.m   -> ASTE583_PrelimBatch_Results.mat
%       with struct prelim_results.X0_batch, prelim_results.P0_batch
%   - Phase2_Prelim_Kalman*.m  -> ASTE583_PrelimKalman_*.mat
%       with struct kf_results.X0_KF, kf_results.P0_KF
%   - load_project_measurements.m
%   - lib_dynamics.m, lib_measurements.m, lib_constants.m

    clc; close all;
    fprintf('=== PHASE 2: PRELIM OD POST-PROCESSING (BATCH + KALMAN) ===\n');

    %% 1. Setup & constants
    try
        init_project();
    catch
    end
    const = lib_constants();
    t0_et = const.t_detect_et;

    %% 2. Load measurements and restrict to 0–6 days
    fprintf('Loading measurements...\n');
    meas_full = load_project_measurements();   % your helper

    t_sec  = meas_full.time_sec;              % seconds since detection
    st_id  = meas_full.station_id;
    rng_km = meas_full.range_km;
    rr_kmps = meas_full.rr_kmps;

    DCO_days_prelim = 6.0;
    use_idx = (t_sec <= DCO_days_prelim * const.day2sec);

    t_sec   = t_sec(use_idx);
    st_id   = st_id(use_idx);
    rng_km  = rng_km(use_idx);
    rr_kmps = rr_kmps(use_idx);

    [t_sec, sort_idx] = sort(t_sec);
    st_id   = st_id(sort_idx);
    rng_km  = rng_km(sort_idx);
    rr_kmps = rr_kmps(sort_idx);

    n_meas   = numel(t_sec);
    t_et_all = t0_et + t_sec;
    t_days   = t_sec / const.day2sec;

    fprintf('Using %d measurements up to %.3f days.\n', n_meas, max(t_days));

    %% 3. Load Prelim Batch & Kalman results
    fprintf('Loading prelim Batch & Kalman MAT files...\n');

    % --- Batch ---
    Sbatch = load('ASTE583_PrelimBatch_Results.mat');
    prelim_results = Sbatch.prelim_results;
    X0_batch = prelim_results.X0_batch;   % 10x1
    P0_batch = prelim_results.P0_batch;   % 10x10

    % --- Kalman ---
    if exist('ASTE583_PrelimKalman_Results.mat','file')
        Skal = load('ASTE583_PrelimKalman_Results.mat');
    else
        Skal = load('ASTE583_PrelimKalman_Debug.mat');
    end
    kf_results = Skal.kf_results;
    X0_kalman  = kf_results.X0_KF;        % 10x1
    P0_kalman  = kf_results.P0_KF;        % 10x10

    % Reference initial state (for PREFIT only)
    stat4 = const.stations(4);
    X0_ref = [ const.X0_ref;
               const.k_SRP_0;
               0.0;
               stat4.lat;
               stat4.lon ];

    %% 4. Helper: propagate 10-state from t0 to all measurement times
    function X_all = propagate_10state(X0_init)
        % X0_init: 10x1 at t0
        % X_all:   10 x n_meas at t_et_all (same order)

        [t_unique, ~, idx_back] = unique(t_et_all);
        tspan = [t0_et; t_unique(:)];

        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [T_int, X_int] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                               tspan, X0_init, opts);

        X_unique = zeros(10, numel(t_unique));
        for j = 1:numel(t_unique)
            tj  = t_unique(j);
            idx = find(T_int <= tj, 1, 'last');
            if T_int(idx) == tj || idx == numel(T_int)
                X_unique(:,j) = X_int(idx,:).';
            else
                alpha = (tj - T_int(idx)) / (T_int(idx+1) - T_int(idx));
                X_unique(:,j) = (1-alpha)*X_int(idx,:).' + alpha*X_int(idx+1,:).';
            end
        end

        X_all = X_unique(:, idx_back);
    end

    %% 5. Propagate reference, Batch and Kalman trajectories
    fprintf('Propagating reference, Batch, and Kalman trajectories...\n');
    X_ref_all   = propagate_10state(X0_ref);
    X_batch_all = propagate_10state(X0_batch);
    X_kal_all   = propagate_10state(X0_kalman);

    %% 6. Compute prefit & postfit Doppler residuals
    fprintf('Computing computed measurements and residuals...\n');

    rr_obs = rr_kmps;

    rr_prefit_batch  = zeros(n_meas,1);
    rr_postfit_batch = zeros(n_meas,1);
    rr_prefit_kal    = zeros(n_meas,1);
    rr_postfit_kal   = zeros(n_meas,1);

    for k = 1:n_meas
        t_k = t_et_all(k);
        st  = st_id(k);

        % Reference for PREFIT
        [Y_ref, ~] = lib_measurements(t_k, X_ref_all(:,k), st, const);
        rr_ref = Y_ref(2);

        % Batch postfit
        [Y_b, ~] = lib_measurements(t_k, X_batch_all(:,k), st, const);
        rr_b = Y_b(2);

        % Kalman postfit
        [Y_k, ~] = lib_measurements(t_k, X_kal_all(:,k), st, const);
        rr_kal = Y_k(2);

        rr_prefit_batch(k)  = rr_obs(k) - rr_ref;
        rr_postfit_batch(k) = rr_obs(k) - rr_b;

        rr_prefit_kal(k)    = rr_obs(k) - rr_ref;
        rr_postfit_kal(k)   = rr_obs(k) - rr_kal;
    end

    %% 7. Plot prefit / postfit RR residuals – Batch
    fprintf('Plotting Batch prefit/postfit RR residuals...\n');

    figure('Name','Prelim Batch: Range-Rate Residuals', ...
           'Color','w','Position',[50 100 1000 400]);

    subplot(1,2,1); hold on; grid on;
    scatter(t_days, 1e3*rr_prefit_batch, 4, 'filled');
    xlabel('Time since detection (days)');
    ylabel('Prefit \rhȯ residual (mm/s)');
    title('Prelim Batch - Prefit Range-Rate Residuals');
    xlim([0 max(t_days)]);

    subplot(1,2,2); hold on; grid on;
    scatter(t_days, 1e3*rr_postfit_batch, 4, 'filled');
    xlabel('Time since detection (days)');
    ylabel('Postfit \rhȯ residual (mm/s)');
    title('Prelim Batch - Postfit Range-Rate Residuals');
    xlim([0 max(t_days)]);

    %% 8. Plot prefit / postfit RR residuals – Kalman
    fprintf('Plotting Kalman prefit/postfit RR residuals...\n');

    figure('Name','Prelim Kalman: Range-Rate Residuals', ...
           'Color','w','Position',[50 550 1000 400]);

    subplot(1,2,1); hold on; grid on;
    scatter(t_days, 1e3*rr_prefit_kal, 4, 'filled');
    xlabel('Time since detection (days)');
    ylabel('Prefit \rhȯ residual (mm/s)');
    title('Prelim Kalman - Prefit Range-Rate Residuals');
    xlim([0 max(t_days)]);

    subplot(1,2,2); hold on; grid on;
    scatter(t_days, 1e3*rr_postfit_kal, 4, 'filled');
    xlabel('Time since detection (days)');
    ylabel('Postfit \rhȯ residual (mm/s)');
    title('Prelim Kalman - Postfit Range-Rate Residuals');
    xlim([0 max(t_days)]);

    %% 9. 3-σ covariance in RTN: three separate planes (R-T, R-N, T-N), no reference
    fprintf('Plotting 3σ covariance ellipses in RTN planes (Batch vs Kalman)...\n');

    figure('Name','Prelim OD: 3\sigma Covariance Ellipses in RTN', ...
           'Color','w','Position',[1100 200 1400 450]);

    % --- Helper: RTN transform from inertial r, v ---
    function T_rtn = compute_T_RTN(r_vec, v_vec)
        r_hat = r_vec / norm(r_vec);
        h_vec = cross(r_vec, v_vec);
        h_hat = h_vec / norm(h_vec);
        t_hat = cross(h_hat, r_hat);
        T_rtn = [r_hat.'; t_hat.'; h_hat.'];   % rows: R, T, N
    end

    % --- Helper: 2D 3σ covariance ellipse ---
    function plot_cov_ellipse(P2, x0, y0, style)
        [V,D] = eig(P2);
        s = 3 * sqrt(diag(D));              % 3-sigma eigenlengths
        th = linspace(0, 2*pi, 200);
        circ = [cos(th); sin(th)];
        ell  = V * (diag(s) * circ);
        plot(x0 + ell(1,:), y0 + ell(2,:), style, 'LineWidth', 1.5);
    end

    % Position & covariance at detection
    r0_batch = X0_batch(1:3);
    v0_batch = X0_batch(4:6);
    r0_kal   = X0_kalman(1:3);

    Pxyz_batch = P0_batch(1:3,1:3);
    Pxyz_kal   = P0_kalman(1:3,1:3);

    % RTN frame from Batch solution
    T_rtn = compute_T_RTN(r0_batch, v0_batch);

    % Transform covariances
    P_RTN_batch = T_rtn * Pxyz_batch * T_rtn.';
    P_RTN_kal   = T_rtn * Pxyz_kal   * T_rtn.';

    % Mean position difference (Kalman relative to Batch) in RTN
    dr_inertial = r0_kal - r0_batch;
    dRTN        = T_rtn * dr_inertial;
    dR = dRTN(1);
    dT = dRTN(2);
    dN = dRTN(3);

    % 2x2 sub-blocks
    P_RT_batch = P_RTN_batch(1:2,1:2);
    P_RT_kal   = P_RTN_kal(1:2,1:2);

    P_RN_batch = P_RTN_batch([1 3],[1 3]);
    P_RN_kal   = P_RTN_kal([1 3],[1 3]);

    P_TN_batch = P_RTN_batch(2:3,2:3);
    P_TN_kal   = P_RTN_kal(2:3,2:3);

    % --- Subplot 1: R vs T ---
    subplot(1,3,1); hold on; grid on; axis equal;
    plot(0, 0, 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);        % Batch
    plot(dR, dT, 'rs', 'MarkerFaceColor','r', 'MarkerSize',6);      % Kalman
    plot_cov_ellipse(P_RT_batch, 0, 0, 'b-');
    plot_cov_ellipse(P_RT_kal,   dR, dT, 'r-');
    xlabel('\DeltaR (km)');
    ylabel('\DeltaT (km)');
    title('\DeltaR-\DeltaT (RT plane)');
    legend('Batch','Kalman','Batch 3\sigma','Kalman 3\sigma', ...
           'Location','best');

    % --- Subplot 2: R vs N ---
    subplot(1,3,2); hold on; grid on; axis equal;
    plot(0, 0, 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
    plot(dR, dN, 'rs', 'MarkerFaceColor','r', 'MarkerSize',6);
    plot_cov_ellipse(P_RN_batch, 0, 0, 'b-');
    plot_cov_ellipse(P_RN_kal,   dR, dN, 'r-');
    xlabel('\DeltaR (km)');
    ylabel('\DeltaN (km)');
    title('\DeltaR-\DeltaN (RN plane)');

    % --- Subplot 3: T vs N ---
    subplot(1,3,3); hold on; grid on; axis equal;
    plot(0, 0, 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
    plot(dT, dN, 'rs', 'MarkerFaceColor','r', 'MarkerSize',6);
    plot_cov_ellipse(P_TN_batch, 0, 0, 'b-');
    plot_cov_ellipse(P_TN_kal,   dT, dN, 'r-');
    xlabel('\DeltaT (km)');
    ylabel('\DeltaN (km)');
    title('\DeltaT-\DeltaN (TN plane)');

    fprintf('Done. Prefit/postfit + RTN 3σ plots are ready for the Prelim OD section.\n');
end
