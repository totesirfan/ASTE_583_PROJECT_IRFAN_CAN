function Phase3_CleanupStats()
% PHASE3_CLEANUP_STATISTICS
% Cleanup Maneuver Statistics (20 pts) using FINAL Kalman OD.
%
% Implements handout requirements:
%  - Sample a-posteriori covariance from FINAL OD (Kalman).
%  - For each sample, propagate from detection to LCM, including LTM,
%    with maneuver execution error applied to LTM.
%  - At LCM, record |r - r_ref| and |v - v_ref|.
%  - For each sample, use a linear feedback law at LCM to target the
%    reference trajectory 30 days after LCM and compute required LCM dV.
%  - Plot histograms of position/velocity dispersion at LCM, LCM dV,
%    and total dV = |LTM| + |LCM|.
%  - Compute DeltaV_99 (99th percentile of total dV) and compare to
%    total dV requirement (const.dV_budget, in km/s).
%
% Additional diagnostics (added here):
%  - Lunar flyby closest-approach sanity check using the FINAL OD
%    reference and a subset of Monte Carlo trajectories with cleanup.

    clear; clc; close all;
    fprintf('=== Cleanup maneuver statistics using FINAL Kalman OD ===\n');
    tic

    %% 1. Setup, constants, epochs
    try
        init_project();
    catch
    end
    const  = lib_constants();
    t0_et  = const.t_detect_et;
    day    = const.day2sec;

    % Maneuver epochs
    t_LTM  = cspice_str2et(const.LTM.date_utc);
    t_LCM  = cspice_str2et(const.LCM.date_utc);
    t_TGT  = t_LCM + 30*day;      % 30 days after LCM

    opts   = odeset('RelTol',1e-12,'AbsTol',1e-12);

    %% 2. Load FINAL Kalman OD solution (matches save in Phase2_Final_Kalman)
    if ~exist('ASTE583_FinalKalman_Results.mat','file')
        error('ASTE583_FinalKalman_Results.mat not found. Run Phase2_Final_Kalman first.');
    end
    Sfinal = load('ASTE583_FinalKalman_Results.mat');

    if ~isfield(Sfinal,'X0_final_kalman') || ~isfield(Sfinal,'P0_final_kalman')
        error('ASTE583_FinalKalman_Results.mat must contain X0_final_kalman and P0_final_kalman.');
    end

    X0_KF = Sfinal.X0_final_kalman(:);   % 10x1, FINAL OD mean at detection
    P0_KF = Sfinal.P0_final_kalman;      % 10x10, FINAL OD covariance at detection

    if numel(X0_KF) ~= 10 || any(size(P0_KF) ~= 10)
        error('FINAL Kalman state must be 10x1 and covariance 10x10.');
    end

    %% 3. Reference trajectory (FINAL OD mean + deterministic LTM)
    % IMPORTANT CHANGE: use X0_KF as the reference initial state so that
    % the reference trajectory uses the same k_SRP, Doppler bias, and
    % Antarctica coordinates as the navigation solution.
    X_ref0 = X0_KF;     % 10x1: [r; v; k_SRP; bias; phi4; lambda4]

    % Detection -> LTM (reference)
    [~, Xref_seg1] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t0_et t_LTM], X_ref0, opts);
    X_ref_LTMm = Xref_seg1(end,:).';

    % Apply deterministic LTM (impulsive)
    dV_LTM_nom = const.LTM.dV(:);        % km/s
    X_ref_LTMp = X_ref_LTMm;
    X_ref_LTMp(4:6) = X_ref_LTMp(4:6) + dV_LTM_nom;

    % LTM -> LCM (reference)
    [~, Xref_seg2] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LTM t_LCM], X_ref_LTMp, opts);
    X_ref_LCM = Xref_seg2(end,:).';

    % LCM -> target (no cleanup, reference)
    [~, Xref_seg3] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LCM t_TGT], X_ref_LCM, opts);
    X_ref_TGT = Xref_seg3(end,:).';

    r_ref_LCM = X_ref_LCM(1:3);
    v_ref_LCM = X_ref_LCM(4:6);
    r_ref_TGT = X_ref_TGT(1:3);

    %% 4. Feedback guidance law sensitivity S = dr_TGT / d(DeltaV_LCM)
    dv_step = 1e-3;   % 1 m/s in km/s
    S = zeros(3,3);
    for j = 1:3
        dv = zeros(3,1);
        dv(j) = dv_step;

        X_LCM_pert = X_ref_LCM;
        X_LCM_pert(4:6) = X_LCM_pert(4:6) + dv;

        [~, Xpert_seg] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                               [t_LCM t_TGT], X_LCM_pert, opts);
        r_pert = Xpert_seg(end,1:3).';
        S(:,j) = (r_pert - r_ref_TGT) / dv_step;
    end

    %% 5. Monte Carlo setup
    Nmc = 1000;
    rng(583);  % reproducible

    % Cholesky factor of posterior covariance
    P_sym = 0.5*(P0_KF + P0_KF.');
    [Lchol,p] = chol(P_sym,'lower');
    if p > 0
        jitter = 1e-10*max(1, max(diag(P_sym)));
        [Lchol,p] = chol(P_sym + jitter*eye(10),'lower');
        if p > 0
            error('Final covariance is not positive definite, even after jitter.');
        end
    end

    % Storage for statistics
    pos_err_LCM = zeros(Nmc,1);   % km
    vel_err_LCM = zeros(Nmc,1);   % km/s
    dV_LTM_ms   = zeros(Nmc,1);   % m/s
    dV_LCM_ms   = zeros(Nmc,1);   % m/s
    dV_tot_ms   = zeros(Nmc,1);   % m/s

    % Storage for later lunar-flyby check (store all, we can subsample)
    LCM_state_store = zeros(10, Nmc);   % state at LCM (before cleanup)
    dV_LCM_store    = zeros(3,  Nmc);   % cleanup DeltaV vector (km/s)

    fprintf('Running %d Monte Carlo samples...\n\n', Nmc);

    sigma_exec = const.LTM.sigma_exec;   % 5e-3 km/s (5 m/s, spherical)

    %% 6. Monte Carlo loop
    for k = 1:Nmc
        % 6.1 Sample initial state at detection from final posterior
        xi      = randn(10,1);
        X0_samp = X0_KF + Lchol*xi;

        % Detection -> LTM
        [~, Xseg1] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t0_et t_LTM], X0_samp, opts);
        X_LTMm = Xseg1(end,:).';

        % LTM execution error (spherical 1-sigma)
        dir_vec = randn(3,1); 
        dir_vec = dir_vec / norm(dir_vec);
        mag_err = sigma_exec * randn();   % km/s
        dV_err  = mag_err * dir_vec;      % km/s

        dV_LTM_actual = dV_LTM_nom + dV_err;    % km/s

        X_LTMp = X_LTMm;
        X_LTMp(4:6) = X_LTMp(4:6) + dV_LTM_actual;

        % LTM -> LCM
        [~, Xseg2] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LTM t_LCM], X_LTMp, opts);
        X_LCM_samp = Xseg2(end,:).';

        r_LCM = X_LCM_samp(1:3);
        v_LCM = X_LCM_samp(4:6);

        % Dispersion at LCM
        pos_err_LCM(k) = norm(r_LCM - r_ref_LCM);   % km
        vel_err_LCM(k) = norm(v_LCM - v_ref_LCM);   % km/s

        % LCM -> target without cleanup (free drift)
        [~, Xseg3] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LCM t_TGT], X_LCM_samp, opts);
        X_TGT_free = Xseg3(end,:).';
        r_TGT_free = X_TGT_free(1:3);

        % Feedback law: choose DeltaV_LCM to kill position error at TGT
        dr_TGT     = r_TGT_free - r_ref_TGT;   % 3x1
        dV_LCM_vec = -S \ dr_TGT;              % km/s

        % Store magnitudes for histograms
        dV_LCM_ms(k) = 1000*norm(dV_LCM_vec);      % m/s
        dV_LTM_ms(k) = 1000*norm(dV_LTM_actual);   % m/s
        dV_tot_ms(k) = dV_LTM_ms(k) + dV_LCM_ms(k);

        % Store full vectors for later lunar-flyby check
        LCM_state_store(:,k) = X_LCM_samp;
        dV_LCM_store(:,k)    = dV_LCM_vec;
    end

    %% 7. Scalar statistics and requirement check
    mean_dV_LCM = mean(dV_LCM_ms);
    std_dV_LCM  = std(dV_LCM_ms);
    min_dV_LCM  = min(dV_LCM_ms);
    max_dV_LCM  = max(dV_LCM_ms);

    dV_tot_sorted = sort(dV_tot_ms);
    idx99  = max(1, ceil(0.99*Nmc));
    dV99   = dV_tot_sorted(idx99);           % m/s

    dV_budget_ms = const.dV_budget * 1000;   % requirement in m/s
    frac_within  = 100*mean(dV_tot_ms <= dV_budget_ms);

    fprintf('--- Cleanup maneuver statistics (LCM) ---\n');
    fprintf('Samples              : %d\n', Nmc);
    fprintf('Mean LCM dV          : %.2f m/s\n', mean_dV_LCM);
    fprintf('Std dev LCM dV       : %.2f m/s\n', std_dV_LCM);
    fprintf('Min / Max LCM dV     : %.2f / %.2f m/s\n', ...
            min_dV_LCM, max_dV_LCM);
    fprintf('dV_budget (LTM+LCM)  : %.1f m/s\n', dV_budget_ms);
    fprintf('DeltaV_99 (total)    : %.2f m/s\n', dV99);
    fprintf('Frac(total <= budget): %.1f %%\n', frac_within);

    %% 8. Plots: dispersion at LCM
    figure('Name','LCM dispersion (position & velocity)', ...
           'Color','w','Position',[100 100 1000 400]);

    subplot(1,2,1); hold on; grid on;
    histogram(pos_err_LCM, 40);
    xlabel('|r - r_{ref}| at LCM (km)','Interpreter','tex');
    ylabel('Count');
    title('Position magnitude error at LCM');

    subplot(1,2,2); hold on; grid on;
    histogram(1000*vel_err_LCM, 40);
    xlabel('|v - v_{ref}| at LCM (m/s)','Interpreter','tex');
    ylabel('Count');
    title('Velocity magnitude error at LCM');

    %% 9. Plots: LCM dV and total dV
    figure('Name','Cleanup maneuver DeltaV statistics', ...
           'Color','w','Position',[150 550 1000 400]);

    subplot(1,2,1); hold on; grid on;
    histogram(dV_LCM_ms, 40);
    xlabel('LCM \DeltaV (m/s)','Interpreter','tex');
    ylabel('Count');
    title('LCM \DeltaV magnitude','Interpreter','tex');

    subplot(1,2,2); hold on; grid on;
    histogram(dV_tot_ms, 40);
    xline(dV_budget_ms,'r--','Budget','LineWidth',1.5);
    xline(dV99,'k--','\DeltaV_{99}','LineWidth',1.5);
    xlabel('Total \DeltaV (LTM + LCM) (m/s)', 'Interpreter','tex');
    ylabel('Count');
    title('Total \DeltaV distribution','Interpreter','tex');
    legend('Samples','Budget','\DeltaV_{99}','Location','best');

    %% 10. Lunar flyby closest-approach sanity check
    % NOTE: This is *not* part of the handout requirements, but is a useful
    % diagnostic to verify that the post-LCM + cleanup trajectories still
    % produce a reasonable lunar flyby in the far future.
    fprintf('\nPerforming lunar flyby closest-approach sanity check...\n');

    % Choose a horizon long enough to cover the flyby (adjust as needed).
    % Example: 200 days after LCM.
    t_fly_end = t_LCM + 200*day;

    % 10.1 Reference trajectory: from LTM (after nominal LTM) to t_fly_end
    [T_ref_fly, X_ref_fly] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                                   [t_LTM t_fly_end], X_ref_LTMp, opts);

    n_ref = numel(T_ref_fly);
    dist_ref = zeros(n_ref,1);
    for i = 1:n_ref
        % Moon state in J2000, then rotate to EMO2000
        [stateM, ~] = cspice_spkezr('MOON', T_ref_fly(i), ...
                                    'J2000', 'NONE', 'SUN');
        rM_EMO = const.R_EME_EMO * stateM(1:3);  % km, EMO2000
        r_sc   = X_ref_fly(i,1:3).';
        dist_ref(i) = norm(r_sc - rM_EMO);
    end
    [ca_ref, idx_ref] = min(dist_ref);
    t_ca_ref = T_ref_fly(idx_ref);
    utc_ca_ref = cspice_et2utc(t_ca_ref,'ISOC',3);

    fprintf('Reference closest approach to Moon (LTM->+200 d): %.1f km at %s UTC\n', ...
            ca_ref, utc_ca_ref);

    % 10.2 Monte Carlo subset with cleanup applied at LCM
    max_check = min(50, Nmc);  % check first 50 samples to keep runtime reasonable
    ca_mc = zeros(max_check,1);

    for k = 1:max_check
        % State just before cleanup
        X_LCM_pre = LCM_state_store(:,k);

        % Apply cleanup DeltaV at LCM
        X_LCM_post        = X_LCM_pre;
        X_LCM_post(4:6)   = X_LCM_post(4:6) + dV_LCM_store(:,k);

        % Propagate from LCM to t_fly_end with cleanup
        [T_mc, X_mc] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                             [t_LCM t_fly_end], X_LCM_post, opts);

        n_mc = numel(T_mc);
        dist_mc_k = zeros(n_mc,1);
        for i = 1:n_mc
            [stateM, ~] = cspice_spkezr('MOON', T_mc(i), ...
                                        'J2000', 'NONE', 'SUN');
            rM_EMO = const.R_EME_EMO * stateM(1:3);
            r_sc   = X_mc(i,1:3).';
            dist_mc_k(i) = norm(r_sc - rM_EMO);
        end
        ca_mc(k) = min(dist_mc_k);
    end

    fprintf('MC mean closest approach (first %d samples): %.1f km\n', ...
            max_check, mean(ca_mc));
    fprintf('MC std  closest approach                 : %.1f km\n', ...
            std(ca_mc));

    fprintf('\nCleanup maneuver statistics and lunar-flyby check complete.\n');
    toc
end
