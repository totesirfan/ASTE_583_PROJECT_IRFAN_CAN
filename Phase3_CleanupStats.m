function Phase3_CleanupStats()
% PHASE2_CLEANUP_STATISTICS
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

    X0_KF = Sfinal.X0_final_kalman(:);   % 10x1
    P0_KF = Sfinal.P0_final_kalman;      % 10x10

    if numel(X0_KF) ~= 10 || any(size(P0_KF) ~= 10)
        error('FINAL Kalman state must be 10x1 and covariance 10x10.');
    end

    %% 3. Reference trajectory (Prime Nav reference + deterministic LTM)
    stat4 = const.stations(4);
    X_ref0 = [const.X0_ref;
              const.k_SRP_0;
              0.0;              % bias
              stat4.lat;
              stat4.lon];       % Antarctica "truth"

    % Detection -> LTM
    [~, Xref_seg1] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t0_et t_LTM], X_ref0, opts);
    X_ref_LTMm = Xref_seg1(end,:).';

    % Apply deterministic LTM
    dV_LTM_nom = const.LTM.dV(:);        % km/s
    X_ref_LTMp = X_ref_LTMm;
    X_ref_LTMp(4:6) = X_ref_LTMp(4:6) + dV_LTM_nom;

    % LTM -> LCM
    [~, Xref_seg2] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LTM t_LCM], X_ref_LTMp, opts);
    X_ref_LCM = Xref_seg2(end,:).';

    % LCM -> target (no cleanup)
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

    % Storage
    pos_err_LCM = zeros(Nmc,1);   % km
    vel_err_LCM = zeros(Nmc,1);   % km/s
    dV_LTM_ms   = zeros(Nmc,1);   % m/s
    dV_LCM_ms   = zeros(Nmc,1);   % m/s
    dV_tot_ms   = zeros(Nmc,1);   % m/s

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
        dir_vec = randn(3,1); dir_vec = dir_vec / norm(dir_vec);
        mag_err = sigma_exec * randn();          % km/s
        dV_err  = mag_err * dir_vec;

        dV_LTM_actual = dV_LTM_nom + dV_err;     % km/s

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

        % LCM -> target without cleanup
        [~, Xseg3] = ode45(@(t,x) lib_dynamics(t,x,const), ...
                           [t_LCM t_TGT], X_LCM_samp, opts);
        X_TGT_free = Xseg3(end,:).';
        r_TGT_free = X_TGT_free(1:3);

        % Feedback law: choose DeltaV_LCM to kill position error at TGT
        dr_TGT     = r_TGT_free - r_ref_TGT;   % 3x1
        dV_LCM_vec = -S \ dr_TGT;              % km/s

        dV_LCM_ms(k) = 1000*norm(dV_LCM_vec);      % m/s
        dV_LTM_ms(k) = 1000*norm(dV_LTM_actual);   % m/s
        dV_tot_ms(k) = dV_LTM_ms(k) + dV_LCM_ms(k);
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
    xlabel('Total \DeltaV (LTM + LCM) (m/s)', 'Interpreter','tex');
    ylabel('Count');
    title('Total \DeltaV distribution','Interpreter','tex');
    legend('Samples','Budget','Location','best');

    fprintf('\nCleanup maneuver statistics plotting complete.\n');
    toc
end
