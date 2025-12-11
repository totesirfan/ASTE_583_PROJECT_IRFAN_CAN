function DEPRECATED_Phase3_Cleanup_Statistics_LCMExec()
% DEPERECATED_Phase2_Final_EKF
% Cleanup Maneuver Statistics with LCM execution error + post-LCM propagation.
%
% Based on Phase2_Cleanup_Statistics, but:
%   - LTM has 5 m/s spherical execution error (const.LTM.sigma_exec).
%   - LCM also has 5 m/s spherical execution error (assumed same sigma).
%   - At LCM:
%       * record |r - r_ref| and |v - v_ref|.
%   - At TGT = 30 days after LCM:
%       * free-run (no cleanup): r_TGT_free
%       * commanded cleanup (perfect execution): r_TGT_cmd
%       * executed cleanup (with LCM error): r_TGT_exec
%   - Use linear feedback law at LCM (via sensitivity S = dr_TGT/d(ΔV_LCM))
%     to compute ΔV_LCM_cmd.
%   - Histograms:
%       * position/velocity dispersion at LCM
%       * position dispersion at TGT (cmd vs exec)
%       * LCM ΔV (commanded) and total ΔV = |LTM_actual| + |LCM_cmd|
%   - Compute ΔV_99 (99th percentile of total commanded ΔV) vs budget.

    clear; clc; close all;
    fprintf('=== Cleanup maneuver statistics (LTM + LCM exec error) ===\n');

    %% 1. Setup / constants / epochs
    try, init_project(); end
    const  = lib_constants();
    t0_et  = const.t_detect_et;
    day    = const.day2sec;

    t_LTM  = cspice_str2et(const.LTM.date_utc);
    t_LCM  = cspice_str2et(const.LCM.date_utc);
    t_TGT  = t_LCM + 30*day;   % 30 days after LCM

    opts   = odeset('RelTol',1e-12,'AbsTol',1e-12);

    %% 2. Load FINAL Kalman OD state/covariance
    fn = 'ASTE583_FinalKalman_Results.mat';
    if ~exist(fn,'file')
        error('%s not found. Run Phase2_Final_Kalman first.', fn);
    end
    S = load(fn);
    if ~isfield(S,'X0_final_kalman') || ~isfield(S,'P0_final_kalman')
        error('%s must contain X0_final_kalman and P0_final_kalman.', fn);
    end

    X0_KF = S.X0_final_kalman(:);   % 10×1
    P0_KF = S.P0_final_kalman;      % 10×10

    %% 3. Reference trajectory (Prime Nav + deterministic LTM)
    stat4  = const.stations(4);
    X_ref0 = [ const.X0_ref;
               const.k_SRP_0;
               0.0;                % bias
               stat4.lat;
               stat4.lon ];

    % Detect → LTM
    [~, Xref1] = ode45(@(t,x) lib_dynamics(t,x,const), [t0_et t_LTM], X_ref0, opts);
    X_ref_LTMm = Xref1(end,:).';

    % Apply nominal LTM
    dV_LTM_nom   = const.LTM.dV(:);   % km/s
    X_ref_LTMp   = X_ref_LTMm;
    X_ref_LTMp(4:6) = X_ref_LTMp(4:6) + dV_LTM_nom;

    % LTM → LCM
    [~, Xref2] = ode45(@(t,x) lib_dynamics(t,x,const), [t_LTM t_LCM], X_ref_LTMp, opts);
    X_ref_LCM  = Xref2(end,:).';

    % LCM → TGT (no cleanup)
    [~, Xref3] = ode45(@(t,x) lib_dynamics(t,x,const), [t_LCM t_TGT], X_ref_LCM, opts);
    X_ref_TGT  = Xref3(end,:).';

    r_ref_LCM = X_ref_LCM(1:3);
    v_ref_LCM = X_ref_LCM(4:6);
    r_ref_TGT = X_ref_TGT(1:3);

    %% 4. Sensitivity S = dr_TGT / d(ΔV_LCM)
    dv_step = 1e-3;   % 1 m/s in km/s
    Ssens   = zeros(3,3);
    for j = 1:3
        dv = zeros(3,1); dv(j) = dv_step;
        X_LCM_pert       = X_ref_LCM;
        X_LCM_pert(4:6)  = X_LCM_pert(4:6) + dv;
        [~, Xpert]       = ode45(@(t,x) lib_dynamics(t,x,const), [t_LCM t_TGT], X_LCM_pert, opts);
        r_pert           = Xpert(end,1:3).';
        Ssens(:,j)       = (r_pert - r_ref_TGT)/dv_step;
    end

    %% 5. Monte Carlo setup
    Nmc = 1000;
    rng(583); % reproducible

    % Make covariance symmetric and get Cholesky
    P = 0.5*(P0_KF + P0_KF.');
    [L,p] = chol(P,'lower');
    if p>0
        jit = 1e-10*max(1,max(diag(P)));
        [L,p] = chol(P + jit*eye(10),'lower');
        if p>0, error('Final covariance not PD, even after jitter.'); end
    end

    % Exec sigmas (km/s)
    sigma_LTM = const.LTM.sigma_exec;          % 5e-3 km/s
    sigma_LCM = sigma_LTM;                     % assume same for LCM; change if needed

    % Storage
    pos_err_LCM      = zeros(Nmc,1);  % |r - r_ref| at LCM (km)
    vel_err_LCM      = zeros(Nmc,1);  % |v - v_ref| at LCM (km/s)
    pos_err_TGT_cmd  = zeros(Nmc,1);  % |r - r_ref| at TGT (perfect LCM)
    pos_err_TGT_exec = zeros(Nmc,1);  % |r - r_ref| at TGT (LCM exec error)

    dV_LTM_ms        = zeros(Nmc,1);  % actual LTM ΔV magnitude (m/s)
    dV_LCM_cmd_ms    = zeros(Nmc,1);  % commanded cleanup ΔV magnitude (m/s)
    dV_LCM_exec_ms   = zeros(Nmc,1);  % executed cleanup ΔV magnitude (m/s)
    dV_tot_ms        = zeros(Nmc,1);  % total commanded ΔV (LTM actual + LCM cmd) (m/s)

    fprintf('Running %d Monte Carlo samples...\n', Nmc);

    %% 6. Monte Carlo loop
    for k = 1:Nmc
        % 6.1 Sample initial state at detection
        xi      = randn(10,1);
        X0_samp = X0_KF + L*xi;

        % Detect → LTM (no burn)
        [~, X1]   = ode45(@(t,x) lib_dynamics(t,x,const), [t0_et t_LTM], X0_samp, opts);
        X_LTMm    = X1(end,:).';

        % 6.2 LTM execution error (spherical)
        dirL = randn(3,1); dirL = dirL/norm(dirL);
        magL = sigma_LTM * randn();         % km/s
        dV_LTM_err  = magL * dirL;
        dV_LTM_act  = dV_LTM_nom + dV_LTM_err;

        X_LTMp = X_LTMm;
        X_LTMp(4:6) = X_LTMp(4:6) + dV_LTM_act;

        % LTM → LCM
        [~, X2]    = ode45(@(t,x) lib_dynamics(t,x,const), [t_LTM t_LCM], X_LTMp, opts);
        X_LCM_samp = X2(end,:).';

        r_LCM = X_LCM_samp(1:3);
        v_LCM = X_LCM_samp(4:6);

        % Dispersion at LCM
        pos_err_LCM(k) = norm(r_LCM - r_ref_LCM);
        vel_err_LCM(k) = norm(v_LCM - v_ref_LCM);

        % 6.3 Free-run to TGT with NO cleanup
        [~, X3]       = ode45(@(t,x) lib_dynamics(t,x,const), [t_LCM t_TGT], X_LCM_samp, opts);
        X_TGT_free    = X3(end,:).';
        r_TGT_free    = X_TGT_free(1:3);
        dr_TGT        = r_TGT_free - r_ref_TGT;

        % Feedback law: ΔV_LCM_cmd to drive r_TGT → r_ref_TGT (linearized)
        dV_LCM_cmd    = -Ssens \ dr_TGT;    % km/s
        dV_LCM_cmd_ms(k)  = 1000*norm(dV_LCM_cmd);

        % 6.4 Perfect LCM execution case (for comparison)
        X_LCM_cmd     = X_LCM_samp;
        X_LCM_cmd(4:6)= X_LCM_cmd(4:6) + dV_LCM_cmd;
        [~, X4]       = ode45(@(t,x) lib_dynamics(t,x,const), [t_LCM t_TGT], X_LCM_cmd, opts);
        X_TGT_cmd     = X4(end,:).';
        pos_err_TGT_cmd(k) = norm(X_TGT_cmd(1:3) - r_ref_TGT);

        % 6.5 LCM execution error (spherical)
        dirC = randn(3,1); dirC = dirC/norm(dirC);
        magC = sigma_LCM * randn();         % km/s
        dV_LCM_err  = magC * dirC;
        dV_LCM_exec = dV_LCM_cmd + dV_LCM_err;

        dV_LCM_exec_ms(k) = 1000*norm(dV_LCM_exec);

        % Apply executed LCM and propagate again to TGT
        X_LCM_exec      = X_LCM_samp;
        X_LCM_exec(4:6) = X_LCM_exec(4:6) + dV_LCM_exec;
        [~, X5]         = ode45(@(t,x) lib_dynamics(t,x,const), [t_LCM t_TGT], X_LCM_exec, opts);
        X_TGT_exec      = X5(end,:).';
        pos_err_TGT_exec(k) = norm(X_TGT_exec(1:3) - r_ref_TGT);

        % 6.6 ΔV bookkeeping (budget uses commanded LCM, actual LTM)
        dV_LTM_ms(k) = 1000*norm(dV_LTM_act);
        dV_tot_ms(k) = dV_LTM_ms(k) + dV_LCM_cmd_ms(k);
    end

    %% 7. Scalar stats and requirement check
    dV_budget_ms = const.dV_budget * 1000;   % km/s → m/s

    dV_tot_sorted = sort(dV_tot_ms);
    idx99         = max(1, ceil(0.99*Nmc));
    dV99          = dV_tot_sorted(idx99);

    fprintf('\n--- Cleanup maneuver statistics (with LCM exec error) ---\n');
    fprintf('Samples                 : %d\n', Nmc);
    fprintf('Mean LCM dV (cmd)       : %.2f m/s\n', mean(dV_LCM_cmd_ms));
    fprintf('Std  LCM dV (cmd)       : %.2f m/s\n', std(dV_LCM_cmd_ms));
    fprintf('Min / Max LCM dV (cmd)  : %.2f / %.2f m/s\n', ...
            min(dV_LCM_cmd_ms), max(dV_LCM_cmd_ms));
    fprintf('Mean LCM dV (exec)      : %.2f m/s\n', mean(dV_LCM_exec_ms));
    fprintf('Mean LTM dV (exec)      : %.2f m/s\n', mean(dV_LTM_ms));
    fprintf('dV_budget (LTM+LCM)     : %.1f m/s\n', dV_budget_ms);
    fprintf('DeltaV_99 (total, cmd)  : %.2f m/s\n', dV99);
    fprintf('Frac(total <= budget)   : %.1f %%\n', ...
            100*mean(dV_tot_ms <= dV_budget_ms));

    fprintf('\nPosition dispersion:\n');
    fprintf('  Mean |r - r_ref| @LCM     : %.2f km\n', mean(pos_err_LCM));
    fprintf('  Mean |r - r_ref| @TGT cmd : %.2f km\n', mean(pos_err_TGT_cmd));
    fprintf('  Mean |r - r_ref| @TGT exec: %.2f km\n', mean(pos_err_TGT_exec));
    fprintf('  99%% |r - r_ref| @TGT exec : %.2f km\n', prctile(pos_err_TGT_exec,99));

    %% 8. Plots – dispersion at LCM
    figure('Name','LCM dispersion', ...
           'Color','w','Position',[100 100 1000 400]);

    subplot(1,2,1); hold on; grid on;
    histogram(pos_err_LCM, 40);
    xlabel('|r - r_ref| at LCM (km)','Interpreter','none');
    ylabel('Count');
    title('Position error at LCM','Interpreter','none');

    subplot(1,2,2); hold on; grid on;
    histogram(1000*vel_err_LCM, 40);
    xlabel('|v - v_ref| at LCM (m/s)','Interpreter','none');
    ylabel('Count');
    title('Velocity error at LCM','Interpreter','none');

    %% 9. Plots – dispersion at TGT (30 days after LCM)
    figure('Name','TGT dispersion (30 days after LCM)', ...
           'Color','w','Position',[150 550 1000 400]);

    subplot(1,2,1); hold on; grid on;
    histogram(pos_err_TGT_cmd, 40);
    xlabel('|r - r_ref| at TGT (km) – perfect LCM','Interpreter','none');
    ylabel('Count');
    title('TGT error (no LCM exec error)','Interpreter','none');

    subplot(1,2,2); hold on; grid on;
    histogram(pos_err_TGT_exec, 40);
    xlabel('|r - r_ref| at TGT (km) – with LCM exec error','Interpreter','none');
    ylabel('Count');
    title('TGT error (with LCM exec error)','Interpreter','none');

    %% 10. Plots – ΔV statistics
    figure('Name','Cleanup ΔV statistics', ...
           'Color','w','Position',[200 1000 1000 400]);

    subplot(1,2,1); hold on; grid on;
    histogram(dV_LCM_cmd_ms, 40);
    xlabel('LCM \DeltaV commanded (m/s)','Interpreter','none');
    ylabel('Count');
    title('LCM commanded \DeltaV','Interpreter','none');

    subplot(1,2,2); hold on; grid on;
    histogram(dV_tot_ms, 40);
    xline(dV_budget_ms,'r--','Budget','LineWidth',1.5,'Interpreter','none');
    xlabel('Total \DeltaV (LTM exec + LCM cmd) (m/s)','Interpreter','none');
    ylabel('Count');
    title('Total commanded \DeltaV vs budget','Interpreter','none');
    legend({'Samples','Budget'},'Location','best');

    fprintf('\nCleanup maneuver statistics (with LCM exec error) complete.\n');
end
