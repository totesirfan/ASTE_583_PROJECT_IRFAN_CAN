function Phase2_Prelim_Plots()
% PHASE2_PRELIM_PLOTS
% Prelim OD (0–6 d):
%   - Prefit/postfit RR residuals (Batch & Kalman)
%   - 3σ RTN covariance ellipses (Batch vs Kalman)

    clc; close all;
    fprintf('=== Prelim OD plots (0–6 d) ===\n');

    %% 1. Constants + measurements (0–6 d)
    try, init_project(); end
    C   = lib_constants();
    t0  = C.t_detect_et;

    meas = load_project_measurements();
    idx  = meas.time_sec <= 6*C.day2sec;

    t_sec = meas.time_sec(idx);          % s since detection
    st_id = meas.station_id(idx);
    rr    = meas.rr_kmps(idx);           % km/s

    [t_sec, I] = sort(t_sec);
    st_id = st_id(I);
    rr    = rr(I);

    t_et  = t0 + t_sec;
    t_day = t_sec / C.day2sec;
    n     = numel(t_sec);

    fprintf('Using %d measurements up to %.3f days.\n', n, max(t_day));

    %% 2. Load prelim Batch & Kalman solutions
    Sbatch = load('ASTE583_PrelimBatch_Results.mat');
    pr     = Sbatch.prelim_results;
    Xb0    = pr.X0_batch(:);      % 10×1
    Pb0    = pr.P0_batch;         % 10×10

    if exist('ASTE583_PrelimKalman_Results.mat','file')
        Skal = load('ASTE583_PrelimKalman_Results.mat');
    else
        Skal = load('ASTE583_PrelimKalman_Debug.mat');
    end
    kfr = Skal.kf_results;
    Xk0 = kfr.X0_KF(:);           % 10×1
    Pk0 = kfr.P0_KF;              % 10×10

    % Nominal reference at t0 (for PREFIT)
    s4    = C.stations(4);
    Xref0 = [C.X0_ref; C.k_SRP_0; 0; s4.lat; s4.lon];

    %% 3. Helper: propagate 10-state to all measurement times
    function X_all = prop10(X0)
        [t_u, ~, ib] = unique(t_et);          % unique ETs, back-map
        tspan = [t0; t_u(:)];
        opt   = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [T, X] = ode45(@(t,x) lib_dynamics(t,x,C), tspan, X0, opt);

        Xu = zeros(10, numel(t_u));
        for j = 1:numel(t_u)
            tj = t_u(j);
            k  = find(T <= tj, 1, 'last');
            if T(k) == tj || k == numel(T)
                Xu(:,j) = X(k,:).';
            else
                a = (tj - T(k)) / (T(k+1) - T(k));
                Xu(:,j) = (1-a)*X(k,:).' + a*X(k+1,:).';
            end
        end
        X_all = Xu(:, ib);                    % 10×n in original order
    end

    Xref = prop10(Xref0);
    Xb   = prop10(Xb0);
    Xk   = prop10(Xk0);

    %% 4. RR prefit/postfit residuals (km/s)
    rr_pre_b  = zeros(n,1);
    rr_post_b = zeros(n,1);
    rr_pre_k  = zeros(n,1);
    rr_post_k = zeros(n,1);

    for i = 1:n
        ti = t_et(i);
        si = st_id(i);

        [Y_ref, ~] = lib_measurements(ti, Xref(:,i), si, C);
        [Y_b,   ~] = lib_measurements(ti, Xb(:,i),   si, C);
        [Y_k,   ~] = lib_measurements(ti, Xk(:,i),   si, C);

        % Prefit (same ref for both)
        rr_pre_b(i)  = rr(i) - Y_ref(2);
        rr_pre_k(i)  = rr_pre_b(i);

        % Postfit
        rr_post_b(i) = rr(i) - Y_b(2);
        rr_post_k(i) = rr(i) - Y_k(2);
    end

    %% 5. Plot helper (RR in mm/s) with fixed axes
    function plot_rr(fig_name, t, pre, post)
        rr_pre_mm  = 1e6 * pre;    % km/s -> mm/s
        rr_post_mm = 1e6 * post;

        figure('Name',fig_name, 'Color','w', 'Position',[50 100 1000 400]);

        % Prefit
        subplot(1,2,1); hold on; grid on;
        scatter(t, rr_pre_mm, 4, 'filled');
        yline(0,'k:');
        xlabel('Time since detection (days)');
        ylabel('Prefit range-rate residual (mm/s)');
        title([fig_name ' - prefit']);
        xlim([0 max(t)]);
        ylim([4500 6500]);   % <-- your preferred prefit scale

        % Postfit
        subplot(1,2,2); hold on; grid on;
        scatter(t, rr_post_mm, 4, 'filled');
        yline(0,'k:');
        xlabel('Time since detection (days)');
        ylabel('Postfit range-rate residual (mm/s)');
        title([fig_name ' - postfit']);
        xlim([0 max(t)]);
        ylim([-6 4]);        % <-- your preferred postfit scale (mm/s)
    end

    plot_rr('Prelim Batch RR',  t_day, rr_pre_b, rr_post_b);
    plot_rr('Prelim Kalman RR', t_day, rr_pre_k, rr_post_k);

    %% 6. RTN 3σ covariance ellipses (position block only)
    fprintf('Computing RTN 3σ ellipses at t0...\n');

    rB = Xb0(1:3); vB = Xb0(4:6);
    rK = Xk0(1:3);
    Pxyz_b = Pb0(1:3,1:3);
    Pxyz_k = Pk0(1:3,1:3);

    % RTN transform from inertial r,v (rows: R, T, N)
    function T = TRTN(r, v)
        rhat = r / norm(r);
        h    = cross(r, v);
        nhat = h / norm(h);
        that = cross(nhat, rhat);
        T    = [rhat.'; that.'; nhat.'];
    end

    T      = TRTN(rB, vB);
    PRTN_b = T * Pxyz_b * T.';
    PRTN_k = T * Pxyz_k * T.';

    dRtn = T * (rK - rB);
    dR = dRtn(1); dT = dRtn(2); dN = dRtn(3);

    % 2×2 blocks
    P_RT_b = PRTN_b(1:2,1:2);     P_RT_k = PRTN_k(1:2,1:2);
    P_RN_b = PRTN_b([1 3],[1 3]); P_RN_k = PRTN_k([1 3],[1 3]);
    P_TN_b = PRTN_b(2:3,2:3);     P_TN_k = PRTN_k(2:3,2:3);

    % 3σ ellipse plot helper
    function ellipse2(P2, x0, y0, style)
        [V, D] = eig(P2);
        D      = max(D,0);                 % guard tiny negatives
        s      = 3*sqrt(diag(D));          % 3σ lengths
        th     = linspace(0, 2*pi, 200);
        circ   = [cos(th); sin(th)];
        e      = V * (diag(s) * circ);
        plot(x0 + e(1,:), y0 + e(2,:), style, 'LineWidth',1.5);
    end

    figure('Name','Prelim RTN 3sigma', ...
           'Color','w','Position',[1100 200 1400 450]);

    % RT plane (dR vs dT)
    subplot(1,3,1); hold on; grid on; axis equal;
    plot(0,  0,  'bo', 'MarkerFaceColor','b');
    plot(dR, dT, 'rs', 'MarkerFaceColor','r');
    ellipse2(P_RT_b, 0,  0,  'b-');
    ellipse2(P_RT_k, dR, dT, 'r-');
    xlabel('\DeltaR (km)','Interpreter','Tex');
    ylabel('\DeltaT (km)','Interpreter','Tex');
    title('RT plane');
    legend('Batch','Kalman','Batch 3\sigma','Kalman 3\sigma','Location','best');

    % RN plane (dR vs dN)
    subplot(1,3,2); hold on; grid on; axis equal;
    plot(0,  0,  'bo', 'MarkerFaceColor','b');
    plot(dR, dN, 'rs', 'MarkerFaceColor','r');
    ellipse2(P_RN_b, 0,  0,  'b-');
    ellipse2(P_RN_k, dR, dN, 'r-');
    xlabel('\DeltaR (km)','Interpreter','Tex');
    ylabel('\DeltaN (km)','Interpreter','Tex');
    title('RN plane');

    % TN plane (dT vs dN)
    subplot(1,3,3); hold on; grid on; axis equal;
    plot(0,  0,  'bo', 'MarkerFaceColor','b');
    plot(dT, dN, 'rs', 'MarkerFaceColor','r');
    ellipse2(P_TN_b, 0,  0,  'b-');
    ellipse2(P_TN_k, dT, dN, 'r-');
    xlabel('\DeltaT (km)','Interpreter','Tex');
    ylabel('\DeltaN (km)','Interpreter','Tex');
    title('TN plane');

    fprintf('Prelim plots complete.\n');
end
