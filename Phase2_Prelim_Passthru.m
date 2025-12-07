function Phase2_Prelim_Passthru()
% Passthru test for PRELIM OD (Batch & Kalman) on 6–14 day data.
% Uses:
%   ASTE583_PrelimBatch_Results.mat  -> prelim_results.X0_batch
%   ASTE583_PrelimKalman_Results.mat -> kf_results.X0_KF

    clear; clc; close all;
    fprintf('=== PRELIM PASSTHRU: 6–14 DAYS ===\n');

    % Adjustable y-limits (set [] for auto)
    rho_ylim = [2580 2640];   % km
    rr_ylim  = [64 76];       % mm/s

    % --- Setup ---
    try, init_project(); end
    const = lib_constants();
    t0_et = const.t_detect_et;

    % --- Load Batch ---
    bf = 'ASTE583_PrelimBatch_Results.mat';
    if ~isfile(bf), error('%s missing. Run Prelim Batch.', bf); end
    Sb = load(bf);
    if ~isfield(Sb,'prelim_results') || ~isfield(Sb.prelim_results,'X0_batch')
        error('prelim_results.X0_batch missing in %s', bf);
    end
    X0_batch = Sb.prelim_results.X0_batch(:);

    % --- Load Kalman ---
    kf = 'ASTE583_PrelimKalman_Results.mat';
    if ~isfile(kf), error('%s missing. Run Prelim Kalman.', kf); end
    Sk = load(kf);
    if ~isfield(Sk,'kf_results') || ~isfield(Sk.kf_results,'X0_KF')
        error('kf_results.X0_KF missing in %s', kf);
    end
    X0_kf = Sk.kf_results.X0_KF(:);

    if numel(X0_batch) ~= 10 || numel(X0_kf) ~= 10
        error('Prelim states must be 10x1.');
    end

    % --- Measurements: 6–14 d subset ---
    fprintf('Loading measurements...\n');
    meas    = load_project_measurements();
    t_sec   = meas.time_sec(:);
    st_id   = meas.station_id(:);
    rho_km  = meas.range_km(:);
    rr_kmps = meas.rr_kmps(:);

    t_days_all = t_sec / const.day2sec;
    t_et_all   = t0_et + t_sec;

    idx    = (t_days_all > 6) & (t_days_all <= 14);
    t_days = t_days_all(idx);
    t_et   = t_et_all(idx);
    st_id  = st_id(idx);
    rho_km = rho_km(idx);
    rr_kmps= rr_kmps(idx);

    [t_days, I] = sort(t_days);
    t_et   = t_et(I);
    st_id  = st_id(I);
    rho_km = rho_km(I);
    rr_kmps= rr_kmps(I);

    n_meas = numel(t_et);
    fprintf('Using %d meas from %.3f to %.3f d.\n', ...
            n_meas, min(t_days), max(t_days));

    % --- Solution list ---
    sols(1).name = 'Batch';   sols(1).X0 = X0_batch;
    sols(2).name = 'Kalman';  sols(2).X0 = X0_kf;

    st_names = {'Goldstone','Canberra','Madrid','Antarctica'};
    opts     = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % --- Propagate + residuals ---
    for s = 1:2
        label = sols(s).name;
        X0    = sols(s).X0;

        fprintf('\n-- %s passthru --\n', label);
        t_end   = max(t_et);
        span_d  = (t_end - t0_et)*const.sec2day;
        fprintf('Prop to %.3f d (%.1f h)\n', span_d, 24*span_d);

        [T, X]   = ode45(@(t,x) lib_dynamics(t,x,const), [t0_et t_end], X0, opts);
        X_meas   = interp1(T, X, t_et, 'spline');  % n×10
        rho_res  = NaN(n_meas,1);
        rr_res   = NaN(n_meas,1);

        for k = 1:n_meas
            Xk = X_meas(k,:).';
            [Y, ~] = lib_measurements(t_et(k), Xk, st_id(k), const);
            if ~isnan(rho_km(k)), rho_res(k) = rho_km(k) - Y(1); end
            rr_res(k) = rr_kmps(k) - Y(2);
        end

        rr_res_mm = rr_res*1e6;

        v_rho = ~isnan(rho_res);
        v_rr  = ~isnan(rr_res_mm);

        rho_rms = sqrt(mean(rho_res(v_rho).^2));
        rr_rms  = sqrt(mean(rr_res_mm(v_rr).^2));
        fprintf('RMS rho = %.4f km, RMS rr = %.2f mm/s\n', rho_rms, rr_rms);

        for st = 1:4
            j_r  = (st_id == st) & v_rho;
            j_rr = (st_id == st) & v_rr;

            rms_r  = sqrt(mean(rho_res(j_r).^2));
            rms_rr = sqrt(mean(rr_res_mm(j_rr).^2));
            fprintf('  St %d %-10s: rho = %7.4f km, rr = %7.2f mm/s\n', ...
                    st, st_names{st}, rms_r, rms_rr);
        end

        sols(s).rho_res = rho_res;
        sols(s).rr_res  = rr_res;
    end

    % --- Combined plots: Batch vs Kalman ---
    rho_b = sols(1).rho_res;
    rho_k = sols(2).rho_res;
    rr_b  = sols(1).rr_res*1e6;
    rr_k  = sols(2).rr_res*1e6;

    vb  = ~isnan(rho_b); vk = ~isnan(rho_k);
    vrb = ~isnan(rr_b);  vrk = ~isnan(rr_k);

    % Range
    figure('Name','Passthru Range Residuals 6-14 d', ...
           'Color','w','Position',[100 100 900 400]);
    hold on; grid on;
    scatter(t_days(vb), rho_b(vb), 8, 'b', 'filled');
    scatter(t_days(vk), rho_k(vk), 8, 'r', 'filled');
    yline(0,'k-');
    xlabel('Time since detection (days)');
    ylabel('Range O - C (km)');
    title('Passthru Range Residuals 6-14 d');
    legend({'Batch','Kalman'},'Location','best');
    if ~isempty(rho_ylim), ylim(rho_ylim); end

    % Range-rate
    figure('Name','Passthru Doppler Residuals 6-14 d', ...
           'Color','w','Position',[100 550 900 400]);
    hold on; grid on;
    scatter(t_days(vrb), rr_b(vrb), 8, 'b', 'filled');
    scatter(t_days(vrk), rr_k(vrk), 8, 'r', 'filled');
    yline(0,'k-');
    xlabel('Time since detection (days)');
    ylabel('Range-rate O - C (mm/s)');
    title('Passthru Doppler Residuals 6-14 d');
    legend({'Batch','Kalman'},'Location','best');
    if ~isempty(rr_ylim), ylim(rr_ylim); end

    fprintf('\nPrelim passthru complete.\n');
end
