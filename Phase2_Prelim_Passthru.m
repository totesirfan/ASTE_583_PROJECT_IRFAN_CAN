function Phase2_Prelim_Passthru()
% Phase2_Prelim_Passthru
% Passthru test for PRELIM OD solutions (Batch & Kalman) using 6–14 day data.

    clear; clc; close all;
    fprintf('=== PHASE 2: PRELIM OD PASSTHRU (6–14 DAYS) ===\n');

    %% 1. Setup
    try, init_project(); end
    const = lib_constants();
    t0_et = const.t_detect_et;

    %% 2. Load prelim Batch & Kalman

    fprintf('Loading preliminary Batch and Kalman results...\n');

    % ---- Batch ----
    if ~exist('ASTE583_PrelimBatch_Results.mat','file')
        error('ASTE583_PrelimBatch_Results.mat not found.');
    end
    Sb = load('ASTE583_PrelimBatch_Results.mat');
    if isfield(Sb,'prelim_results') && isfield(Sb.prelim_results,'X0_batch')
        X0_batch = Sb.prelim_results.X0_batch(:);
    elseif isfield(Sb,'X0_batch')
        X0_batch = Sb.X0_batch(:);
    else
        error('X0_batch not found in ASTE583_PrelimBatch_Results.mat.');
    end

    % ---- Kalman ----
    if exist('ASTE583_PrelimKalman_Results.mat','file')
        Sk = load('ASTE583_PrelimKalman_Results.mat');
    elseif exist('ASTE583_PrelimKalman_Debug.mat','file')
        Sk = load('ASTE583_PrelimKalman_Debug.mat');
    else
        error('No prelim Kalman MAT file found.');
    end
    if isfield(Sk,'kf_results') && isfield(Sk.kf_results,'X0_KF')
        X0_kf = Sk.kf_results.X0_KF(:);
    elseif isfield(Sk,'X0_KF')
        X0_kf = Sk.X0_KF(:);
    else
        error('X0_KF not found in prelim Kalman MAT file.');
    end

    if numel(X0_batch) ~= 10 || numel(X0_kf) ~= 10
        error('Both preliminary solutions must be 10x1 vectors.');
    end

    %% 3. Load measurements and select 6–14 day window

    fprintf('Loading measurements (0–14 days)...\n');
    meas = load_project_measurements();

    t_sec   = meas.time_sec(:);      % s since detection
    st_id   = meas.station_id(:);
    rho_km  = meas.range_km(:);      % km (NaN if missing)
    rr_kmps = meas.rr_kmps(:);       % km/s

    t_days_all = t_sec / const.day2sec;
    t_et_all   = t0_et + t_sec;

    idx = (t_days_all > 6.0) & (t_days_all <= 14.0);
    t_days = t_days_all(idx);
    t_et   = t_et_all(idx);
    st_id  = st_id(idx);
    rho_km = rho_km(idx);
    rr_kmps= rr_kmps(idx);

    [t_days, sort_idx] = sort(t_days);
    t_et   = t_et(sort_idx);
    st_id  = st_id(sort_idx);
    rho_km = rho_km(sort_idx);
    rr_kmps= rr_kmps(sort_idx);

    n_meas = numel(t_et);
    fprintf('Passthru: %d measurements between %.3f and %.3f days.\n', ...
            n_meas, min(t_days), max(t_days));

    %% 4. Solution list

    sols(1).name   = 'Prelim Batch';
    sols(1).X0_ref = X0_batch;

    sols(2).name   = 'Prelim Kalman';
    sols(2).X0_ref = X0_kf;

    st_colors = {'r','g','b','k'};
    st_names  = {'Goldstone','Canberra','Madrid','Antarctica'};

    %% 5. Loop over solutions

    for s = 1:numel(sols)
        label = sols(s).name;
        X0    = sols(s).X0_ref;

        fprintf('\n--- PASSTHRU: %s ---\n', label);

        % 5.1 Propagate t0 -> max(t_et)
        t_end   = max(t_et);
        dt_days = (t_end - t0_et) * const.sec2day;
        fprintf('Propagating from t0 to %.3f days (%.1f h)...\n', ...
                dt_days, dt_days*24);

        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [T_grid, X_grid] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                                 [t0_et t_end], X0, opts);

        X_meas = interp1(T_grid, X_grid, t_et, 'spline');  % n_meas x 10

        % 5.2 Compute residuals
        rho_res = NaN(n_meas,1);
        rr_res  = NaN(n_meas,1);

        for k = 1:n_meas
            Xk = X_meas(k,:).';
            [Y, ~] = lib_measurements(t_et(k), Xk, st_id(k), const);  % [rho; rho_dot]
            if ~isnan(rho_km(k))
                rho_res(k) = rho_km(k) - Y(1);
            end
            rr_res(k) = rr_kmps(k) - Y(2);
        end

        rr_res_mm = rr_res * 1e6;  % km/s -> mm/s

        % 5.3 RMS metrics (overall)
        rho_valid = ~isnan(rho_res);
        rr_valid  = ~isnan(rr_res_mm);

        if any(rho_valid)
            rho_rms = sqrt(mean(rho_res(rho_valid).^2));
        else
            rho_rms = NaN;
        end
        if any(rr_valid)
            rr_rms = sqrt(mean(rr_res_mm(rr_valid).^2));
        else
            rr_rms = NaN;
        end

        fprintf('Overall RMS Range      : %.4f km\n',  rho_rms);
        fprintf('Overall RMS Range-rate : %.2f mm/s\n', rr_rms);

        % 5.3 Per-station RMS
        for st = 1:4
            j_r  = (st_id == st) & rho_valid;
            j_rr = (st_id == st) & rr_valid;

            if any(j_r)
                rms_r = sqrt(mean(rho_res(j_r).^2));
            else
                rms_r = NaN;
            end

            if any(j_rr)
                rms_rr = sqrt(mean(rr_res_mm(j_rr).^2));
            else
                rms_rr = NaN;
            end

            fprintf('  Station %d (%-10s): RMS Range = %7.4f km,  RMS RR = %7.2f mm/s\n', ...
                    st, st_names{st}, rms_r, rms_rr);
        end

        % 5.4 Plots
        % Range
        figure('Name',['Passthru Range Residuals - ' label], ...
               'Color','w','Position',[100 100 900 400]);
        hold on; grid on;
        for st = 1:4
            j = (st_id == st) & rho_valid;
            if any(j)
                scatter(t_days(j), rho_res(j), 10, st_colors{st}, 'filled');
            end
        end
        yline(0,'k-');
        xlabel('Time since detection (days)');
        ylabel('Range O - C (km)');
        title(sprintf('Passthru Range Residuals (6-14 d) - %s', label), ...
              'Interpreter','none');
        legend(st_names,'Location','best');

        % Range-rate
        figure('Name',['Passthru Doppler Residuals - ' label], ...
               'Color','w','Position',[100 550 900 400]);
        hold on; grid on;
        for st = 1:4
            j = (st_id == st) & rr_valid;
            if any(j)
                scatter(t_days(j), rr_res_mm(j), 10, st_colors{st}, 'filled');
            end
        end
        yline(0,'k-');
        xlabel('Time since detection (days)');
        ylabel('Range-rate O - C (mm/s)');
        title(sprintf('Passthru Doppler Residuals (6-14 d) - %s', label), ...
              'Interpreter','none');
        legend(st_names,'Location','best');

        sols(s).rho_res = rho_res;
        sols(s).rr_res  = rr_res;
    end

    fprintf('\nPassthru for preliminary Batch and Kalman solutions complete.\n');
end
