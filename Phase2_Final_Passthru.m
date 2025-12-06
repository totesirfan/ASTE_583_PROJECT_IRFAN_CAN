function Phase2_Final_Passthru()
% PHASE2_FINAL_PASSTHRU
% Passthru test for Preliminary OD solutions (Batch & Kalman) using
% 6–14 day data that were NOT included in the prelim fits.
%
% For each solution:
%   - Propagate 10-state from detection epoch through last passthru time
%   - Compute range and range-rate O-C residuals
%   - Print overall and per-station RMS metrics
%   - Plot residuals vs time (6–14 days)

    clear; clc; close all;
    fprintf('=== PHASE 2: FINAL OD PASSTHRU USING PRELIM SOLUTIONS ===\n');

    %% 1. Setup & constants
    try
        init_project();
    catch
    end
    const = lib_constants();
    t0_et = const.t_detect_et;

    %% 2. Load prelim Batch & Kalman solutions
    fprintf('Loading preliminary Batch and Kalman results...\n');

    % ----- Batch -----
    if ~exist('ASTE583_PrelimBatch_Results.mat','file')
        error('ASTE583_PrelimBatch_Results.mat not found.');
    end
    Sbatch = load('ASTE583_PrelimBatch_Results.mat');
    if isfield(Sbatch,'prelim_results') && isfield(Sbatch.prelim_results,'X0_batch')
        X0_batch = Sbatch.prelim_results.X0_batch;
        P0_batch = Sbatch.prelim_results.P0_batch;
    elseif isfield(Sbatch,'X0_batch')
        X0_batch = Sbatch.X0_batch;
        P0_batch = Sbatch.P0_batch;
    else
        error('Could not find X0_batch in ASTE583_PrelimBatch_Results.mat');
    end

    % ----- Kalman -----
    if exist('ASTE583_PrelimKalman_Results.mat','file')
        Skal = load('ASTE583_PrelimKalman_Results.mat');
    elseif exist('ASTE583_PrelimKalman_Debug.mat','file')
        Skal = load('ASTE583_PrelimKalman_Debug.mat');
    else
        error('No prelim Kalman MAT file found.');
    end
    if isfield(Skal,'kf_results') && isfield(Skal.kf_results,'X0_KF')
        X0_kalman = Skal.kf_results.X0_KF;
        P0_kalman = Skal.kf_results.P0_KF;
    elseif isfield(Skal,'X0_KF')
        X0_kalman = Skal.X0_KF;
        P0_kalman = Skal.P0_KF;
    else
        error('Could not find X0_KF in prelim Kalman MAT file.');
    end

    if numel(X0_batch) ~= 10 || numel(X0_kalman) ~= 10
        error('Both preliminary solutions must be 10-state vectors.');
    end

    %% 3. Load full measurements and select 6–14 day passthru window
    fprintf('Loading measurements (0–14 days)...\n');
    meas = load_project_measurements();  % must include both CSVs

    t_sec  = meas.time_sec;          % seconds since detection
    st_id  = meas.station_id;        % station IDs
    rho_km = meas.range_km;          % NaN if no range
    rr_kmps= meas.rr_kmps;           % km/s

    t_days_all = t_sec / const.day2sec;
    t_et_all   = t0_et + t_sec;

    % Passthru = post-prelim data: 6–14 days
    idx_passthru = (t_days_all > 6.0) & (t_days_all <= 14.0);

    t_days = t_days_all(idx_passthru);
    t_et   = t_et_all(idx_passthru);
    st_id  = st_id(idx_passthru);
    rho_km = rho_km(idx_passthru);
    rr_kmps= rr_kmps(idx_passthru);

    % Sort by time just in case
    [t_days, sort_idx] = sort(t_days);
    t_et   = t_et(sort_idx);
    st_id  = st_id(sort_idx);
    rho_km = rho_km(sort_idx);
    rr_kmps= rr_kmps(sort_idx);

    n_meas = numel(t_et);
    fprintf('Passthru: using %d measurements between %.3f and %.3f days.\n', ...
            n_meas, min(t_days), max(t_days));

    %% 4. Define solutions array (Batch & Kalman) for loop
    solutions(1).name   = 'Prelim Batch';
    solutions(1).X0_ref = X0_batch;
    solutions(1).P0_ref = P0_batch;

    solutions(2).name   = 'Prelim Kalman';
    solutions(2).X0_ref = X0_kalman;
    solutions(2).P0_ref = P0_kalman;

    % Station labels/colors for plots
    st_colors = {'r','g','b','k'};
    st_names  = {'Goldstone','Canberra','Madrid','Antarctica'};

    %% 5. Loop over both solutions and perform passthru
    for s = 1:numel(solutions)
        sol_label = solutions(s).name;
        X0_ref    = solutions(s).X0_ref;

        fprintf('\n--- PASSTHRU USING %s ---\n', sol_label);

        % 5.1 Propagate 10-state from t0 to last passthru time
        t_end = max(t_et);
        fprintf('Propagating from t0 to %.3f days (%.1f hours)...\n', ...
                (t_end - t0_et)/const.day2sec, (t_end - t0_et)/3600);

        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [T_grid, X_grid] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                                 [t0_et, t_end], X0_ref, opts);

        % Interpolate state at measurement times
        X_meas = interp1(T_grid, X_grid, t_et, 'spline');  % n_meas x 10

        % 5.2 Compute O-C residuals
        rho_resid = NaN(n_meas,1);
        rr_resid  = NaN(n_meas,1);

        for k = 1:n_meas
            Xk = X_meas(k,:).';
            tk = t_et(k);
            st = st_id(k);

            [Y_comp, ~] = lib_measurements(tk, Xk, st, const); % [rho; rho_dot]

            if ~isnan(rho_km(k))
                rho_resid(k) = rho_km(k) - Y_comp(1);
            end
            if ~isnan(rr_kmps(k))
                rr_resid(k)  = rr_kmps(k) - Y_comp(2);
            end
        end

        %% 5.3 Print useful scalar metrics (overall + per station)
        rr_resid_mmps = rr_resid * 1e6; % km/s -> mm/s

        % Overall RMS (ignoring NaN)
        overall_rho_rms = sqrt(mean(rho_resid(~isnan(rho_resid)).^2));
        overall_rr_rms  = sqrt(mean(rr_resid_mmps(~isnan(rr_resid_mmps)).^2));

        fprintf('Overall RMS Range Residual:      %.4f km\n',  overall_rho_rms);
        fprintf('Overall RMS Range-Rate Residual: %.2f mm/s\n', overall_rr_rms);

        % Per-station RMS
        for st = 1:4
            idx_st_rho = (st_id == st) & ~isnan(rho_resid);
            idx_st_rr  = (st_id == st) & ~isnan(rr_resid_mmps);

            if any(idx_st_rho)
                rms_rho_st = sqrt(mean(rho_resid(idx_st_rho).^2));
            else
                rms_rho_st = NaN;
            end

            if any(idx_st_rr)
                rms_rr_st = sqrt(mean(rr_resid_mmps(idx_st_rr).^2));
            else
                rms_rr_st = NaN;
            end

            fprintf('  Station %d (%-10s): RMS Range = %7.4f km,  RMS RR = %7.2f mm/s\n', ...
                    st, st_names{st}, rms_rho_st, rms_rr_st);
        end

        %% 5.4 Plots: range & Doppler residuals vs time (6–14 d)

        % ---- Range residuals ----
        figure('Name',['Passthru Range Residuals – ' sol_label], ...
               'Color','w','Position',[100 100 900 400]);
        hold on; grid on;
        for st = 1:4
            idx_st = (st_id == st) & ~isnan(rho_resid);
            if any(idx_st)
                scatter(t_days(idx_st), rho_resid(idx_st), 10, st_colors{st}, 'filled');
            end
        end
        yline(0,'k-');
        xlabel('Time since detection (days)');
        ylabel('Range Residual O - C (km)');
        title(sprintf('Passthru Range Residuals (6–14 days) – %s', sol_label));
        legend(st_names,'Location','best');
        hold off;

        % ---- Range-rate residuals ----
        figure('Name',['Passthru Doppler Residuals – ' sol_label], ...
               'Color','w','Position',[100 550 900 400]);
        hold on; grid on;
        for st = 1:4
            idx_st = (st_id == st) & ~isnan(rr_resid_mmps);
            if any(idx_st)
                scatter(t_days(idx_st), rr_resid_mmps(idx_st), 10, st_colors{st}, 'filled');
            end
        end
        yline(0,'k-');
        xlabel('Time since detection (days)');
        ylabel('Range-Rate Residual O - C (mm/s)');
        title(sprintf('Passthru Doppler Residuals (6–14 days) – %s', sol_label));
        legend(st_names,'Location','best');
        hold off;

        % Store results if you want to save later
        solutions(s).rho_resid = rho_resid;
        solutions(s).rr_resid  = rr_resid;
    end

    fprintf('\nPassthru for both preliminary solutions complete.\n');
end
