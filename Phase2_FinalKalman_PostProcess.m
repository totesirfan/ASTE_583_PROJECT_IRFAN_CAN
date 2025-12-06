function Phase2_FinalKalman_PostProcess()
% Phase2_FinalKalman_PostProcess.m
% Post-processing for FINAL Kalman OD:
%   1) Plot converged prefit/postfit residuals (range & range rate)
%   2) Propagate covariance to LTM and compute RTN sigmas

    clear; clc; close all;
    fprintf('=== PHASE 2: FINAL KALMAN POST-PROCESSING (Residuals + LTM Covariance) ===\n');

    %% 1. Setup & load constants
    try
        init_project();
    catch
    end
    const  = lib_constants();
    t0_et  = const.t_detect_et;

    %% 2. Load FINAL Kalman results
    fname = 'ASTE583_FinalKalman_Results.mat';
    if ~isfile(fname)
        error('File %s not found. Run Phase2_Final_Kalman first.', fname);
    end

    S = load(fname);

    requiredFields = {'X0_final_kalman','P0_final_kalman', ...
                      't_et_all','st_id','has_range', ...
                      'rho_prefit','rho_postfit', ...
                      'rr_prefit','rr_postfit'};
    for f = requiredFields
        if ~isfield(S, f{1})
            error('Field "%s" missing from %s. Re-run Phase2_Final_Kalman.', f{1}, fname);
        end
    end

    X0_final_kalman = S.X0_final_kalman(:);    % 10x1
    P0_final_kalman = S.P0_final_kalman;       % 10x10
    t_et_all  = S.t_et_all(:);                 % Nx1 ET times
    st_id     = S.st_id(:);                    % Nx1 station ids
    has_range = S.has_range(:) ~= 0;           % logical

    rho_prefit  = S.rho_prefit(:);            % Nx1 (km)
    rho_postfit = S.rho_postfit(:);           % Nx1 (km)
    rr_prefit   = S.rr_prefit(:);             % Nx1 (km/s)
    rr_postfit  = S.rr_postfit(:);            % Nx1 (km/s)

    n_meas = numel(t_et_all);

    % Time in days from detection epoch
    t_days = (t_et_all - t0_et) * const.sec2day;

    fprintf('Loaded FINAL Kalman solution: %d measurements spanning %.3f days.\n', ...
             n_meas, max(t_days));

    %% 3. Plot converged MMR range-rate residuals (prefit & postfit)

    fprintf('Plotting converged range-rate prefits/postfits...\n');

    figure('Name','FINAL Kalman RR Residuals','Color','w','Position',[100 100 1000 600]);

    % Colours by station
    colors = lines(4);
    station_ids = 1:4;
    station_names = {const.stations.name};

    % --- Prefit RR ---
    subplot(2,1,1); hold on; grid on;
    for i = 1:numel(station_ids)
        sid = station_ids(i);
        idx = (st_id == sid);
        if any(idx)
            plot(t_days(idx), rr_prefit(idx)*1e6, '.', ...
                 'MarkerSize', 4, 'Color', colors(i,:));
        end
    end
    yline(0,'k-','LineWidth',1);
    xlabel('Time since detection (days)');
    ylabel('RR Prefit (mm/s)');
    title('FINAL Kalman: Range-Rate Prefit Residuals');
    legend(station_names{:}, 'Location','bestoutside');

    % --- Postfit RR ---
    subplot(2,1,2); hold on; grid on;
    for i = 1:numel(station_ids)
        sid = station_ids(i);
        idx = (st_id == sid);
        if any(idx)
            plot(t_days(idx), rr_postfit(idx)*1e6, '.', ...
                 'MarkerSize', 4, 'Color', colors(i,:));
        end
    end
    yline(0,'k-','LineWidth',1);
    xlabel('Time since detection (days)');
    ylabel('RR Postfit (mm/s)');
    title('FINAL Kalman: Range-Rate Postfit Residuals');

    %% 4. Plot converged range residuals (only where range is present)

    fprintf('Plotting converged range prefits/postfits...\n');

    hasR = has_range & ~isnan(rho_prefit);

    figure('Name','FINAL Kalman Range Residuals','Color','w','Position',[150 150 1000 600]);

    % --- Prefit Range ---
    subplot(2,1,1); hold on; grid on;
    for i = 1:numel(station_ids)
        sid = station_ids(i);
        idx = hasR & (st_id == sid);
        if any(idx)
            plot(t_days(idx), rho_prefit(idx), '.', ...
                 'MarkerSize', 4, 'Color', colors(i,:));
        end
    end
    yline(0,'k-','LineWidth',1);
    xlabel('Time since detection (days)');
    ylabel('Range Prefit (km)');
    title('FINAL Kalman: Range Prefit Residuals');
    legend(station_names{:}, 'Location','bestoutside');

    % --- Postfit Range ---
    subplot(2,1,2); hold on; grid on;
    for i = 1:numel(station_ids)
        sid = station_ids(i);
        idx = hasR & (st_id == sid);
        if any(idx)
            plot(t_days(idx), rho_postfit(idx), '.', ...
                 'MarkerSize', 4, 'Color', colors(i,:));
        end
    end
    yline(0,'k-','LineWidth',1);
    xlabel('Time since detection (days)');
    ylabel('Range Postfit (km)');
    title('FINAL Kalman: Range Postfit Residuals');

    %% 5. Propagate covariance to LTM and compute RTN sigmas

    fprintf('Propagating FINAL Kalman covariance to LTM and computing RTN sigmas...\n');

    % LTM epoch (ET)
    t_LTM_et = cspice_str2et(const.LTM.date_utc);

    n_state = 10;
    Phi0 = eye(n_state);
    Z0   = [X0_final_kalman(:); Phi0(:)];

    options_ode = odeset('RelTol',1e-10,'AbsTol',1e-9);

    [T_prop, Z_prop] = ode45(@(t,Z) lib_dynamics(t, Z, const), ...
                             [t0_et, t_LTM_et], Z0, options_ode);

    Z_LTM = Z_prop(end,:).';
    X_LTM = Z_LTM(1:n_state);
    Phi_flat = Z_LTM(n_state+1:end);
    Phi_LTM  = reshape(Phi_flat, n_state, n_state);

    % Propagate covariance: P(LTM) = Phi * P0 * Phi^T
    P0 = P0_final_kalman;
    P_LTM = Phi_LTM * P0 * Phi_LTM.';

    % Extract 6x6 position/velocity block
    P_LTM_state = P_LTM(1:6,1:6);

    % Build RTN frame at LTM (EMO2000)
    r = X_LTM(1:3);
    v = X_LTM(4:6);

    r_hat = r / norm(r);
    h     = cross(r, v);
    n_hat = h / norm(h);
    t_hat = cross(n_hat, r_hat);

    C_RTN_I = [r_hat.'; t_hat.'; n_hat.'];      % rows: R, T, N in inertial coords
    R6 = blkdiag(C_RTN_I, C_RTN_I);

    % Covariance in RTN
    P_RTN = R6 * P_LTM_state * R6.';

    sigma_R = sqrt(P_RTN(1,1));
    sigma_T = sqrt(P_RTN(2,2));
    sigma_N = sqrt(P_RTN(3,3));

    fprintf('\n--- LTM RTN POSITION UNCERTAINTY (1σ) ---\n');
    fprintf('Sigma_R (radial)      : %.3f km (3σ = %.3f km)\n', sigma_R, 3*sigma_R);
    fprintf('Sigma_T (along-track) : %.3f km (3σ = %.3f km)\n', sigma_T, 3*sigma_T);
    fprintf('Sigma_N (cross-track) : %.3f km (3σ = %.3f km)\n', sigma_N, 3*sigma_N);
    fprintf('LTM epoch: %s\n', const.LTM.date_utc);

    fprintf('\nPost-processing complete. Use the RR/Range plots + sigma_T at LTM in your write-up.\n');
end
