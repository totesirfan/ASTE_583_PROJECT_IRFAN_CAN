function Phase2_FinalKalman_PostProcess()
% Post-process FINAL CKF:
%  1) Plot pre/post range & range-rate residuals
%  2) Propagate covariance to LTM and compute RTN sigmas

    clear; clc; close all;
    fprintf('=== FINAL KALMAN POST-PROCESS (Residuals + LTM RTN Covariance) ===\n');

    % 1. Setup
    try, init_project(); catch, end
    const = lib_constants();
    t0_et = const.t_detect_et;

    % 2. Load FINAL CKF file
    fname = 'ASTE583_FinalKalman_Results.mat';
    if ~isfile(fname)
        error('File %s not found. Run Phase2_Final_Kalman first.', fname);
    end
    S = load(fname);

    % Expect fields saved by Phase2_Final_Kalman
    X0_final_kalman = S.X0_final_kalman(:);   % 10x1
    P0_final_kalman = S.P0_final_kalman;      % 10x10
    t_et_all  = S.t_et(:);                    % Nx1 ET times
    st_id     = S.st_id(:);                   % Nx1 station ids
    has_range = logical(S.has_range(:));      % Nx1

    rho_prefit  = S.rho_prefit(:);            % km
    rho_postfit = S.rho_postfit(:);           % km
    rr_prefit   = S.rr_prefit(:);             % km/s
    rr_postfit  = S.rr_postfit(:);            % km/s

    n_meas = numel(t_et_all);
    t_days = (t_et_all - t0_et) * const.sec2day;

    fprintf('Loaded FINAL CKF: %d measurements, span %.3f days.\n', ...
            n_meas, max(t_days));

    %% 3. Range-rate residuals (mm/s), per station

    fprintf('Plotting range-rate residuals...\n');
    figure('Name','FINAL Kalman RR Residuals','Color','w','Position',[100 100 1000 600]);

    colors        = lines(4);
    station_ids   = 1:4;
    station_names = {const.stations.name};

    % Prefit
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

    % Postfit
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

    %% 4. Range residuals (km), only where range exists

    fprintf('Plotting range residuals...\n');
    hasR = has_range & ~isnan(rho_prefit);

    figure('Name','FINAL Kalman Range Residuals','Color','w','Position',[150 150 1000 600]);

    % Prefit
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

    % Postfit
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

    %% 5. Covariance propagation to LTM and RTN sigmas

    fprintf('Propagating covariance to LTM and computing RTN sigmas...\n');

    t_LTM_et = cspice_str2et(const.LTM.date_utc);
    n_state  = 10;
    Phi0     = eye(n_state);
    Z0       = [X0_final_kalman(:); Phi0(:)];

    opt_ode = odeset('RelTol',1e-10,'AbsTol',1e-9);
    [~, Z_prop] = ode45(@(t,Z) lib_dynamics(t,Z,const), ...
                        [t0_et, t_LTM_et], Z0, opt_ode);

    Z_LTM   = Z_prop(end,:).';
    X_LTM   = Z_LTM(1:n_state);
    Phi_LTM = reshape(Z_LTM(n_state+1:end), n_state, n_state);

    P_LTM = Phi_LTM * P0_final_kalman * Phi_LTM.';   % 10x10
    P_LTM_state = P_LTM(1:6,1:6);                    % r,v block

    r = X_LTM(1:3);
    v = X_LTM(4:6);

    r_hat = r / norm(r);
    h     = cross(r, v);
    n_hat = h / norm(h);
    t_hat = cross(n_hat, r_hat);

    C_RTN_I = [r_hat.'; t_hat.'; n_hat.'];          % rows: R,T,N
    R6      = blkdiag(C_RTN_I, C_RTN_I);
    P_RTN   = R6 * P_LTM_state * R6.';

    sigma_R = sqrt(P_RTN(1,1));
    sigma_T = sqrt(P_RTN(2,2));
    sigma_N = sqrt(P_RTN(3,3));

    fprintf('\n--- LTM RTN POSITION UNCERTAINTY (1σ) ---\n');
    fprintf('Sigma_R (radial)      : %.3f km (3σ = %.3f km)\n', sigma_R, 3*sigma_R);
    fprintf('Sigma_T (along-track) : %.3f km (3σ = %.3f km)\n', sigma_T, 3*sigma_T);
    fprintf('Sigma_N (cross-track) : %.3f km (3σ = %.3f km)\n', sigma_N, 3*sigma_N);
    fprintf('LTM epoch: %s\n', const.LTM.date_utc);
    fprintf('\nDone. Use RR/Range residual plots + sigma_T at LTM in your write-up.\n');
end
