function Phase3_ManeuverExecution()
% Predict RR signature of the LTM maneuver using FINAL 10-state Kalman OD.
% - Propagate no-LTM and LTM trajectories
% - For each station: compute RR_noLTM, RR_LTM, dRR = RR_LTM - RR_noLTM
% - Print max |dRR| and plot all stations in a single figure.

    clear; clc; close all;
    fprintf('=== LTM maneuver RR residual prediction (FINAL Kalman OD) ===\n');

    %% 1. Setup + final Kalman state
    try
        init_project();
    catch
    end
    const = lib_constants();
    t0_et  = const.t_detect_et;
    t_LTM  = cspice_str2et(const.LTM.date_utc);

    fname = 'ASTE583_FinalKalman_Results.mat';
    if ~isfile(fname)
        error('File %s not found. Run Phase2_Final_Kalman first.', fname);
    end
    S  = load(fname);
    X0 = S.X0_final_kalman(:);
    if numel(X0) ~= 10
        error('Final Kalman state must be 10x1; got %d.', numel(X0));
    end

    %% 2. Time grid around LTM (±1 h in 1 min steps)
    t_span_min = 60;                        % half-width [min]
    dt_min     = 1;                         % step [min]
    t_rel_min  = (-t_span_min:dt_min:t_span_min).';  % [N x 1], minutes
    t_eval     = t_LTM + 60*t_rel_min;      % ET [s]

    t_end = max(t_eval);
    fprintf('Propagating no-LTM trajectory to %.3f days...\n', ...
            (t_end - t0_et)*const.sec2day);

    opts = odeset('RelTol',1e-11,'AbsTol',1e-11);

    %% 3. No-LTM propagation
    [T_no, X_no] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                         [t0_et, t_end], X0, opts);

    %% 4. LTM propagation: t0 -> LTM, apply dV, LTM -> t_end
    fprintf('Propagating with LTM at %s...\n', const.LTM.date_utc);

    [T1, X1] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                     [t0_et, t_LTM], X0, opts);

    X_LTM = X1(end,:).';
    X_LTM(4:6) = X_LTM(4:6) + const.LTM.dV(:);

    [T2, X2] = ode45(@(t,X) lib_dynamics(t,X,const), ...
                     [t_LTM, t_end], X_LTM, opts);

    %% 5. Compute RR_noLTM, RR_LTM for all stations
    nT    = numel(t_eval);
    nSt   = 4;

    % Color convention (by station ID):
    % 1: Goldstone  -> red
    % 2: Canberra   -> magenta
    % 3: Madrid     -> green
    % 4: Antarctica -> black
    colors = [
        0.89 0.10 0.11;   % Goldstone (red)
        0.23 0.49 0.77;   % Canberra (magenta)
        0.20 0.63 0.17;   % Madrid (green)
        0.60 0.31 0.64    % Antarctica (black)
    ];

    RR_no_all  = zeros(nT, nSt);   % km/s
    RR_LTM_all = zeros(nT, nSt);   % km/s

    for st = 1:nSt
        name = const.stations(st).name;

        for k = 1:nT
            tk = t_eval(k);

            % no-LTM state
            X_no_k = interp1(T_no, X_no, tk, 'pchip').';

            % LTM state (piecewise)
            if tk <= t_LTM
                X_LTM_k = interp1(T1, X1, tk, 'pchip').';
            else
                X_LTM_k = interp1(T2, X2, tk, 'pchip').';
            end

            [Y_no,  ~] = lib_measurements(tk, X_no_k,  st, const);
            [Y_LTM, ~] = lib_measurements(tk, X_LTM_k, st, const);

            RR_no_all(k,st)  = Y_no(2);
            RR_LTM_all(k,st) = Y_LTM(2);
        end

        % ΔRR in mm/s for this station
        dRR_mmps = (RR_LTM_all(:,st) - RR_no_all(:,st))*1e6;
        [dRR_max_abs, idx_max] = max(abs(dRR_mmps));
        t_max_min = t_rel_min(idx_max);

        fprintf('  Station %d (%s)\n', st, name);
        fprintf('    Max |RR residual| = %.2f mm/s at t = %+5.2f min from LTM.\n', ...
                dRR_max_abs, t_max_min);
    end

    dRR_all_mmps = (RR_LTM_all - RR_no_all)*1e6;  % mm/s

    %% 6. Single figure: all stations
    figure('Name','LTM RR signature', ...
           'Color','w','Position',[100 100 900 600]);

    % (a) RR_noLTM and RR_LTM
    subplot(2,1,1); hold on; grid on;
    for st = 1:nSt
        plot(t_rel_min, RR_no_all(:,st)*1e6,  '--',  ...
             'LineWidth',1.1,'Color',colors(st,:));
        plot(t_rel_min, RR_LTM_all(:,st)*1e6, '-', ...
             'LineWidth',1.1,'Color',colors(st,:));
    end
    yline(0,'k:');
    xlabel('Time from LTM (min)','Interpreter','Tex');
    ylabel('Range-Rate Residuals (mm/s)','Interpreter','none');
    title('Range-Rate Residuals noLTM vs LTM','Interpreter','Tex');
    leg = cell(1,2*nSt);
    for st = 1:nSt
        name = const.stations(st).name;
        leg{2*st-1} = sprintf('%s noLTM', name);
        leg{2*st}   = sprintf('%s LTM'  , name);
    end
    legend(leg,'Location','best','Interpreter','Tex');

    % (b) ΔRR per station
    subplot(2,1,2); hold on; grid on;
    for st = 1:nSt
        plot(t_rel_min, dRR_all_mmps(:,st), '-', ...
             'LineWidth',1.1,'Color',colors(st,:));
    end
    yline(0,'k:');
    xlabel('Time from LTM (min)','Interpreter','none');
    ylabel('\Delta Range-Rate Residuals (mm/s)','Interpreter','Tex');
    title('\Delta Range-Rate Residuals (LTM - noLTM)','Interpreter','Tex');
    leg_d = cell(1,nSt);
    for st = 1:nSt
        leg_d{st} = const.stations(st).name;
    end
    legend(leg_d,'Location','best','Interpreter','none');

    fprintf('Maneuver execution prediction (FINAL Kalman OD) complete.\n');
end
