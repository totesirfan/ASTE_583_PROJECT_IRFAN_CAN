% Phase1_Debug_Check.m
% Verifies trajectory/measurement model against Professor's email debug values.
% 
% Reference Values from Email:
% Time: 6.060763888888889 days
% Range: 5.894311564735891e+07 km
% Rate:  1.885427030618933 km/s

clear; clc;

% 1. Initialize
try
    init_project();
catch
    warning('init_project failed or already initialized.');
end
const = lib_constants();

fprintf('--- DEBUG VERIFICATION ---\n');

% 2. Define Debug Time
dt_days = 6.060763888888889;
t_debug = const.t_detect_et + dt_days * 86400;

% 3. Propagate Reference Trajectory to Debug Time
% State: [r; v; k_SRP=1.0; bias=0.0] - Matches Phase 1 State
X0 = [const.X0_ref; 1.0; 0.0]; 

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

fprintf('Propagating to t = %.9f days...\n', dt_days);
[~, X_out] = ode45(@(t,x) lib_dynamics(t, x, const), [const.t_detect_et t_debug], X0, options);

X_debug = X_out(end, :)';

% 4. Compute Measurement for Station 1 (Goldstone)
stn_id = 1; 
[G_x, ~] = lib_measurements(t_debug, X_debug, stn_id, const);

range_calc = G_x(1);
rate_calc  = G_x(2);

% 5. Compare with Professor's Values
range_ref = 5.894311564735891e+07;
rate_ref  = 1.885427030618933;

range_diff = range_calc - range_ref;
rate_diff  = rate_calc - rate_ref;

fprintf('\nRESULTS:\n');
fprintf('Parameter      | Your Calc          | Reference          | Diff\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('Range (km)     | %18.9e | %18.9e | %.2e\n', range_calc, range_ref, range_diff);
fprintf('Rate (km/s)    | %18.9f | %18.9f | %.2e\n', rate_calc, rate_ref, rate_diff);

% 6. Auto-Check
if abs(range_diff) < 50 && abs(rate_diff) < 1e-4
    fprintf('\n[PASS] Your physics model matches the reference solution!\n');
    fprintf('       (Small differences are expected due to integration tolerances)\n');
else
    fprintf('\n[FAIL] Discrepancy detected. Check GST angle logic or Station Coordinates.\n');
end