function Phase1_Test_Integrity()
    % PHASE1_TEST_INTEGRITY Validates physics and math logic.
    % Passes if H-Matrix matches Finite Difference within tolerance.
    
    clc; close all;
    fprintf('PHASE 1: SYSTEM INTEGRITY CHECK\n');
    
    try
        init_project(); 
    catch
    end
    const = lib_constants();
    t0 = cspice_str2et('2025 DEC 01');

    % Define Test State (Sun-Centered EMO2000)
    st_earth = cspice_spkezr('EARTH', t0, 'J2000', 'NONE', 'SUN');
    r_E_EMO = const.R_EME_EMO * st_earth(1:3);
    v_E_EMO = const.R_EME_EMO * st_earth(4:6);
    
    % Place SC ~400,000 km from Earth (Lunar distance)
    r_sc = r_E_EMO + [300000; 100000; 50000]; 
    v_sc = v_E_EMO + [0.5; 0.5; 0.1];
    X_nom = [r_sc; v_sc; 1.2; 1e-5];
    
    [~, H_anal] = lib_measurements(t0, X_nom, 1, const);
    
    % Finite Difference
    % Large delta (1e-3) used to overcome machine precision limits
    % when subtracting large Sun-Centered vectors.
    delta = 1e-3; 
    H_num = zeros(2,8);
    for i=1:8
        X_p = X_nom; X_p(i) = X_p(i) + delta;
        X_m = X_nom; X_m(i) = X_m(i) - delta;
        [Y_p, ~] = lib_measurements(t0, X_p, 1, const);
        [Y_m, ~] = lib_measurements(t0, X_m, 1, const);
        H_num(:,i) = (Y_p - Y_m)/(2*delta);
    end
    
    diff = max(max(abs(H_anal - H_num)));
    fprintf('Max H-Matrix Error: %e\n', diff);
    
    if diff < 1e-4
        fprintf('RESULT: PASS\n');
    else
        fprintf('RESULT: FAIL\n');
        disp(H_anal - H_num);
    end
end