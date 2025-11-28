function Phase2_Test_STM()
    % PHASE2_TEST_STM Validates the Jacobian (A-Matrix) logic.
    % Compares Analytical A-matrix vs. Finite Difference.
    
    clc; close all;
    fprintf('PHASE 2: STM DYNAMICS INTEGRITY CHECK\n');
    
    % 1. Setup
    try
        init_project(); 
    catch
    end
    const = lib_constants();
    t0 = cspice_str2et(const.epoch_utc_str);
    
    % 2. Define Nominal State (Sun-Centered)
    % We use a state slightly perturbed from X0 to ensure gradients are non-zero
    X_nom = [const.X0_ref; 
             1.2;   % k_SRP
             0.5;   % bias
             const.stations(4).lat; 
             const.stations(4).lon];
         
    % Create Augmented State for Identity STM
    Phi_identity = eye(10);
    X_aug = [X_nom; Phi_identity(:)];
    
    % 3. Compute Analytical dX and A
    % We call lib_dynamics. The derivative of the STM part is A*Phi.
    % Since Phi is Identity, dPhi = A.
    dX_aug_anal = lib_dynamics(t0, X_aug, const);
    
    % Extract the A matrix from the augmented derivative
    % The last 100 elements are dPhi(:). Reshape to 10x10.
    A_anal = reshape(dX_aug_anal(11:end), 10, 10);
    
    % 4. Compute Numerical A (Finite Difference)
    delta = 1e-4;
    A_num = zeros(10, 10);
    
    % Loop through 10 state elements
    for i = 1:10
        % Perturb +
        X_p = X_nom; X_p(i) = X_p(i) + delta;
        dX_aug_p = lib_dynamics(t0, [X_p; Phi_identity(:)], const);
        dX_state_p = dX_aug_p(1:10);
        
        % Perturb -
        X_m = X_nom; X_m(i) = X_m(i) - delta;
        dX_aug_m = lib_dynamics(t0, [X_m; Phi_identity(:)], const);
        dX_state_m = dX_aug_m(1:10);
        
        % Central Difference
        A_num(:, i) = (dX_state_p - dX_state_m) / (2 * delta);
    end
    
    % 5. Compare
    % We only care about rows 1-6 (Dynamics). Rows 7-10 are parameters (0=0).
    A_diff = A_anal(1:6, :) - A_num(1:6, :);
    max_err = max(max(abs(A_diff)));
    
    fprintf('Max Jacobian Error (Dynamics Rows): %e\n', max_err);
    
    % Tolerance check
    if max_err < 1e-6
        fprintf('RESULT: PASS (Jacobian is consistent)\n');
    else
        fprintf('RESULT: FAIL\n');
        disp('Difference Matrix (Top 6 rows):');
        disp(A_diff);
        
        fprintf('\nAnalytical A (Top 6 rows):\n');
        disp(A_anal(1:6, :));
        
        fprintf('\nNumerical A (Top 6 rows):\n');
        disp(A_num(1:6, :));
    end
end