function setup_project()
    % SETUP_PROJECT Configures MICE and loads specific kernels for Irfan's setup
    
    % --- 1. Add MICE to MATLAB Path ---
    mice_root = 'C:\Users\irfan\OneDrive\Documents\MATLAB\mice';
    
    if exist(mice_root, 'dir')
        addpath(genpath(mice_root));
        savepath; % Saves this path for future sessions
        fprintf('SUCCESS: MICE added to MATLAB path.\n');
    else
        error('ERROR: Could not find MICE at: %s', mice_root);
    end
    
    % --- 2. Load SPICE Kernels ---
    % Define the specific kernels you downloaded
    kernels = {
        'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\de440s.bsp';
        'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\naif0012.tls.pc'
    };
    
    % Load them
    fprintf('Loading Kernels...\n');
    for i = 1:length(kernels)
        k = kernels{i};
        if isfile(k)
            cspice_furnsh(k);
            fprintf('  LOADED: %s\n', k);
        else
            warning('  MISSING: %s', k);
        end
    end
    
    % --- 3. Verify Installation ---
    try
        % Try to convert a date to check if it works
        et = cspice_str2et('2025 DEC 01');
        fprintf('\nVERIFICATION PASSED: MICE is working. ET = %f\n', et);
    catch ME
        fprintf('\nVERIFICATION FAILED: %s\n', ME.message);
    end
end