function init_project()
    % INIT_PROJECT Configures MICE and loads kernels for Irfan's setup.
    % Run this once at the start of every MATLAB session.
    
    % 1. Add MICE to Path
    mice_root = 'C:\Users\irfan\OneDrive\Documents\MATLAB\mice';
    if exist(mice_root, 'dir')
        addpath(genpath(mice_root));
        fprintf('SUCCESS: MICE added to MATLAB path.\n');
    else
        error('ERROR: Could not find MICE at: %s', mice_root);
    end
    
    % 2. Load Kernels
    kernels = {
        'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\de440s.bsp';
        'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\naif0012.tls.pc'
    };
    
    fprintf('Kernals Loaded\n');
    for i = 1:length(kernels)
        k = kernels{i};
        if isfile(k)
            cspice_furnsh(k);
        else
            warning('MISSING: %s', k);
        end
    end
end