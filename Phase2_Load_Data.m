function obs_data = Phase2_Load_Data()
    % PHASE2_LOAD_DATA  Parses measurement CSVs for the Lunar Trailblazer project.
    % Updates: Converts relative time (sec) to Ephemeris Time (ET).
    
    %% 1. Setup
    try
        init_project(); 
    catch
    end
    
    % Define Epoch (Matches lib_constants)
    epoch_str = '2025 DEC 01 00:00:00.00';
    t0_et = cspice_str2et(epoch_str);
    fprintf('Reference Epoch: %s (ET: %.4f)\n', epoch_str, t0_et);

    % Measurement noise sigmas
    sig_range_dsn = 1e-3;      % 1 m
    sig_doppler_dsn = 1e-7;    % 0.1 mm/s
    sig_range_ant = 10e-3;     % 10 m
    sig_doppler_ant = 1e-6;    % 1 mm/s
    
    %% 2. Load File 1: Days 0-6 (Safe Mode)
    filename1 = 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv';
    if isfile(filename1)
        fprintf('Loading %s...\n', filename1);
        raw1 = readmatrix(filename1);
        
        t1   = raw1(:, 1);
        st1  = raw1(:, 2);
        rr1  = raw1(:, 3); % Range Rate
        
        n1 = length(t1);
        data1 = repmat(struct('t',0,'ID',0,'type','', 'value',0, 'sigma',0, 'bias_flag',0), n1, 1);
        
        for k = 1:n1
            % FIX: Add Epoch to Time
            data1(k).t = t0_et + t1(k); 
            data1(k).ID = st1(k);
            data1(k).type = 'Doppler';
            data1(k).value = rr1(k);
            data1(k).bias_flag = 1; 
            
            if st1(k) == 4
                data1(k).sigma = sig_doppler_ant;
            else
                data1(k).sigma = sig_doppler_dsn;
            end
        end
    else
        data1 = [];
    end
    
    %% 3. Load File 2: Days 6-14 (Normal Ops)
    filename2 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';
    if isfile(filename2)
        fprintf('Loading %s...\n', filename2);
        raw2 = readmatrix(filename2);
        
        t2   = raw2(:, 1);
        st2  = raw2(:, 2);
        r2   = raw2(:, 3); 
        rr2  = raw2(:, 4); 
        
        n2 = length(t2);
        data2 = repmat(struct('t',0,'ID',0,'type','', 'value',0, 'sigma',0, 'bias_flag',0), 2*n2, 1);
        
        idx = 1;
        for k = 1:n2
            % A. Range
            data2(idx).t = t0_et + t2(k); % FIX: Add Epoch
            data2(idx).ID = st2(k);
            data2(idx).type = 'Range';
            data2(idx).value = r2(k);
            data2(idx).bias_flag = 0;
            
            if st2(k) == 4
                data2(idx).sigma = sig_range_ant;
            else
                data2(idx).sigma = sig_range_dsn;
            end
            idx = idx + 1;
            
            % B. Doppler
            data2(idx).t = t0_et + t2(k); % FIX: Add Epoch
            data2(idx).ID = st2(k);
            data2(idx).type = 'Doppler';
            data2(idx).value = rr2(k);
            data2(idx).bias_flag = 0;
            
            if st2(k) == 4
                data2(idx).sigma = sig_doppler_ant;
            else
                data2(idx).sigma = sig_doppler_dsn;
            end
            idx = idx + 1;
        end
    else
        data2 = [];
    end
    
    %% 4. Combine
    obs_data = [data1; data2];
    [~, sort_idx] = sort([obs_data.t]);
    obs_data = obs_data(sort_idx);
    
    fprintf('SUCCESS: Loaded %d measurements.\n', length(obs_data));
end