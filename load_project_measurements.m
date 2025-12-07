function meas = load_project_measurements()
% LOAD_PROJECT_MEASUREMENTS  Merge 0–6 d and 6–14 d LTB radiometric data.

    f0 = 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv';
    f1 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';

    % 0–6 days: Antarctica Doppler only (RR in col 3)
    T0 = readtable(f0,'VariableNamingRule','preserve');
    t0  = T0{:,1};
    st0 = T0{:,2};
    rr0 = T0{:,3};
    rng0 = nan(height(T0),1);        % no range in this span

    % 6–14 days: DSN + Antarctica, range + RR
    T1 = readtable(f1,'VariableNamingRule','preserve');
    t1  = T1{:,1};
    st1 = T1{:,2};
    rng1 = T1{:,3};
    rr1  = T1{:,4};

    % Stack and sort by time
    t  = [t0;  t1];
    st = [st0; st1];
    rng = [rng0; rng1];
    rr  = [rr0;  rr1];

    [t, idx] = sort(t);
    st  = st(idx);
    rng = rng(idx);
    rr  = rr(idx);

    meas.time_sec   = t;
    meas.station_id = st;
    meas.range_km   = rng;
    meas.rr_kmps    = rr;
end
