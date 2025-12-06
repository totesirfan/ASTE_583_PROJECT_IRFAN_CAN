function meas = load_project_measurements()
% LOAD_PROJECT_MEASUREMENTS
%   Load LTB radiometric data from the two project CSVs and return a
%   unified measurement struct.
%
% Files (must be in the current folder or on the MATLAB path):
%   - ASTE583_Project_LTB_Measurements_0-6D_Truth.csv
%   - ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv
%
% Output struct fields:
%   meas.time_sec    : [N×1] time since detection epoch, seconds
%   meas.station_id  : [N×1] station ID (1,2,3,4)
%   meas.range_km    : [N×1] two-way range (km), NaN where not available
%   meas.rr_kmps     : [N×1] two-way range-rate (km/s)

    file_0_6  = 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv';
    file_6_14 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';

    % --- 0–6 days: Antarctica only, Doppler-only, col 3 is RR (km/s) ---
    T0 = readtable(file_0_6,  'VariableNamingRule','preserve');

    t0   = T0{:,1};   % seconds from detection epoch
    st0  = T0{:,2};   % station ID (should be 4)
    rr0  = T0{:,3};   % actually range-rate [km/s]
    rng0 = nan(size(rr0));  % no range in this file

    % --- 6–14 days: DSN + Antarctica, proper Range + Range-rate ---
    T1  = readtable(file_6_14, 'VariableNamingRule','preserve');

    t1   = T1{:,1};   % seconds from detection epoch
    st1  = T1{:,2};   % station ID
    rng1 = T1{:,3};   % range [km]
    rr1  = T1{:,4};   % range-rate [km/s]

    % --- Combine into a single struct and sort by time ---
    meas.time_sec   = [t0;  t1];
    meas.station_id = [st0; st1];
    meas.range_km   = [rng0; rng1];
    meas.rr_kmps    = [rr0;  rr1];

    [meas.time_sec, idx] = sort(meas.time_sec);
    meas.station_id = meas.station_id(idx);
    meas.range_km   = meas.range_km(idx);
    meas.rr_kmps    = meas.rr_kmps(idx);
end
