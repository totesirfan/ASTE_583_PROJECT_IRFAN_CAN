function [Y_comp, H] = lib_measurements(t, X, station_id, const)
% LIB_MEASUREMENTS  Range / range-rate and H-matrix in Sun-centered EMO2000.

%% 1) Unpack state
nX   = numel(X);
r_sc = X(1:3);
v_sc = X(4:6);

% k_SRP present but not used in measurement geometry
if nX >= 7, k_SRP = X(7); %#ok<NASGU>
else,       k_SRP = 1.0;  %#ok<NASGU>
end

% Doppler bias
if nX >= 8, bias = X(8); else, bias = 0.0; end

% Are we estimating station-4 lat/lon?
est_station = (nX >= 10) && (station_id == 4);

%% 2) Earth ephemeris & Earth-fixed rotation
stE     = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
r_E_Sun = const.R_EME_EMO * stE(1:3);
v_E_Sun = const.R_EME_EMO * stE(4:6);

phi_G = const.phi_G_detect + const.we * (t - const.t_detect_et);
c = cos(phi_G); s = sin(phi_G);

R_ECF_ECI = [ c -s  0;
              s  c  0;
              0  0  1];
R_ECF_EMO = const.R_EME_EMO * R_ECF_ECI;

%% 3) Station position / velocity
if est_station
    lat = X(9);  lon = X(10);
else
    st  = const.stations(station_id);
    lat = st.lat; lon = st.lon;
end

cl = cos(lat); sl = sin(lat);
cd = cos(lon); sd = sin(lon);

r_sta_ecf = const.R_E * [cl*cd; cl*sd; sl];
w_vec     = [0; 0; const.we];
v_sta_ecf = cross(w_vec, r_sta_ecf);

r_sta_emo = R_ECF_EMO * r_sta_ecf;
v_sta_emo = R_ECF_EMO * v_sta_ecf;

r_sta = r_E_Sun + r_sta_emo;
v_sta = v_E_Sun + v_sta_emo;

%% 4) Observations (range, range-rate + bias)
rho_vec   = r_sc - r_sta;
rho       = norm(rho_vec);
u_rho     = rho_vec / rho;
rel_vel   = v_sc - v_sta;
rho_dot   = dot(rho_vec, rel_vel) / rho;

t_bias_end = const.t_detect_et + 6*const.day2sec;
if t < t_bias_end
    y_bias = bias;
    H_b    = 1;
else
    y_bias = 0;
    H_b    = 0;
end

Y_comp = [rho;
          rho_dot + y_bias];

%% 5) H-matrix (state partials)
% Range wrt spacecraft state
H_r   = u_rho.';          % 1x3
H_v   = zeros(1,3);       % 1x3

% Range-rate wrt spacecraft state
H_rr_r = (rel_vel.' - rho_dot*u_rho.') / rho;
H_rr_v = u_rho.';

% Base 2x8: [r(3) v(3) k_SRP bias]
H = [H_r,    H_v,    0, 0;
     H_rr_r, H_rr_v, 0, H_b];

%% 6) Station partials (lat / lon)
if est_station
    % dr/d(lat,lon) in ECF
    dr_dlat = const.R_E * [-sl*cd; -sl*sd; cl];
    dr_dlon = const.R_E * [-cl*sd;  cl*cd;  0];

    % Transform to EMO and get dv/d(lat,lon)
    drI_dlat = R_ECF_EMO * dr_dlat;
    drI_dlon = R_ECF_EMO * dr_dlon;
    dvI_dlat = R_ECF_EMO * cross(w_vec, dr_dlat);
    dvI_dlon = R_ECF_EMO * cross(w_vec, dr_dlon);

    % Range partials
    H_rho_lat = -u_rho.' * drI_dlat;
    H_rho_lon = -u_rho.' * drI_dlon;

    % Range-rate partials
    d_rr_dr_sta = -H_rr_r;
    d_rr_dv_sta = -u_rho.';
    H_rr_lat    = d_rr_dr_sta*drI_dlat + d_rr_dv_sta*dvI_dlat;
    H_rr_lon    = d_rr_dr_sta*drI_dlon + d_rr_dv_sta*dvI_dlon;

    H = [H, [H_rho_lat; H_rr_lat], [H_rho_lon; H_rr_lon]];
elseif nX >= 10
    % Keep 2x10 shape when state has lat/lon but station is not estimated
    H = [H, zeros(2,2)];
end
end
