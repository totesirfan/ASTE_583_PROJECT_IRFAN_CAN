function const = lib_constants()
% LIB_CONSTANTS  Project constants (Sun-centered EMO2000).

%% 1. Time / math
d2r = pi/180;  r2d = 180/pi;
const.deg2rad = d2r;
const.rad2deg = r2d;
const.sec2day = 1/86400;
const.day2sec = 86400;
const.c       = 299792.458;      % km/s
const.AU      = 149597870.7;     % km

%% 2. Bodies
const.mu_E  = 398600;            % km^3/s^2
const.R_E   = 6378;              % km
const.J2    = 0.0010826;
const.we    = 7.292116e-5;       % rad/s
const.mu_S  = 1.32712e11;        % km^3/s^2
const.P_SRP = 4.54e-6;           % N/m^2

%% 3. Spacecraft
const.m_sc    = 200;             % kg
const.A_sc    = 1.5;             % m^2
const.rho_r   = 0.1;
const.k_SRP_0 = 1.0;
const.X0_ref  = [ ...
    1.067623147085261e8; ...
    1.148757045773147e8; ...
   -0.000321627221208e8; ...
  -22.148505873534173;  ...
   18.814312217049999;  ...
   -0.098774507382220];

%% 4. EME2000 -> EMO2000 rotation
eps = (23 + 26/60 + 21.448/3600) * d2r;
c = cos(eps); s = sin(eps);
const.R_EME_EMO = [1 0 0; 0 c s; 0 -s c];

%% 5. Stations (lat/lon in rad)
ids   = 1:4;
names = {'Goldstone','Canberra','Madrid','Antarctica'};
lat_d = [ 35.244352, -35.220919,  40.241355, -80.0 ];
lon_d = [-116.889538, 148.981267, -4.2480085,   0.0 ];

for k = 1:4
    const.stations(k).id   = ids(k);
    const.stations(k).name = names{k};
    const.stations(k).lat  = lat_d(k) * d2r;
    const.stations(k).lon  = lon_d(k) * d2r;
end

%% 6. Time system / Earth rotation
const.epoch_utc_str = '2025 DEC 01 00:00:00.00';
const.t_detect_et   = 817819269.183122;        % ET at detection
gst_deg             = 0 + 10/60 + 43/3600;     % 00:10:43
const.phi_G_detect  = gst_deg * d2r;

%% 7. Maneuvers
const.LTM.date_utc   = '2025 DEC 16 00:00:00.00';
const.LTM.dV         = [0.6931075; -0.8462091; 0.0956979];  % km/s
const.LTM.sigma_exec = 5e-3;                                % km/s

const.LCM.date_utc = '2025 DEC 25 00:00:00.00';
const.dV_budget    = 1.139;                                % km/s

%% 8. Initial covariance (state-only)
sigma_r = 100;          % km
sigma_v = 1e-3;         % km/s
const.P0_state = diag([sigma_r^2 * ones(3,1); ...
                       sigma_v^2 * ones(3,1)]);

%% 9. Augmented covariance (10-state filters)
sigma_kSRP = 1/3;
sigma_bias = 1.0;
sigma_stat = 1 * d2r;

const.P0_aug = blkdiag(const.P0_state, ...
                       sigma_kSRP^2, ...
                       sigma_bias^2, ...
                       sigma_stat^2, ...
                       sigma_stat^2);
end
