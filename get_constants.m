function const = get_constants()
    % GET_CONSTANTS Returns a struct containing all physical constants, 
    % spacecraft parameters, station data, and initial conditions required 
    % for the ASTE 583 Lunar Trailblazer project.
    %
    % Sources: 
    %   - ASTE583_Project.pdf (Project Description)
    %   - ASTE583_Orbital Perturbations_250916.pdf (Lecture Slides)
    
    %% 1. Time and Math Constants
    const.deg2rad = pi / 180;
    const.rad2deg = 180 / pi;
    const.sec2day = 1 / 86400;
    const.day2sec = 86400;
    
    % Speed of Light (Required for Doppler) - Standard Value
    const.c = 299792.458; % km/s
    
    % Astronomical Unit (AU) - Standard Value used for SRP scaling
    const.AU = 149597870.7; % km

    %% 2. Earth Parameters 
    const.mu_E = 398600;       % km^3/s^2
    const.R_E  = 6378;         % km (Earth Radius)
    const.J2   = 0.0010826;    % Earth J2
    const.we   = 7.292116e-5;  % rad/s (Earth rotation rate)
    
    %% 3. Sun Parameters 
    const.mu_S = 132712000000; % km^3/s^2
    
    % Solar Radiation Pressure at 1 AU 
    % Source: Slide 14, "Orbital Perturbations"
    const.P_SRP = 4.54e-6;     % N/m^2
    
    %% 4. Spacecraft Parameters 
    const.m_sc = 200;          % kg (Spacecraft Mass)
    const.A_sc = 1.5;          % m^2 (Cross-sectional Area)
    
    % Reflectivity coefficient
    % Note: Project text implies simple model. Lecture Slide 14 suggests 
    % rho = 0.1 is reasonable for solar arrays.
    const.rho_r = 0.1;         
    
    % Apriori SRP Scale Factor
    const.k_SRP_0 = 1.0;       % Nominal value
    const.sigma_kSRP = 1/3;    % 1-sigma uncertainty (3-sigma = 100%)

    %% 5. Coordinate Transformations 
    % Obliquity of the Ecliptic (epsilon)
    % Value: 23 deg 26' 21.448"
    deg = 23;
    arcmin = 26;
    arcsec = 21.448;
    epsilon_deg = deg + arcmin/60 + arcsec/3600;
    const.epsilon = epsilon_deg * const.deg2rad; % Radians
    
    % Rotation Matrix from EME2000 (ECI) to EMO2000
    % R_EME_EMO = [1 0 0; 0 cos(eps) sin(eps); 0 -sin(eps) cos(eps)]
    c_eps = cos(const.epsilon);
    s_eps = sin(const.epsilon);
    const.R_EME_EMO = [1,     0,      0;
                       0, c_eps,  s_eps;
                       0, -s_eps, c_eps];
    
    % Note: To go from EMO to EME, use transpose(const.R_EME_EMO)

    %% 6. Ground Stations
    % Coordinates given in Latitude / East Longitude (Degrees)
    % Stored in Radians for computation
    
    % Station 1: Goldstone
    const.stations(1).id = 1;
    const.stations(1).name = 'Goldstone';
    const.stations(1).lat = 35.244352 * const.deg2rad;
    const.stations(1).lon = -116.889538 * const.deg2rad;
    
    % Station 2: Canberra
    const.stations(2).id = 2;
    const.stations(2).name = 'Canberra';
    const.stations(2).lat = -35.220919 * const.deg2rad;
    const.stations(2).lon = 148.981267 * const.deg2rad;
    
    % Station 3: Madrid
    const.stations(3).id = 3;
    const.stations(3).name = 'Madrid';
    const.stations(3).lat = 40.241355 * const.deg2rad;
    const.stations(3).lon = -4.2480085 * const.deg2rad;
    
    % Station 4: Antarctica (Secret Location)
    const.stations(4).id = 4;
    const.stations(4).name = 'Antarctica';
    const.stations(4).lat = -80 * const.deg2rad;
    const.stations(4).lon = 0 * const.deg2rad;
    
    % Station Noise Sigmas
    % DSN (Stations 1-3)
    const.sigma_range_DSN = 1e-3;      % km (1 m)
    const.sigma_rr_DSN    = 0.1e-6;    % km/s (0.1 mm/s)
    
    % Antarctica (Station 4)
    const.sigma_range_Ant = 10e-3;     % km (10 m)
    const.sigma_rr_Ant    = 1e-6;      % km/s (1 mm/s)

    %% 7. Initial Conditions (at Detection) 
    % Epoch: 1 DEC 2025 00:00:00.00 UTC
    % NOTE: You must calculate the ET (Ephemeris Time) for this string 
    % using MICE (cspice_str2et) in your main script.
    const.epoch_utc_str = '2025 DEC 01 00:00:00.00';
    
    % Initial Greenwich Sidereal Time (phi_G0)
    % Value: 00:10:43 GST
    % Conversion to radians: (10 min + 43/60 min) / 60 min/hr = 0.1786 hr
    % 0.1786 hr * 15 deg/hr = 2.6791 deg
    hrs = 0; 
    mins = 10; 
    secs = 43;
    decimal_hours = hrs + mins/60 + secs/3600;
    const.phi_G0 = decimal_hours * 15 * const.deg2rad; % Radians
    
    % Reference State (Sun-Centered EMO2000)
    r0 = [1.067623147085261; 1.148757045773147; -0.000321627221208] * 1e8; % km
    v0 = [-22.148505873534173; 18.814312217049999; -0.098774507382220];     % km/s
    const.X0_ref = [r0; v0];
    
    % Apriori Covariance 
    % 1-sigma values
    sigma_r = 100;   % km
    sigma_v = 1e-3;  % km/s (1 m/s)
    
    % P0 Matrix (6x6 for state) - Diagonal
    const.P0_state = diag([sigma_r^2 * ones(3,1); sigma_v^2 * ones(3,1)]);

    %% 8. Maneuvers 
    % Lunar Targeting Maneuver (LTM)
    const.LTM.date_utc = '2025 DEC 16 00:00:00.00';
    % Delta-V in EMO2000 (km/s)
    const.LTM.dV = [0.6931075; -0.8462091; 0.0956979]; 
    % Execution Error (1-sigma spherical)
    const.LTM.sigma_exec = 5e-3; % km/s (5 m/s)
    
    % Lunar Cleanup Maneuver (LCM)
    const.LCM.date_utc = '2025 DEC 25 00:00:00.00';
    
    % Delta-V Budget
    const.dV_budget = 1.139; % km/s (1139 m/s)

end