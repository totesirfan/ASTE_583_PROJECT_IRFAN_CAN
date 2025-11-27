function const = get_constants()
    % GET_CONSTANTS Database of physics and project parameters.
    % Frame Definition: Sun-Centered EMO2000.
    
    %% 1. Time and Math
    const.deg2rad = pi / 180;
    const.rad2deg = 180 / pi;
    const.sec2day = 1 / 86400;
    const.day2sec = 86400;
    const.c = 299792.458; % km/s
    const.AU = 149597870.7; % km

    %% 2. Bodies
    const.mu_E = 398600;       % km^3/s^2
    const.R_E  = 6378;         % km
    const.J2   = 0.0010826;    
    const.we   = 7.292116e-5;  % rad/s
    const.mu_S = 132712000000; % km^3/s^2
    const.P_SRP = 4.54e-6;     % N/m^2
    
    %% 3. Spacecraft
    const.m_sc = 200;          % kg
    const.A_sc = 1.5;          % m^2
    const.rho_r = 0.1;         
    const.k_SRP_0 = 1.0;       
    const.X0_ref = [1.067623147085261e8; 1.148757045773147e8; -0.000321627221208e8; ...
                    -22.148505873534173; 18.814312217049999; -0.098774507382220];

    %% 4. Coordinate Rotations 
    % Rotation from EME2000 (J2000) -> EMO2000 (Integration Frame)
    deg = 23; arcmin = 26; arcsec = 21.448;
    eps = (deg + arcmin/60 + arcsec/3600) * const.deg2rad;
    c = cos(eps); s = sin(eps);
    const.R_EME_EMO = [1, 0, 0; 0, c, s; 0, -s, c];

    %% 5. Stations (Lat/Lon in Radians)
    const.stations(1).id = 1; const.stations(1).name = 'Goldstone';
    const.stations(1).lat = 35.244352 * const.deg2rad;  
    const.stations(1).lon = -116.889538 * const.deg2rad;
    
    const.stations(2).id = 2; const.stations(2).name = 'Canberra';
    const.stations(2).lat = -35.220919 * const.deg2rad; 
    const.stations(2).lon = 148.981267 * const.deg2rad;
    
    const.stations(3).id = 3; const.stations(3).name = 'Madrid';
    const.stations(3).lat = 40.241355 * const.deg2rad;  
    const.stations(3).lon = -4.2480085 * const.deg2rad;
    
    const.stations(4).id = 4; const.stations(4).name = 'Antarctica';
    const.stations(4).lat = -80 * const.deg2rad;        
    const.stations(4).lon = 0 * const.deg2rad;

    % Initial Epoch String
    const.epoch_utc_str = '2025 DEC 01 00:00:00.00';
    
    % GST at J2000 (for rotation calculation)
    const.phi_G_J2000 = 280.46061837 * const.deg2rad;
    
    %% 8. Maneuvers 
    % This is the section that was causing your error
    
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
