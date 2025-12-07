function dX_out = lib_dynamics(t, X, const)
% LIB_DYNAMICS  Dynamics (state + optional STM) in Sun-centered EMO2000.

%% 1) Parse state
n_total  = numel(X);
calc_stm = (n_total > 20);      % 10-state + 10x10 STM => 110 elements

if calc_stm
    n_state = 10;
    state   = X(1:n_state);
    Phi     = reshape(X(n_state+1:end), n_state, n_state);
else
    n_state = n_total;
    state   = X;
end

r_sc  = state(1:3);
v_sc  = state(4:6);
k_SRP = state(7);    % 8: bias, 9–10: station lat/lon (no dynamics)

%% 2) Ephemeris (Earth in EMO2000)
st_earth = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
r_E_EMO  = const.R_EME_EMO * st_earth(1:3);

%% 3) Accelerations
r_mag     = norm(r_sc);
r_rel     = r_sc - r_E_EMO;
r_rel_mag = norm(r_rel);
r_E_mag   = norm(r_E_EMO);

% Sun point-mass
a_sun   = -const.mu_S * r_sc   / r_mag^3;

% Earth third-body: a = -mu_E( r_rel/|r_rel|^3 + r_E/|r_E|^3 )
a_earth = -const.mu_E * (r_rel / r_rel_mag^3 + r_E_EMO / r_E_mag^3);

% SRP cannonball (scale factor k_SRP)
C_srp = const.P_SRP * const.AU^2 * (1+const.rho_r) * const.A_sc ...
        / const.m_sc / 1000;                   % km/s^2
a_srp = (k_SRP * C_srp) * r_sc / r_mag^3;

a_tot = a_sun + a_earth + a_srp;

%% 4) State derivatives
dX_state        = zeros(n_state,1);
dX_state(1:3)   = v_sc;
dX_state(4:6)   = a_tot;   % params frozen

%% 5) STM derivatives (if requested)
if calc_stm
    A        = zeros(10,10);
    A(1:3,4:6) = eye(3);

    I3 = eye(3);

    % Gravity gradients: G = mu/r^5(3rr^T - r^2 I)
    G_sun   = (const.mu_S / r_mag^5)     * (3*(r_sc*r_sc.')   - r_mag^2*I3);
    G_earth = (const.mu_E / r_rel_mag^5) * (3*(r_rel*r_rel.') - r_rel_mag^2*I3);

    % SRP gradient: C/r^5( r^2 I - 3rr^T )
    G_srp = (k_SRP * C_srp / r_mag^5) * (r_mag^2*I3 - 3*(r_sc*r_sc.'));

    A(4:6,1:3) = G_sun + G_earth + G_srp;

    % da/dk_SRP = a_srp / k_SRP  (guard k_SRP ~ 0)
    if abs(k_SRP) > 0
        A(4:6,7) = a_srp / k_SRP;
    else
        A(4:6,7) = 0;
    end

    dPhi   = A * Phi;
    dX_out = [dX_state; dPhi(:)];
else
    dX_out = dX_state;
end
end
