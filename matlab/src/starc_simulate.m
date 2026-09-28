function sim = starc_simulate(cfg, p, dt, record)
%STARC_SIMULATE  Periodic steady-state orbit of the one-node radiator model.
%
%   sim = starc_simulate(cfg, p, dt, record)
%
%   Energy balance of a single lumped radiator node (one radiating face,
%   back side adiabatic):
%
%     C dT/dt = Q_int + Q_solar + Q_albedo + Q_earthIR + Q_heater
%               - eps(T) sigma A T^4
%
%     Q_solar   = alpha G A F_sun             (sunlit only)
%     Q_albedo  = alpha a G A F_earth         (sunlit only)
%     Q_earthIR = eps(T) q_IR A F_earth       (whole orbit)
%
%   Kirchhoff's law: for the (grey) IR band, absorptance = emittance, so
%   the Earth-IR term uses the same eps(T) as emission. A cold STARC
%   surface therefore also ABSORBS less Earth IR, which partly offsets
%   its reduced emission. Solar absorptance alpha is a separate property.
%
%   Heater: ideal thermostat. After each RK4 step, if T fell below T_min
%   the heater adds exactly the energy C*(T_min - T) to hold T_min. When
%   the node sits at T_min this equals the instantaneous net heat loss.
%
%   Integration: fixed-step RK4; the sun/eclipse phase is held constant
%   inside a step (dt must divide both the sunlit time and the period, so
%   no step straddles a terminator).
%
%   Periodic steady state: the start-of-orbit temperature T0 is found with
%   a secant iteration on g(T0) = T(end of orbit) - T0; then the final
%   orbit is re-run and reported.
%
%   Inputs
%     cfg    : from STARC_CONFIG (N cases, simulated together, vectorised)
%     p      : from STARC_DEFAULT_PARAMS
%     dt     : time step [s]
%     record : true (default) to keep time histories of the final orbit
%
%   Output struct sim (N-by-1 unless noted)
%     Tpeak_K, Tlow_K         max / min temperature over the orbit
%     E_heater_J              heater energy per orbit
%     E_heater_ecl_J/_sun_J   ... split by orbit phase
%     Qh_peak_W, Qh_avg_W     peak and orbit-average heater power
%     E_in_J, E_out_J         energy in (incl. heater) / radiated per orbit
%     closure                 (E_in - E_out)/E_out  (0 at periodic steady state)
%     converged, iterations   periodic-solve status
%     t_min (1-by-n+1), T_C (N-by-n+1)      temperature history   (record)
%     tq_min (1-by-n), Qh_W (N-by-n)        heater power per step (record)

if nargin < 4, record = true; end

P    = p.orbit.period_min * 60;
tSun = p.orbit.sun_min * 60;
nStep = round(P / dt);
assert(abs(nStep * dt - P) < 1e-9 && abs(round(tSun / dt) * dt - tSun) < 1e-9, ...
    'starc:simulate:dt', 'dt must divide both the sunlit time and the orbit period.');

k = coefficients(cfg, p);

% Initial guess: orbit-average radiative equilibrium at the mid emissivity
fs   = tSun / P;
Qavg = fs * k.QinSun + (1 - fs) * k.QinEcl;
emid = 0.5 * (cfg.eps_c + cfg.eps_h);
Tg   = ((Qavg + emid .* k.IR) ./ (emid .* k.rad)).^0.25;

% Secant iteration on the start-of-orbit temperature
tol = p.num.periodic_tol_K;
Ta = max(Tg, cfg.Tmin);
Pa = one_orbit(Ta, cfg, k, dt, nStep, tSun, false);  ga = Pa - Ta;
Tb = Pa;
Pb = one_orbit(Tb, cfg, k, dt, nStep, tSun, false);  gb = Pb - Tb;
it = 2;
while any(abs(gb) > tol) && it < p.num.periodic_maxit
    denom = gb - ga;
    Tc = Tb - gb .* (Tb - Ta) ./ denom;
    bad = ~isfinite(Tc) | abs(denom) < 1e-12;
    Tc(bad) = Pb(bad);                     % fall back to plain fixed-point step
    Tc = min(max(Tc, cfg.Tmin), 2000);     % keep iterates physical
    Ta = Tb;  ga = gb;
    Tb = Tc;
    Pb = one_orbit(Tb, cfg, k, dt, nStep, tSun, false);  gb = Pb - Tb;
    it = it + 1;
end

% Report the converged orbit
[~, sim] = one_orbit(Tb, cfg, k, dt, nStep, tSun, record);
sim.converged  = abs(gb) <= tol;
sim.iterations = it;
sim.closure    = (sim.E_in_J - sim.E_out_J) ./ sim.E_out_J;
sim.Tpeak_C    = sim.Tpeak_K - p.T0K;
sim.Tlow_C     = sim.Tlow_K  - p.T0K;
sim.Qh_avg_W   = sim.E_heater_J / P;
if record
    sim.t_min  = (0:nStep) * dt / 60;
    sim.tq_min = ((1:nStep) - 0.5) * dt / 60;
    sim.T_C    = sim.T_K - p.T0K;
end
end

% -------------------------------------------------------------------------
function k = coefficients(cfg, p)
% Temperature-independent parts of each heat term [W]
env = p.env;
Qsol = cfg.alpha .* cfg.A .* env.G .* env.F_sun;                 % direct sun
Qalb = cfg.alpha .* cfg.A .* env.albedo .* env.G .* env.F_earth;  % albedo
k.QinSun = cfg.Qsun + Qsol + Qalb;       % Q_int + solar + albedo, sunlit
k.QinEcl = cfg.Qecl;                     % Q_int only, eclipse
k.IR  = env.q_IR .* env.F_earth .* cfg.A;% x eps(T) -> absorbed Earth IR
k.rad = p.sigma .* cfg.A;                % x eps(T) T^4 -> emitted
end

% -------------------------------------------------------------------------
function [Tend, o] = one_orbit(T0, cfg, k, dt, nStep, tSun, record)
N = numel(T0);
C = cfg.C;
T = T0;
Ein = zeros(N, 1);  Eout = Ein;  Eh = Ein;  EhEcl = Ein;
Tpk = T0;  Tlo = T0;  QhPk = Ein;
if record
    Th = zeros(N, nStep + 1);  Th(:, 1) = T0;
    Qh = zeros(N, nStep);
end

for s = 1:nStep
    sunlit = ((s - 0.5) * dt) < tSun;    % phase is constant within the step
    if sunlit, Q0 = k.QinSun; else, Q0 = k.QinEcl; end

    [i1, o1] = flux(T,                           Q0, cfg, k);
    [i2, o2] = flux(T + 0.5 * dt * (i1 - o1) ./ C, Q0, cfg, k);
    [i3, o3] = flux(T + 0.5 * dt * (i2 - o2) ./ C, Q0, cfg, k);
    [i4, o4] = flux(T +       dt * (i3 - o3) ./ C, Q0, cfg, k);
    qin  = (i1 + 2 * i2 + 2 * i3 + i4) / 6;  % RK4-weighted heat in  [W]
    qout = (o1 + 2 * o2 + 2 * o3 + o4) / 6;  % RK4-weighted heat out [W]

    % Same weights for T and for the energy tallies, so the bookkeeping
    % closes exactly: C*(T_end - T_start) = E_in - E_out.
    T    = T + dt * (qin - qout) ./ C;
    Ein  = Ein  + qin  * dt;
    Eout = Eout + qout * dt;

    % Ideal thermostatic heater
    dEh = max(cfg.Tmin - T, 0) .* C;     % [J] this step
    T   = max(T, cfg.Tmin);
    Eh  = Eh + dEh;
    if ~sunlit, EhEcl = EhEcl + dEh; end

    Tpk  = max(Tpk, T);
    Tlo  = min(Tlo, T);
    QhPk = max(QhPk, dEh / dt);
    if record
        Th(:, s + 1) = T;
        Qh(:, s)     = dEh / dt;
    end
end

Tend = T;
if nargout > 1
    o.Tpeak_K        = Tpk;
    o.Tlow_K         = Tlo;
    o.E_heater_J     = Eh;
    o.E_heater_ecl_J = EhEcl;
    o.E_heater_sun_J = Eh - EhEcl;
    o.Qh_peak_W      = QhPk;
    o.E_in_J         = Ein + Eh;
    o.E_out_J        = Eout;
    if record
        o.T_K  = Th;
        o.Qh_W = Qh;
    end
end
end

% -------------------------------------------------------------------------
function [qin, qout] = flux(T, Q0, cfg, k)
e = starc_emissivity(T, cfg.eps_c, cfg.eps_h, cfg.Tt, cfg.w);
qin  = Q0 + e .* k.IR;           % Earth IR absorbed with absorptance eps(T)
qout = e .* k.rad .* T.^4;       % emitted to space
end
