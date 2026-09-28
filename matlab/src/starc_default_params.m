function p = starc_default_params()
%STARC_DEFAULT_PARAMS  Every model input in one struct.
%
%   p = starc_default_params() returns the default inputs for the STARC
%   heater-power study. Edit values here, or override fields after the
%   call (e.g. p.coat.STARC.Tt_C = 15) before running the study.
%
%   Source tags used in the comments:
%     [WP-T3]  Aetara White Paper 01 (Sep 2026), Table 3 "Common assumptions"
%     [WP-S1]  White paper, Customer Economics, Step 1 (surface effect)
%     [WP-S2]  White paper, Customer Economics, Step 2 (value of a heater watt)
%     [WP-A]   White paper, STARC properties quoted as company data [A]
%     [STD]    Standard physical constant or textbook LEO value
%     [ASSUME] Modelling choice made for this tool; not from the paper
%
%   FIRST-ORDER PUBLIC-DATA ESTIMATE. Not validated against STARC test data.

%% Physical constants
p.sigma = 5.670374419e-8;   % Stefan-Boltzmann constant [W/m^2/K^4]   [STD]
p.T0K   = 273.15;           % 0 degC in kelvin                         [STD]

%% Orbit (reference LEO)
p.orbit.period_min = 95;    % orbit period [min]                       [WP-T3]
p.orbit.sun_min    = 60;    % sunlit part [min]; eclipse = rest (35)   [WP-T3]

%% Thermal environment
p.env.G       = 1361;       % solar constant [W/m^2]                   [STD]
p.env.albedo  = 0.30;       % Earth albedo [-]                         [STD]
p.env.q_IR    = 237;        % Earth outgoing longwave IR [W/m^2]       [STD]
p.env.F_sun   = 0.25;       % effective solar projection factor on the
                            % radiator face while sunlit [-]           [ASSUME]
p.env.F_earth = 0.30;       % radiator-to-Earth view factor [-]        [ASSUME]

%% Spacecraft / radiator node
p.sc.panel_areal_mass = 5;  % radiator panel system mass [kg/m^2]      [WP-T3]
p.sc.equip_mass       = 10; % equipment lumped with the node [kg]      [ASSUME]
p.sc.cp               = 900;% specific heat, aluminium [J/kg/K]        [STD]
p.sc.T_max_C          = 40; % hot-case limit used to size area [degC]  [ASSUME]
p.sc.T_min_C          = -10;% cold-case heater set point [degC]        [ASSUME]

%% Internal dissipation Q_int [W] per case and orbit phase
% Hot case: payload on all orbit. Cold case: payload off (bus only).
% Set *_sun and *_ecl differently to model e.g. payload on only in sunlight.
p.hot.Qint_sun  = 250;      %                                          [ASSUME]
p.hot.Qint_ecl  = 250;      %                                          [ASSUME]
p.cold.Qint_sun = 50;       %                                          [ASSUME]
p.cold.Qint_ecl = 50;       %                                          [ASSUME]

%% Coatings
% Every coating uses the same logistic model
%   eps(T) = eps_c + (eps_h - eps_c) / (1 + exp(-(T - Tt)/w)),
% so a fixed coating is just eps_c = eps_h (Tt_C and w are then unused).
% Solar absorptance alpha = 0.10 for all three (spec input).
p.coat.FEP   = mk('Silver FEP', 0.10, 0.78, 0.78,  0, 1); % eps 0.75-0.82 [WP-T3, Sheldahl E3]
p.coat.OSR   = mk('Glass OSR',  0.10, 0.88, 0.88,  0, 1); % typical OSR   [STD]
p.coat.STARC = mk('STARC',      0.10, 0.20, 0.92, 10, 4); % range ~0.20-0.95 [WP-A];
                                                          % Tt, w [ASSUME, user input]

p.baseline = 'OSR';         % fixed coating used as reference in sweep & Monte Carlo

%% Economics (white paper Table 3 / Step 2)
p.econ.battery_Wh_per_kg = 150;        % usable battery specific energy  [WP-T3]
p.econ.array_W_per_kg    = 60;         % solar-array specific power      [WP-T3]
p.econ.recharge_eff      = 0.85;       % array -> battery -> load        [WP-T3]
p.econ.launch_usd_per_kg = 7000;       % rideshare price                 [WP-T3, E1]
p.econ.usd_per_W         = [150 800];  % value of a heater watt:
                                       % constellation .. science        [WP-S2]
p.econ.ft2_per_m2        = 10.7639;    % unit conversion                 [STD]

%% Numerics
p.num.dt             = 1;        % RK4 step for headline results [s]
p.num.dt_study       = 5;        % RK4 step for sweep + Monte Carlo [s]
                                 % (run_starc_study checks it against dt)
p.num.periodic_tol_K = 1e-4;     % |T(end of orbit) - T(start)| to accept
p.num.periodic_maxit = 60;       % max orbits in the periodic solve
p.num.size_tol_K     = 0.005;    % |peak T - T_max| to accept the area
p.num.size_maxit     = 60;       % max sizing iterations
p.num.A_bracket_m2   = [0.02 50];% search range for radiator area [m^2]

%% Studies
p.study.Tt_sweep_C = -30:2:50;   % transition temperatures to sweep [degC]

p.mc.N      = 500;               % Monte Carlo runs
p.mc.seed   = 1;                 % RNG seed (repeatable)
% Uniform ranges [lo hi]. All [ASSUME]; centred on the defaults above.
p.mc.w_K        = [2 8];         % STARC transition width
p.mc.eps_c      = [0.15 0.30];   % STARC cold emissivity
p.mc.eps_h      = [0.85 0.95];   % STARC hot emissivity (paper: 0.90-0.95)
p.mc.Qint_scale = [0.8 1.2];     % multiplier on Q_int (hot and cold drawn separately)
p.mc.dalpha     = [-0.02 0.08];  % change in solar absorptance (EOL darkening
                                 % mostly positive), same draw for both coatings

%% Output
p.out.show_figures = true;       % false = render off-screen, only save PNGs
end

function c = mk(name, alpha, eps_c, eps_h, Tt_C, w)
c = struct('name', name, 'alpha', alpha, 'eps_c', eps_c, 'eps_h', eps_h, ...
           'Tt_C', Tt_C, 'w', w);
end
