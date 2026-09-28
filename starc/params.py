"""Every model input in one place.

Call ``default_params()`` for a fresh, editable copy, e.g.::

    p = default_params()
    p["coat"]["STARC"]["Tt_C"] = 15

Source tags used in the comments:
    [WP-T3]  Aetara White Paper 01 (Sep 2026), Table 3 "Common assumptions"
    [WP-S1]  White paper, Customer Economics, Step 1 (surface effect)
    [WP-S2]  White paper, Customer Economics, Step 2 (value of a heater watt)
    [WP-A]   White paper, STARC properties quoted as company data [A]
    [STD]    Standard physical constant or textbook LEO value
    [ASSUME] Modelling choice made for this tool; not from the paper

FIRST-ORDER PUBLIC-DATA ESTIMATE. Not validated against STARC test data.
"""

import copy


def _coating(name, alpha, eps_c, eps_h, Tt_C, w):
    return {"name": name, "alpha": alpha, "eps_c": eps_c, "eps_h": eps_h,
            "Tt_C": Tt_C, "w": w}


_DEFAULTS = {
    # Physical constants
    "sigma": 5.670374419e-8,   # Stefan-Boltzmann [W/m^2/K^4]            [STD]
    "T0K": 273.15,             # 0 degC in kelvin                         [STD]

    # Orbit (reference LEO)
    "orbit": {
        "period_min": 95,      # orbit period [min]                       [WP-T3]
        "sun_min": 60,         # sunlit part [min]; eclipse = rest (35)   [WP-T3]
    },

    # Thermal environment
    "env": {
        "G": 1361.0,           # solar constant [W/m^2]                   [STD]
        "albedo": 0.30,        # Earth albedo [-]                         [STD]
        "q_IR": 237.0,         # Earth outgoing longwave IR [W/m^2]       [STD]
        "F_sun": 0.25,         # effective solar projection factor on the
                               # radiator face while sunlit [-]           [ASSUME]
        "F_earth": 0.30,       # radiator-to-Earth view factor [-]        [ASSUME]
    },

    # Spacecraft / radiator node
    "sc": {
        "panel_areal_mass": 5.0,  # radiator panel system mass [kg/m^2]   [WP-T3]
        "equip_mass": 10.0,       # equipment lumped with the node [kg]   [ASSUME]
        "cp": 900.0,              # specific heat, aluminium [J/kg/K]     [STD]
        "T_max_C": 40.0,          # hot-case limit used to size area      [ASSUME]
        "T_min_C": -10.0,         # cold-case heater set point            [ASSUME]
    },

    # Internal dissipation Q_int [W] per case and orbit phase.
    # Hot: payload on all orbit. Cold: payload off (bus only).
    # Set *_sun and *_ecl differently to model e.g. payload on only in sun.
    "hot": {"Qint_sun": 250.0, "Qint_ecl": 250.0},    # [ASSUME]
    "cold": {"Qint_sun": 50.0, "Qint_ecl": 50.0},     # [ASSUME]

    # Coatings. Every coating uses the same logistic model
    #   eps(T) = eps_c + (eps_h - eps_c) / (1 + exp(-(T - Tt)/w)),
    # so a fixed coating is just eps_c = eps_h (Tt_C and w then unused).
    "coat": {
        "FEP": _coating("Silver FEP", 0.10, 0.78, 0.78, 0.0, 1.0),  # eps 0.75-0.82 [WP-T3]
        "OSR": _coating("Glass OSR", 0.10, 0.88, 0.88, 0.0, 1.0),   # [STD]
        "STARC": _coating("STARC", 0.10, 0.20, 0.92, 10.0, 4.0),    # eps range [WP-A];
                                                                     # Tt, w [ASSUME]
    },
    "baseline": "OSR",         # fixed coating used as reference in sweep & MC

    # Economics (white paper Table 3 / Step 2)
    "econ": {
        "battery_Wh_per_kg": 150.0,   # usable battery specific energy   [WP-T3]
        "array_W_per_kg": 60.0,       # solar-array specific power       [WP-T3]
        "recharge_eff": 0.85,         # array -> battery -> load         [WP-T3]
        "launch_usd_per_kg": 7000.0,  # rideshare price                  [WP-T3, E1]
        "usd_per_W": (150.0, 800.0),  # heater-watt value, constellation
                                      # .. science                       [WP-S2]
        "ft2_per_m2": 10.7639,        # unit conversion                  [STD]
    },

    # Numerics
    "num": {
        "dt": 1.0,               # RK4 step for headline results [s]
        "dt_study": 5.0,         # RK4 step for sweep + Monte Carlo [s]
                                 # (run_starc_study checks it against dt)
        "periodic_tol_K": 1e-4,  # |T(end of orbit) - T(start)| to accept
        "periodic_maxit": 60,    # max orbits in the periodic solve
        "size_tol_K": 0.005,     # |peak T - T_max| to accept the area
        "size_maxit": 60,        # max sizing iterations
        "A_bracket_m2": (0.02, 50.0),  # search range for radiator area [m^2]
    },

    # Studies
    "study": {"Tt_sweep_C": list(range(-30, 51, 2))},  # [degC]
    "mc": {
        "N": 500,
        "seed": 1,
        # Uniform ranges (lo, hi). All [ASSUME]; centred on the defaults.
        "w_K": (2.0, 8.0),           # STARC transition width
        "eps_c": (0.15, 0.30),       # STARC cold emissivity
        "eps_h": (0.85, 0.95),       # STARC hot emissivity (paper: 0.90-0.95)
        "Qint_scale": (0.8, 1.2),    # multiplier on Q_int (hot, cold drawn separately)
        "dalpha": (-0.02, 0.08),     # change in solar absorptance (EOL darkening
                                     # mostly positive); same draw for both coatings
    },
}


def default_params():
    """Return a fresh deep copy of the default inputs."""
    return copy.deepcopy(_DEFAULTS)
