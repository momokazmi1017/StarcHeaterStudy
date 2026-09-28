"""One-node radiator thermal model, radiator sizing and fair coating comparison.

Energy balance of a single lumped radiator node (one radiating face, back
side adiabatic)::

    C dT/dt = Q_int + Q_solar + Q_albedo + Q_earthIR + Q_heater
              - eps(T) sigma A T^4

    Q_solar   = alpha G A F_sun             (sunlit only)
    Q_albedo  = alpha a G A F_earth         (sunlit only)
    Q_earthIR = eps(T) q_IR A F_earth       (whole orbit)

Kirchhoff's law: in the (grey) IR band absorptance = emittance, so the
Earth-IR term uses the same eps(T) as emission. A cold STARC surface
therefore also ABSORBS less Earth IR, which partly offsets its reduced
emission. Solar absorptance alpha is a separate property.

All cases (coatings, sweep points, Monte Carlo draws) are simulated together
as numpy vectors of length N.
"""

import numpy as np


# ---------------------------------------------------------------------------
def emissivity(T, eps_c, eps_h, Tt, w):
    """Logistic emissivity eps(T) = eps_c + (eps_h-eps_c)/(1+exp(-(T-Tt)/w)).

    T and Tt share units (K or degC); w is in K. Arrays broadcast. With
    eps_c == eps_h this is a fixed coating. This is an idealised curve; the
    real STARC response (shape, hysteresis, angle/wavelength) is not public.
    """
    with np.errstate(over="ignore"):
        return eps_c + (eps_h - eps_c) / (1.0 + np.exp(-(np.asarray(T) - Tt) / w))


def stack_coatings(p, names):
    """Combine named coatings from p['coat'] into one dict of numpy vectors."""
    coat = {f: np.array([p["coat"][n][f] for n in names], dtype=float)
            for f in ("alpha", "eps_c", "eps_h", "Tt_C", "w")}
    coat["name"] = [p["coat"][n]["name"] for n in names]
    coat["key"] = list(names)
    return coat


def config(p, coat, load, A):
    """Build the length-N case arrays used by ``simulate``.

    coat : dict with alpha, eps_c, eps_h, Tt_C, w (scalars or length-N arrays)
    load : dict with Qint_sun, Qint_ecl [W] (e.g. p['hot'] or p['cold'])
    A    : radiator area [m^2]

    The node's thermal capacitance follows the spec:
        C = (panel areal mass * A + equipment mass) * c_p
    """
    raw = {"A": A, "alpha": coat["alpha"], "eps_c": coat["eps_c"],
           "eps_h": coat["eps_h"], "Tt_C": coat["Tt_C"], "w": coat["w"],
           "Qsun": load["Qint_sun"], "Qecl": load["Qint_ecl"],
           "m_eq": p["sc"]["equip_mass"]}
    raw = {k: np.atleast_1d(np.asarray(v, dtype=float)).ravel() for k, v in raw.items()}
    sizes = {v.size for v in raw.values()}
    N = max(sizes)
    if not sizes <= {1, N}:
        raise ValueError("Inputs must be scalars or share one length N.")
    col = {k: np.broadcast_to(v, (N,)).copy() for k, v in raw.items()}

    cfg = {
        "N": N,
        "A": col["A"],
        "C": (p["sc"]["panel_areal_mass"] * col["A"] + col["m_eq"]) * p["sc"]["cp"],
        "alpha": col["alpha"],
        "eps_c": col["eps_c"],
        "eps_h": col["eps_h"],
        "Tt": col["Tt_C"] + p["T0K"],     # transition temperature [K]
        "w": col["w"],                   # transition width [K]
        "Qsun": col["Qsun"],             # internal dissipation, sunlit [W]
        "Qecl": col["Qecl"],             # internal dissipation, eclipse [W]
        "Tmin": np.full(N, p["sc"]["T_min_C"] + p["T0K"]),  # heater set point [K]
    }
    if np.any(cfg["A"] <= 0) or np.any(cfg["w"] <= 0):
        raise ValueError("Area and transition width must be positive.")
    return cfg


# ---------------------------------------------------------------------------
def _coefficients(cfg, p):
    """Temperature-independent parts of each heat term [W]."""
    env = p["env"]
    Qsol = cfg["alpha"] * cfg["A"] * env["G"] * env["F_sun"]                    # direct sun
    Qalb = cfg["alpha"] * cfg["A"] * env["albedo"] * env["G"] * env["F_earth"]  # albedo
    return {
        "QinSun": cfg["Qsun"] + Qsol + Qalb,  # Q_int + solar + albedo, sunlit
        "QinEcl": cfg["Qecl"],                # Q_int only, eclipse
        "IR": env["q_IR"] * env["F_earth"] * cfg["A"],  # x eps(T) -> Earth IR absorbed
        "rad": p["sigma"] * cfg["A"],                   # x eps(T) T^4 -> emitted
    }


def _flux(T, Q0, cfg, k):
    e = emissivity(T, cfg["eps_c"], cfg["eps_h"], cfg["Tt"], cfg["w"])
    qin = Q0 + e * k["IR"]         # Earth IR absorbed with absorptance eps(T)
    qout = e * k["rad"] * T**4     # emitted to space
    return qin, qout


def _one_orbit(T0, cfg, k, dt, n_step, t_sun, record):
    """Integrate one orbit from T0 with RK4 and the ideal thermostat."""
    C, Tmin = cfg["C"], cfg["Tmin"]
    T = T0.copy()
    zeros = np.zeros_like(T0)
    Ein, Eout, Eh, Eh_ecl, Qh_pk = (zeros.copy() for _ in range(5))
    T_pk, T_lo = T0.copy(), T0.copy()
    if record:
        Th = np.empty((T0.size, n_step + 1))
        Th[:, 0] = T0
        Qh = np.empty((T0.size, n_step))

    for s in range(n_step):
        sunlit = (s + 0.5) * dt < t_sun      # phase is constant within the step
        Q0 = k["QinSun"] if sunlit else k["QinEcl"]

        i1, o1 = _flux(T, Q0, cfg, k)
        i2, o2 = _flux(T + 0.5 * dt * (i1 - o1) / C, Q0, cfg, k)
        i3, o3 = _flux(T + 0.5 * dt * (i2 - o2) / C, Q0, cfg, k)
        i4, o4 = _flux(T + dt * (i3 - o3) / C, Q0, cfg, k)
        qin = (i1 + 2 * i2 + 2 * i3 + i4) / 6    # RK4-weighted heat in  [W]
        qout = (o1 + 2 * o2 + 2 * o3 + o4) / 6   # RK4-weighted heat out [W]

        # Same weights for T and for the energy tallies, so the bookkeeping
        # closes exactly: C*(T_end - T_start) = E_in - E_out.
        T = T + dt * (qin - qout) / C
        Ein += qin * dt
        Eout += qout * dt

        # Ideal thermostatic heater: add exactly the energy to hold T_min
        dEh = np.maximum(Tmin - T, 0.0) * C
        T = np.maximum(T, Tmin)
        Eh += dEh
        if not sunlit:
            Eh_ecl += dEh

        np.maximum(T_pk, T, out=T_pk)
        np.minimum(T_lo, T, out=T_lo)
        np.maximum(Qh_pk, dEh / dt, out=Qh_pk)
        if record:
            Th[:, s + 1] = T
            Qh[:, s] = dEh / dt

    out = {
        "Tpeak_K": T_pk, "Tlow_K": T_lo,
        "E_heater_J": Eh, "E_heater_ecl_J": Eh_ecl, "E_heater_sun_J": Eh - Eh_ecl,
        "Qh_peak_W": Qh_pk,
        "E_in_J": Ein + Eh,   # heater counts as energy in
        "E_out_J": Eout,
    }
    if record:
        out["T_K"] = Th
        out["Qh_W"] = Qh
    return T, out


def simulate(cfg, p, dt, record=True):
    """Periodic steady-state orbit of the one-node radiator model.

    Fixed-step RK4; the sun/eclipse phase is held constant within a step
    (dt must divide both the sunlit time and the period, so no step
    straddles a terminator). The start-of-orbit temperature T0 is found by
    a secant iteration on g(T0) = T(end of orbit) - T0; the converged orbit
    is then re-run and reported.

    Heater: ideal thermostat. After each step, if T fell below T_min the
    heater adds exactly C*(T_min - T). When the node sits at T_min this
    equals the instantaneous net heat loss.

    Returns a dict of length-N arrays: Tpeak_C, Tlow_C, E_heater_J (+ _ecl,
    _sun), Qh_peak_W, Qh_avg_W, E_in_J, E_out_J, closure, converged; with
    record=True also t_min, T_C / T_K (N x n+1) and tq_min, Qh_W (N x n).
    """
    P = p["orbit"]["period_min"] * 60.0
    t_sun = p["orbit"]["sun_min"] * 60.0
    n_step = int(round(P / dt))
    if abs(n_step * dt - P) > 1e-9 or abs(round(t_sun / dt) * dt - t_sun) > 1e-9:
        raise ValueError("dt must divide both the sunlit time and the orbit period.")

    k = _coefficients(cfg, p)
    Tmin = cfg["Tmin"]

    # Initial guess: orbit-average radiative equilibrium at mid emissivity
    fs = t_sun / P
    Qavg = fs * k["QinSun"] + (1 - fs) * k["QinEcl"]
    emid = 0.5 * (cfg["eps_c"] + cfg["eps_h"])
    Tg = ((Qavg + emid * k["IR"]) / (emid * k["rad"])) ** 0.25

    # Secant iteration on the start-of-orbit temperature
    tol = p["num"]["periodic_tol_K"]
    orbit = lambda T0: _one_orbit(T0, cfg, k, dt, n_step, t_sun, False)[0]
    Ta = np.maximum(Tg, Tmin)
    Pa = orbit(Ta)
    ga = Pa - Ta
    Tb = Pa
    Pb = orbit(Tb)
    gb = Pb - Tb
    it = 2
    while np.any(np.abs(gb) > tol) and it < p["num"]["periodic_maxit"]:
        denom = gb - ga
        with np.errstate(divide="ignore", invalid="ignore"):
            Tc = Tb - gb * (Tb - Ta) / denom
        bad = ~np.isfinite(Tc) | (np.abs(denom) < 1e-12)
        Tc[bad] = Pb[bad]                         # fall back to fixed-point step
        Tc = np.clip(Tc, Tmin, 2000.0)            # keep iterates physical
        Ta, ga = Tb, gb
        Tb = Tc
        Pb = orbit(Tb)
        gb = Pb - Tb
        it += 1

    _, sim = _one_orbit(Tb, cfg, k, dt, n_step, t_sun, record)
    sim["converged"] = np.abs(gb) <= tol
    sim["iterations"] = it
    sim["closure"] = (sim["E_in_J"] - sim["E_out_J"]) / sim["E_out_J"]
    sim["Tpeak_C"] = sim["Tpeak_K"] - p["T0K"]
    sim["Tlow_C"] = sim["Tlow_K"] - p["T0K"]
    sim["Qh_avg_W"] = sim["E_heater_J"] / P
    if record:
        sim["t_min"] = np.arange(n_step + 1) * dt / 60.0
        sim["tq_min"] = (np.arange(1, n_step + 1) - 0.5) * dt / 60.0
        sim["T_C"] = sim["T_K"] - p["T0K"]
    return sim


# ---------------------------------------------------------------------------
def size_radiator(p, coat, load, dt):
    """Radiator area whose hot-case orbit peaks at T_max.

    Peak temperature falls monotonically with area (more emitting surface,
    more thermal mass), so a bracketed Illinois (modified regula falsi)
    iteration on log(A) is used, all cases in parallel.

    Returns (A, info) with A a length-N array.
    """
    Tmax = p["sc"]["T_max_C"] + p["T0K"]
    N = config(p, coat, load, 1.0)["N"]

    def f(x):
        sim = simulate(config(p, coat, load, np.exp(x)), p, dt, record=False)
        return sim["Tpeak_K"] - Tmax

    lo, hi = p["num"]["A_bracket_m2"]
    a = np.full(N, np.log(lo)); fa = f(a)
    b = np.full(N, np.log(hi)); fb = f(b)
    if not (np.all(fa > 0) and np.all(fb < 0)):
        raise RuntimeError("T_max cannot be met inside A_bracket_m2 for every case "
                           "(too hot even at the largest area, or too cold at the smallest).")

    x, fx = b, fb
    for it in range(1, p["num"]["size_maxit"] + 1):
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)
        if np.all(np.abs(fx) < p["num"]["size_tol_K"]):
            break
        same_b = np.sign(fx) == np.sign(fb)
        fa = np.where(same_b, fa / 2, fb)        # Illinois down-weighting
        a = np.where(same_b, a, b)
        b, fb = x, fx

    info = {"iterations": it, "peak_error_K": fx,
            "converged": np.abs(fx) < p["num"]["size_tol_K"]}
    return np.exp(x), info


def compare(p, names, dt=None):
    """Fair comparison: size each coating on the hot case, run the cold case.

    1. Each coating gets its own radiator area so the hot case (sunlit,
       high Q_int) peaks exactly at T_max.
    2. The cold case (low Q_int) runs on that same area with the thermostat
       holding T >= T_min.
    """
    dt = p["num"]["dt"] if dt is None else dt
    coat = stack_coatings(p, names)
    A, sizing = size_radiator(p, coat, p["hot"], dt)
    return {
        "A": A,
        "sizing": sizing,
        "hot": simulate(config(p, coat, p["hot"], A), p, dt, record=True),
        "cold": simulate(config(p, coat, p["cold"], A), p, dt, record=True),
        "names": coat["name"],
        "keys": coat["key"],
        "dt": dt,
    }
