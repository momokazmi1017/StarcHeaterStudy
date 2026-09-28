"""Economics and white-paper hand-calculation helpers."""

import numpy as np


def radiated_flux(eps, T_C, p):
    """Grey-body emission to deep space, eps*sigma*T^4 [W/m^2] (paper Step 1)."""
    return eps * p["sigma"] * (np.asarray(T_C) + p["T0K"]) ** 4


def datacenter_area(P_compute_W, overhead, T_K, eps, q_env, n_faces, sigma):
    """First-order radiator panel area for an orbital data center [m^2].

    Heat to reject = P_compute * (1 + overhead). Each face radiates
    eps*sigma*T^4 and absorbs q_env [W/m^2] from the environment:

        A = P (1 + overhead) / (n_faces * (eps sigma T^4 - q_env))

    White paper: 1 MW, 10 %, 330 K, eps 0.90, 150 W/m^2, two faces -> ~1,200 m^2.
    """
    return P_compute_W * (1 + overhead) / (n_faces * (eps * sigma * T_K**4 - q_env))


def economics(E_ecl_J, E_sun_J, p):
    """Battery + array mass and dollar value of heater energy per orbit.

    E_ecl_J : heater energy drawn during eclipse, per orbit [J]
    E_sun_J : heater energy drawn in sunlight, per orbit [J]
    (arrays allowed; pass a SAVING to get the value of that saving)

    Follows white paper Step 2 / Table 3:
      battery = eclipse energy / usable specific energy
      array   = power needed over the sunlit period to (a) recharge the
                eclipse energy through the recharge efficiency and (b) run
                the sunlit heater directly, / array specific power
    e.g. 1 W through a 35 min eclipse -> 3.9 g battery + 11.4 g array.

    A "heater watt" in the paper is 1 W held through the whole eclipse.
    W_equiv converts any heater profile into that unit by mass, and the
    paper's $150-$800 per watt is applied to W_equiv.
    """
    ec = p["econ"]
    t_sun_h = p["orbit"]["sun_min"] / 60
    t_ecl_h = (p["orbit"]["period_min"] - p["orbit"]["sun_min"]) / 60

    E_ecl_Wh = np.asarray(E_ecl_J, dtype=float) / 3600
    E_sun_Wh = np.asarray(E_sun_J, dtype=float) / 3600

    e = {}
    e["battery_kg"] = E_ecl_Wh / ec["battery_Wh_per_kg"]
    e["array_W"] = (E_ecl_Wh / ec["recharge_eff"] + E_sun_Wh) / t_sun_h
    e["array_kg"] = e["array_W"] / ec["array_W_per_kg"]
    e["mass_kg"] = e["battery_kg"] + e["array_kg"]
    e["battery_g"] = 1e3 * e["battery_kg"]
    e["array_g"] = 1e3 * e["array_kg"]
    e["mass_g"] = 1e3 * e["mass_kg"]

    # Mass for one eclipse heater watt (white paper: ~15.3 g)
    e["kg_per_eclipse_W"] = (t_ecl_h / ec["battery_Wh_per_kg"]
                             + (t_ecl_h / ec["recharge_eff"]) / t_sun_h / ec["array_W_per_kg"])
    e["W_equiv"] = e["mass_kg"] / e["kg_per_eclipse_W"]

    e["usd_launch"] = e["mass_kg"] * ec["launch_usd_per_kg"]  # launch cost of the mass only
    e["usd_low"] = e["W_equiv"] * ec["usd_per_W"][0]          # constellation value
    e["usd_high"] = e["W_equiv"] * ec["usd_per_W"][-1]        # science-mission value
    return e
