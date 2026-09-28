"""First-order heater-power study for an adaptive-emissivity radiator.

Run:   python run_starc_study.py

Compares a STARC-like adaptive-emissivity coating with fixed silver FEP and
glass OSR on a one-node LEO radiator:
  1. sizes each coating's radiator so the hot case peaks at T_max,
  2. runs the cold case on that area with a thermostatic heater,
  3. converts heater energy saved into battery + array mass and $,
  4. sweeps the transition temperature T_t,
  5. runs a 500-case Monte Carlo over uncertain inputs.
Plots go to results/*.png and key numbers to results/summary.txt.
All inputs live in starc/params.py.

FIRST-ORDER PUBLIC-DATA ESTIMATE. NOT VALIDATED AGAINST STARC TEST DATA.
"""

import time
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import starc

OUT = Path(__file__).resolve().parent / "results"
DEG = "\N{DEGREE SIGN}C"

# Fixed colour per coating (validated categorical slots 1-3) plus a distinct
# line style, so identity never relies on colour alone.
COLOR = {"OSR": "#2a78d6", "FEP": "#eb6834", "STARC": "#1baf7a"}
LINE = {"OSR": "--", "FEP": "-.", "STARC": "-"}
INK = "#52514e"          # secondary text / reference lines
SHADE = "#8c8c8c"

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#6b6b6b", "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e3e3e3", "grid.linewidth": 0.8,
    "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9,
    "legend.frameon": False, "font.size": 10,
})


# ---------------------------------------------------------------------------
def coating_savings(R, iS, p):
    """STARC savings relative to every coating in R (row iS is STARC)."""
    ft2 = p["econ"]["ft2_per_m2"]
    c = R["cold"]
    out = []
    for j in range(len(R["A"])):
        dEecl = c["E_heater_ecl_J"][j] - c["E_heater_ecl_J"][iS]
        dEsun = c["E_heater_sun_J"][j] - c["E_heater_sun_J"][iS]
        e = starc.economics(dEecl, dEsun, p)
        out.append({
            "E_saved_Wh": (dEecl + dEsun) / 3600,
            "Qavg_saved_W": c["Qh_avg_W"][j] - c["Qh_avg_W"][iS],
            "Qpeak_saved_W": c["Qh_peak_W"][j] - c["Qh_peak_W"][iS],
            "econ": e,
            "area_saved_pct": 100 * (R["A"][j] - R["A"][iS]) / R["A"][j],
            "usd_per_ft2": (e["usd_low"] / (R["A"][iS] * ft2), e["usd_high"] / (R["A"][iS] * ft2)),
            "Wequiv_per_ft2": e["W_equiv"] / (R["A"][iS] * ft2),
        })
    return out


def run_sweep(p, Rs, iFix, iB):
    """Re-size and re-run STARC for each transition temperature."""
    dt = p["num"]["dt_study"]
    Tt = np.asarray(p["study"]["Tt_sweep_C"], dtype=float)
    cS = dict(p["coat"]["STARC"])
    cS["Tt_C"] = Tt
    A, _ = starc.size_radiator(p, cS, p["hot"], dt)
    cold = starc.simulate(starc.config(p, cS, p["cold"], A), p, dt, record=False)
    # Heater energy saved vs each fixed coating [Wh/orbit]; columns follow iFix
    saved = (Rs["cold"]["E_heater_J"][iFix][None, :] - cold["E_heater_J"][:, None]) / 3600

    # Best switch point. "Good" T_t values save within 0.5 % (of the baseline
    # heater energy) of the maximum AND need a radiator within 0.1 % of the
    # smallest such area. Often a window, not a point: pick its centre so the
    # choice has the most margin against T_t uncertainty.
    sb = saved[:, iFix.index(iB)]
    tolW = 0.005 * Rs["cold"]["E_heater_J"][iB] / 3600
    cand = np.flatnonzero(sb >= sb.max() - tolW)
    good = cand[A[cand] <= 1.001 * A[cand].min()]
    window = (Tt[good].min(), Tt[good].max())
    best = good[np.argmin(np.abs(Tt[good] - np.mean(window)))]
    return {"Tt_C": Tt, "A": A, "cold": cold, "saved_Wh": saved,
            "window_C": window, "iBest": int(best)}


def run_monte_carlo(p):
    """Paired Monte Carlo: each draw is applied to STARC and to the baseline
    coating (same loads and alpha change); both are re-sized, then compared."""
    rng = np.random.default_rng(p["mc"]["seed"])
    N = p["mc"]["N"]
    U = lambda r: rng.uniform(r[0], r[1], N)
    d = {"w": U(p["mc"]["w_K"]), "eps_c": U(p["mc"]["eps_c"]), "eps_h": U(p["mc"]["eps_h"]),
         "kHot": U(p["mc"]["Qint_scale"]), "kCold": U(p["mc"]["Qint_scale"]),
         "dalpha": U(p["mc"]["dalpha"])}

    hot = {k: p["hot"][k] * d["kHot"] for k in ("Qint_sun", "Qint_ecl")}
    cold = {k: p["cold"][k] * d["kCold"] for k in ("Qint_sun", "Qint_ecl")}
    cS = dict(p["coat"]["STARC"], w=d["w"], eps_c=d["eps_c"], eps_h=d["eps_h"],
              alpha=p["coat"]["STARC"]["alpha"] + d["dalpha"])
    cB = dict(p["coat"][p["baseline"]])
    cB["alpha"] = cB["alpha"] + d["dalpha"]

    dt = p["num"]["dt_study"]
    A_S, _ = starc.size_radiator(p, cS, hot, dt)
    A_B, _ = starc.size_radiator(p, cB, hot, dt)
    colS = starc.simulate(starc.config(p, cS, cold, A_S), p, dt, record=False)
    colB = starc.simulate(starc.config(p, cB, cold, A_B), p, dt, record=False)

    e = starc.economics(colB["E_heater_ecl_J"] - colS["E_heater_ecl_J"],
                        colB["E_heater_sun_J"] - colS["E_heater_sun_J"], p)
    area = 100 * (A_B - A_S) / A_B
    pc = [10, 50, 90]
    return {"A_S": A_S, "A_B": A_B, "econ": e, "W_saved": e["W_equiv"],
            "area_saved_pct": area, "draws": d,
            "pW": np.percentile(e["W_equiv"], pc), "pA": np.percentile(area, pc),
            "pUSD": np.percentile(e["usd_low"], pc), "pUSDhi": np.percentile(e["usd_high"], pc),
            "converged": bool(colS["converged"].all() and colB["converged"].all())}


# ---------------------------------------------------------------------------
def shade_eclipse(ax, p):
    x1, x2 = p["orbit"]["sun_min"], p["orbit"]["period_min"]
    ax.axvspan(x1, x2, color=SHADE, alpha=0.12, lw=0, zorder=0)
    ax.text((x1 + x2) / 2, 0.98, "eclipse", transform=ax.get_xaxis_transform(),
            ha="center", va="top", color=INK, fontsize=9)
    ax.set_xlim(0, p["orbit"]["period_min"])


def ref_line(ax, y, label, horizontal=True, va="bottom", right=False, **kw):
    kw = {"color": INK, "ls": ":", "lw": 1, **kw}
    if horizontal:
        ax.axhline(y, **kw)
        ax.text(0.995 if right else 0.005, y, f" {label} ", transform=ax.get_yaxis_transform(),
                ha="right" if right else "left", va=va, color=kw["color"], fontsize=9)
    else:
        ax.axvline(y, **kw)
        ax.text(y, 0.02, " " + label, transform=ax.get_xaxis_transform(),
                ha="left", va="bottom", color=kw["color"], fontsize=9)


def plot_temperature(R, p, keys, file):
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=True, constrained_layout=True)
    cases = [
        (R["hot"], p["sc"]["T_max_C"], "$T_{max}$",
         f"Hot case: $Q_{{int}}$ = {p['hot']['Qint_sun']:g} W, radiator sized so peak = "
         f"{p['sc']['T_max_C']:g} {DEG}"),
        (R["cold"], p["sc"]["T_min_C"], "$T_{min}$",
         f"Cold case: $Q_{{int}}$ = {p['cold']['Qint_sun']:g} W on the same area, "
         f"heater holds T $\\geq$ {p['sc']['T_min_C']:g} {DEG}"),
    ]
    for ax, (sim, lim, lbl, title) in zip(axes, cases):
        for j, k in enumerate(keys):
            ax.plot(sim["t_min"], sim["T_C"][j], color=COLOR[k], ls=LINE[k], lw=2,
                    label=f"{R['names'][j]}  (A = {R['A'][j]:.2f} m$^2$)")
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(lo, lim - 3), max(hi, lim + 3))
        ref_line(ax, lim, lbl)
        shade_eclipse(ax, p)
        ax.set_title(title, loc="left")
        ax.set_ylabel(f"Radiator temperature ({DEG})")
    axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    axes[-1].set_xlabel("Time in orbit (min)")
    fig.savefig(file, dpi=150)
    plt.close(fig)


def plot_heater(R, p, keys, file):
    fig, ax = plt.subplots(figsize=(10, 4.4), constrained_layout=True)
    c = R["cold"]
    for j, k in enumerate(keys):
        ax.plot(c["tq_min"], c["Qh_W"][j], color=COLOR[k], ls=LINE[k], lw=2,
                label=f"{R['names'][j]}: {c['E_heater_J'][j] / 3600:.1f} Wh/orbit, "
                      f"peak {c['Qh_peak_W'][j]:.0f} W")
    ax.set_ylim(0, max(ax.get_ylim()[1], 1))
    shade_eclipse(ax, p)
    ax.set_title(f"Heater power, cold case ($Q_{{int}}$ = {p['cold']['Qint_sun']:g} W, "
                 f"$T_{{min}}$ = {p['sc']['T_min_C']:g} {DEG})", loc="left")
    ax.set_xlabel("Time in orbit (min)")
    ax.set_ylabel("Heater power (W)")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.savefig(file, dpi=150)
    plt.close(fig)


def plot_emissivity(p, keys, file):
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    T = np.arange(-60, 80.01, 0.25)
    for k in keys:
        c = p["coat"][k]
        e = starc.emissivity(T, c["eps_c"], c["eps_h"], c["Tt_C"], c["w"]) * np.ones_like(T)
        lbl = c["name"]
        if c["eps_c"] != c["eps_h"]:
            lbl += f"  ($T_t$ = {c['Tt_C']:g} {DEG}, w = {c['w']:g} K)"
        ax.plot(T, e, color=COLOR[k], ls=LINE[k], lw=2, label=lbl)
    ref_line(ax, p["sc"]["T_min_C"], "$T_{min}$", horizontal=False)
    ref_line(ax, p["sc"]["T_max_C"], "$T_{max}$", horizontal=False)
    ax.set_ylim(0, 1)
    ax.set_xlim(T[0], T[-1])
    ax.set_title("Infrared emissivity vs temperature (idealised logistic STARC curve)", loc="left")
    ax.set_xlabel(f"Surface temperature ({DEG})")
    ax.set_ylabel("Emissivity $\\epsilon$ (-)")
    ax.legend(loc="center left")
    fig.savefig(file, dpi=150)
    plt.close(fig)


def plot_sweep(sw, Rs, p, keys, iFix, iB, iS, file):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 7.2), sharex=True, constrained_layout=True)
    b = sw["iBest"]
    Tt = sw["Tt_C"]
    ax1.axvspan(*sw["window_C"], color=COLOR["STARC"], alpha=0.12, lw=0, zorder=0,
                label="good $T_t$ window")
    for col, j in enumerate(iFix):
        ax1.plot(Tt, sw["saved_Wh"][:, col], color=COLOR[keys[j]], ls=LINE[keys[j]], lw=2,
                 label=f"vs {Rs['names'][j]}")
    kB = iFix.index(iB)
    ax1.plot(Tt[b], sw["saved_Wh"][b, kB], "o", ms=9, mew=2, mec="white",
             mfc=COLOR["STARC"], label=f"best $T_t$ = {Tt[b]:g} {DEG}", zorder=5)
    ax1.axhline(0, color=INK, lw=1)
    ax1.set_title("Cold-case heater energy saved by STARC (each coating sized for its own hot case)",
                  loc="left")
    ax1.set_ylabel("Heater energy saved (Wh/orbit)")
    ax1.legend(loc="center right")

    ax2.plot(Tt, sw["A"], color=COLOR["STARC"], ls=LINE["STARC"], lw=2, label="STARC")
    top = max(iFix, key=lambda j: Rs["A"][j])
    for j in iFix:
        ref_line(ax2, Rs["A"][j], Rs["names"][j], va="bottom" if j == top else "top", right=True,
                 color=COLOR[keys[j]], ls=LINE[keys[j]], lw=1.5)
    ax2.plot(Tt[b], sw["A"][b], "o", ms=9, mew=2, mec="white", mfc=COLOR["STARC"], zorder=5)
    for ax in (ax1, ax2):
        ax.axvline(Tt[b], color=INK, ls=":", lw=1)
    ax2.set_title("Radiator area needed to hold the hot-case peak at $T_{max}$", loc="left")
    ax2.set_ylabel("Radiator area (m$^2$)")
    ax2.set_xlabel(f"STARC transition temperature $T_t$ ({DEG})")
    ax2.set_xlim(Tt.min(), Tt.max())
    ax2.legend(loc="upper left")
    fig.savefig(file, dpi=150)
    plt.close(fig)


def plot_mc(mc, p, file):
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    ax.hist(mc["W_saved"], bins=30, color=COLOR["STARC"], alpha=0.9, edgecolor="white", lw=1)
    ax.grid(axis="x", visible=False)
    for q, lbl in zip(mc["pW"], ("P10", "P50", "P90")):
        ax.axvline(q, color=INK, ls="--", lw=1.2)
        ax.text(q, 0.98, f" {lbl} = {q:.1f} W", transform=ax.get_xaxis_transform(),
                rotation=90, ha="right", va="top", color=INK, fontsize=9)
    ax.set_title(f"Monte Carlo ({p['mc']['N']} runs): heater watts saved by STARC vs "
                 f"{p['coat'][p['baseline']]['name']}", loc="left")
    ax.set_xlabel("Heater power saved (equivalent eclipse watts, W)")
    ax.set_ylabel("Runs")
    fig.savefig(file, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
def money(x):
    x = float(x)
    x = round(x, -1) if abs(x) >= 100 else round(x)
    return f"{x:,.0f}"


def build_summary(p, R, sav, iS, iB, iFix, sw, mc, dt_check, runtime):
    ft2 = p["econ"]["ft2_per_m2"]
    d = "C"
    c = R["cold"]
    s = p["coat"]["STARC"]
    L = [
        "=" * 72,
        " STARC adaptive-emissivity radiator: first-order heater-power study",
        "=" * 72,
        f" Generated {datetime.now():%Y-%m-%d %H:%M}",
        "",
        " *** FIRST-ORDER PUBLIC-DATA ESTIMATE. NOT VALIDATED AGAINST STARC ***",
        " *** TEST DATA. One-node model; see README.md for all limitations. ***",
        "",
        " SCENARIO",
        f"   LEO orbit {p['orbit']['period_min']:g} min ({p['orbit']['sun_min']:g} min sun, "
        f"{p['orbit']['period_min'] - p['orbit']['sun_min']:g} min eclipse); "
        f"F_sun {p['env']['F_sun']:.2f}, F_earth {p['env']['F_earth']:.2f}",
        f"   Hot case  Q_int {p['hot']['Qint_sun']:g} W (sun) / {p['hot']['Qint_ecl']:g} W (eclipse); "
        f"radiator sized so peak = {p['sc']['T_max_C']:g} {d}",
        f"   Cold case Q_int {p['cold']['Qint_sun']:g} W (sun) / {p['cold']['Qint_ecl']:g} W (eclipse); "
        f"heater holds T >= {p['sc']['T_min_C']:g} {d}",
        f"   STARC eps {s['eps_c']:.2f} -> {s['eps_h']:.2f}, T_t = {s['Tt_C']:g} {d}, "
        f"w = {s['w']:g} K, alpha {s['alpha']:.2f}",
        "",
        f" 1. RADIATOR AREA AND COLD-CASE HEATER DEMAND (dt = {p['num']['dt']:g} s)",
        f"   {'Coating':<11} {'Area m2':>8} {'Area ft2':>8} {'Heater Wh/o':>12} {'eclipse Wh':>11} "
        f"{'Peak W':>8} {'Avg W':>7}  Cold T range",
    ]
    for j, name in enumerate(R["names"]):
        L.append(f"   {name:<11} {R['A'][j]:8.3f} {R['A'][j] * ft2:8.2f} "
                 f"{c['E_heater_J'][j] / 3600:12.2f} {c['E_heater_ecl_J'][j] / 3600:11.2f} "
                 f"{c['Qh_peak_W'][j]:8.1f} {c['Qh_avg_W'][j]:7.1f}  "
                 f"{c['Tlow_C'][j]:.1f} to {c['Tpeak_C'][j]:.1f} {d}")
    L += ["", " 2. WHAT STARC SAVES"]
    for j in iFix:
        sv, e = sav[j], sav[j]["econ"]
        L += [
            f"   vs {R['names'][j]}{'  (baseline)' if j == iB else ''}",
            f"     Heater energy saved          {sv['E_saved_Wh']:8.1f} Wh per orbit",
            f"     Heater power saved           {sv['Qavg_saved_W']:8.1f} W orbit-average, "
            f"{sv['Qpeak_saved_W']:.1f} W peak",
            f"     = equivalent eclipse watts   {float(e['W_equiv']):8.1f} W  (paper Step 2 unit)",
            f"     Battery mass saved           {float(e['battery_g']):8.0f} g",
            f"     Solar-array mass saved       {float(e['array_g']):8.0f} g",
            f"     Total mass saved             {float(e['mass_g']):8.0f} g  "
            f"(launch only: ${money(e['usd_launch'])} at ${money(p['econ']['launch_usd_per_kg'])}/kg)",
            f"     Value at ${p['econ']['usd_per_W'][0]:g}-${p['econ']['usd_per_W'][-1]:g} per W"
            f"      ${money(e['usd_low'])} - ${money(e['usd_high'])} per spacecraft",
            f"     Per sq ft of STARC radiator  ${money(sv['usd_per_ft2'][0])} - "
            f"${money(sv['usd_per_ft2'][1])}  ({float(sv['Wequiv_per_ft2']):.1f} W/sq ft; "
            f"paper credits 3 W/sq ft)",
            f"     Radiator area saved          {sv['area_saved_pct']:8.1f} %",
        ]
    b, kB = sw["iBest"], iFix.index(iB)
    L += [
        "",
        " 3. TRANSITION-TEMPERATURE SWEEP",
        f"   Best T_t for this mission: {sw['Tt_C'][b]:g} {d}  (saves {sw['saved_Wh'][b, kB]:.1f} "
        f"Wh/orbit vs {R['names'][iB]}, area {sw['A'][b]:.3f} m2).",
        f"   Good window: T_t = {sw['window_C'][0]:g} to {sw['window_C'][1]:g} {d} (heater saving "
        f"within 0.5% of max and",
        "   radiator within 0.1% of the smallest); best = centre of that window.",
        "   Too low a T_t keeps eps high when cold (heater still needed); too high",
        "   a T_t leaves eps low when hot, so the radiator must grow.",
        "",
        f" 4. MONTE CARLO ({p['mc']['N']} paired runs vs {R['names'][iB]})",
        f"   Uncertain: w {p['mc']['w_K'][0]:g}-{p['mc']['w_K'][1]:g} K, eps_c "
        f"{p['mc']['eps_c'][0]:.2f}-{p['mc']['eps_c'][1]:.2f}, eps_h {p['mc']['eps_h'][0]:.2f}-"
        f"{p['mc']['eps_h'][1]:.2f}, Q_int x{p['mc']['Qint_scale'][0]:.1f}-{p['mc']['Qint_scale'][1]:.1f},",
        f"   alpha change {p['mc']['dalpha'][0]:+.2f} to {p['mc']['dalpha'][1]:+.2f} "
        f"(uniform distributions).",
        "   Heater watts saved (eclipse-W)   P10 {:.1f}   P50 {:.1f}   P90 {:.1f} W".format(*mc["pW"]),
        f"   Value at ${p['econ']['usd_per_W'][0]:g}/W                  P10 ${money(mc['pUSD'][0])}  "
        f"P50 ${money(mc['pUSD'][1])}  P90 ${money(mc['pUSD'][2])}",
        f"   Value at ${p['econ']['usd_per_W'][-1]:g}/W                  P10 ${money(mc['pUSDhi'][0])}  "
        f"P50 ${money(mc['pUSDhi'][1])}  P90 ${money(mc['pUSDhi'][2])}",
        "   Radiator area saved              P10 {:.1f}%  P50 {:.1f}%  P90 {:.1f}%".format(*mc["pA"]),
        "",
        " 5. NUMERICAL CHECKS",
        f"   Periodic energy closure, last orbit: max |in-out|/out = "
        f"{np.abs(R['hot']['closure']).max():.1e} (hot), {np.abs(R['cold']['closure']).max():.1e} (cold)",
        f"   dt = {p['num']['dt_study']:g} s vs {p['num']['dt']:g} s: area differs by "
        f"{100 * dt_check['A_rel']:.3f}%, heater energy by {100 * dt_check['E_rel']:.3f}% of max",
        f"   Periodic solves converged: "
        f"{'yes' if R['hot']['converged'].all() and R['cold']['converged'].all() else 'NO'};  "
        f"sizing converged: {'yes' if R['sizing']['converged'].all() else 'NO'};  "
        f"Monte Carlo: {'yes' if mc['converged'] else 'NO'}",
        f"   Run time {runtime:.0f} s",
        "",
        " Heater watts are expressed in the white paper's unit: the battery + array",
        " mass of 1 W held through the 35 min eclipse (~15.3 g). $ values apply the",
        " paper's $150 (constellation) to $800 (science) per watt; they are gross",
        " savings, before the price difference of STARC vs the incumbent coating.",
    ]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    p = starc.default_params()
    keys = ["OSR", "FEP", "STARC"]          # fixed order = fixed plot colours
    iS = keys.index("STARC")
    iB = keys.index(p["baseline"])
    iFix = [j for j, k in enumerate(keys) if k != "STARC"]

    print("\nSTARC heater-power study\n"
          "  First-order public-data estimate - NOT validated against STARC data.\n")

    print(f"[1/5] Sizing radiators and running hot/cold cases (dt = {p['num']['dt']:g} s) ...")
    R = starc.compare(p, keys, p["num"]["dt"])
    sav = coating_savings(R, iS, p)
    # Check that the coarser step used for the sweep and Monte Carlo agrees
    Rs = starc.compare(p, keys, p["num"]["dt_study"])
    Escale = max(R["cold"]["E_heater_J"].max(), 1.0)
    dt_check = {"A_rel": np.max(np.abs(Rs["A"] - R["A"]) / R["A"]),
                "E_rel": np.max(np.abs(Rs["cold"]["E_heater_J"] - R["cold"]["E_heater_J"])) / Escale}

    print("[2/5] Plotting orbit histories and eps(T) ...")
    plot_temperature(R, p, keys, OUT / "fig1_temperature_orbit.png")
    plot_heater(R, p, keys, OUT / "fig2_heater_power_cold.png")
    plot_emissivity(p, keys, OUT / "fig3_starc_emissivity.png")

    print(f"[3/5] Sweeping transition temperature T_t (dt = {p['num']['dt_study']:g} s) ...")
    sw = run_sweep(p, Rs, iFix, iB)
    plot_sweep(sw, Rs, p, keys, iFix, iB, iS, OUT / "fig4_transition_temperature_sweep.png")

    print(f"[4/5] Monte Carlo, {p['mc']['N']} runs vs {p['baseline']} "
          f"(dt = {p['num']['dt_study']:g} s) ...")
    mc = run_monte_carlo(p)
    plot_mc(mc, p, OUT / "fig5_monte_carlo_heater_watts_saved.png")

    print("[5/5] Writing summary ...\n")
    txt = build_summary(p, R, sav, iS, iB, iFix, sw, mc, dt_check, time.time() - t0)
    (OUT / "summary.txt").write_text(txt, encoding="utf-8")
    print(txt)
    print(f"Results written to {OUT}")


if __name__ == "__main__":
    main()
