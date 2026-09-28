"""Validation suite: checks the model against the white paper and itself.

Run:  python -m pytest -v

Hand-calc checks reproduce numbers printed in Aetara White Paper 01; the
rest check the solver (limit case, energy closure, time-step convergence).
None of them validate the model against STARC hardware data, which is not
public.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import starc  # noqa: E402

NAMES = ["FEP", "OSR", "STARC"]


@pytest.fixture(scope="module")
def p():
    return starc.default_params()


@pytest.fixture(scope="module")
def R1(p):
    return starc.compare(p, NAMES, 1.0)


def test_step1_radiated_flux(p):
    # Paper Step 1: at -10 C, eps 0.80 radiates ~217 W/m^2 and eps 0.20 ~54 W/m^2
    assert starc.radiated_flux(0.80, -10, p) == pytest.approx(217, abs=2)
    assert starc.radiated_flux(0.20, -10, p) == pytest.approx(54, abs=1)


def test_heater_watt_cost(p):
    # Paper Step 2: 1 W through a 35 min eclipse -> ~3.9 g battery + ~11.4 g array = ~15.3 g
    e = starc.economics(1 * 35 * 60, 0, p)
    assert e["battery_g"] == pytest.approx(3.9, abs=0.05)
    assert e["array_g"] == pytest.approx(11.4, abs=0.1)
    assert e["mass_g"] == pytest.approx(15.3, abs=0.1)
    assert e["W_equiv"] == pytest.approx(1.0, rel=1e-12)


def test_orbital_datacenter_area(p):
    # Paper: 330 K, eps 0.90, two faces, 150 W/m^2 per face, 1 MW + 10 % -> ~1,200 m^2
    A = starc.datacenter_area(1e6, 0.10, 330, 0.90, 150, 2, p["sigma"])
    assert A == pytest.approx(1200, rel=0.02)


def test_logistic_curve_shape(p):
    c = p["coat"]["STARC"]
    f = lambda T: starc.emissivity(T, c["eps_c"], c["eps_h"], c["Tt_C"], c["w"])
    assert f(c["Tt_C"]) == pytest.approx((c["eps_c"] + c["eps_h"]) / 2, abs=1e-12)
    assert f(c["Tt_C"] - 20 * c["w"]) == pytest.approx(c["eps_c"], abs=1e-6)
    assert f(c["Tt_C"] + 20 * c["w"]) == pytest.approx(c["eps_h"], abs=1e-6)
    assert np.all(np.diff(f(np.arange(-50, 51))) > 0)


def test_equal_emissivities_give_fixed_coating(p):
    # STARC with eps_c = eps_h = 0.88 (same alpha) must behave exactly like glass OSR
    osr = p["coat"]["OSR"]
    s = dict(p["coat"]["STARC"], eps_c=osr["eps_c"], eps_h=osr["eps_h"], alpha=osr["alpha"])
    coat = {f: np.array([osr[f], s[f]], dtype=float) for f in ("alpha", "eps_c", "eps_h", "Tt_C", "w")}
    dt = p["num"]["dt_study"]
    A, _ = starc.size_radiator(p, coat, p["hot"], dt)
    cold = starc.simulate(starc.config(p, coat, p["cold"], A), p, dt, record=True)
    assert A[1] == pytest.approx(A[0], rel=1e-9)
    assert cold["E_heater_J"][1] == pytest.approx(cold["E_heater_J"][0], rel=1e-9)
    np.testing.assert_allclose(cold["T_C"][1], cold["T_C"][0], atol=1e-9)
    e = starc.emissivity(np.arange(-50, 81), s["eps_c"], s["eps_h"], s["Tt_C"], s["w"])
    np.testing.assert_allclose(e, 0.88, atol=1e-15)


def test_periodic_energy_closure(p, R1):
    # Over the reported (last) orbit, energy in = energy out within 1 %
    assert R1["hot"]["converged"].all() and R1["cold"]["converged"].all()
    assert np.abs(R1["hot"]["closure"]).max() < 0.01
    assert np.abs(R1["cold"]["closure"]).max() < 0.01
    # Independent check: re-integrate emitted power from the recorded
    # temperature history with the trapezoid rule
    cfg = starc.config(p, starc.stack_coatings(p, R1["keys"]), p["cold"], R1["A"])
    TK = R1["cold"]["T_K"]
    e = starc.emissivity(TK, cfg["eps_c"][:, None], cfg["eps_h"][:, None],
                         cfg["Tt"][:, None], cfg["w"][:, None])
    Eout = np.trapezoid(e * p["sigma"] * cfg["A"][:, None] * TK**4, R1["cold"]["t_min"] * 60, axis=1)
    np.testing.assert_allclose(Eout, R1["cold"]["E_out_J"], rtol=1e-3)
    # Sized hot case really peaks at T_max
    np.testing.assert_allclose(R1["hot"]["Tpeak_C"], p["sc"]["T_max_C"], atol=0.01)


def test_timestep_convergence(p, R1):
    # Halving dt (1 s -> 0.5 s) changes area, heater energy and savings by < 0.5 %
    R2 = starc.compare(p, NAMES, 0.5)
    assert np.max(np.abs(R2["A"] - R1["A"]) / R1["A"]) < 0.005
    Escale = max(R1["cold"]["E_heater_J"].max(), 1.0)
    assert np.max(np.abs(R2["cold"]["E_heater_J"] - R1["cold"]["E_heater_J"])) / Escale < 0.005
    s1 = R1["cold"]["E_heater_J"][:2] - R1["cold"]["E_heater_J"][2]
    s2 = R2["cold"]["E_heater_J"][:2] - R2["cold"]["E_heater_J"][2]
    assert np.max(np.abs(s2 - s1)) / max(np.abs(s1).max(), 1.0) < 0.005
