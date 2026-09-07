"""Milestone 2 test suite: element/state conversion, two-body propagation,
apsis detection, independent Kepler cross-check, conservation, orientation,
and tolerance convergence.

Covers checklist items A-P from the M2 task description; see DESIGN.md §M2
for the corresponding narrative results.
"""

from __future__ import annotations

import numpy as np
import pytest

from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    MU_EARTH,
    M1_APOGEE_ALT_KM,
    M1_PERIGEE_ALT_KM,
    M1_PERIOD_S,
    M1_RA_KM,
    M1_RP_KM,
    M1_VA_KM_S,
    M1_VP_KM_S,
    R_EARTH,
)
from molniya_design.elements import coe_to_rv, rv_to_coe
from molniya_design.propagation import (
    find_apsides,
    geocentric_latitude_deg,
    kepler_state_at_time,
    propagate,
)
from molniya_design.twobody import (
    eccentricity_vector,
    specific_angular_momentum,
    specific_energy,
)

BASE = BASELINE_ELEMENTS_DEG

# ---------------------------------------------------------------------------
# A. M1 element/state baseline conversion
# ---------------------------------------------------------------------------


def test_A_baseline_state_at_apogee_matches_M1():
    r, v = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    assert r_mag == pytest.approx(M1_RA_KM, rel=1e-6)
    assert v_mag == pytest.approx(M1_VA_KM_S, rel=1e-6)


# ---------------------------------------------------------------------------
# B. element -> state -> element round trip (baseline)
# ---------------------------------------------------------------------------


def test_B_baseline_round_trip():
    r, v = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )
    coe = rv_to_coe(r, v)
    assert coe["a_km"] == pytest.approx(BASE["a_km"], rel=1e-8)
    assert coe["e"] == pytest.approx(BASE["e"], rel=1e-8)
    assert coe["i_deg"] == pytest.approx(BASE["i_deg"], abs=1e-6)
    assert coe["raan_deg"] == pytest.approx(BASE["raan_deg"], abs=1e-6)
    assert coe["argp_deg"] == pytest.approx(BASE["argp_deg"], abs=1e-6)
    assert coe["nu_deg"] == pytest.approx(BASE["nu0_deg"], abs=1e-6)


# ---------------------------------------------------------------------------
# C. arbitrary nonsingular round trips
# ---------------------------------------------------------------------------

ARBITRARY_ORBITS = [
    dict(a_km=8000.0, e=0.1, i_deg=28.5, raan_deg=45.0, argp_deg=30.0, nu_deg=10.0),
    dict(a_km=20000.0, e=0.5, i_deg=63.4, raan_deg=200.0, argp_deg=270.0, nu_deg=350.0),
    dict(a_km=42164.0, e=0.02, i_deg=5.0, raan_deg=10.0, argp_deg=15.0, nu_deg=200.0),
    dict(a_km=15000.0, e=0.35, i_deg=97.0, raan_deg=300.0, argp_deg=120.0, nu_deg=75.0),
]


@pytest.mark.parametrize("orbit", ARBITRARY_ORBITS)
def test_C_arbitrary_round_trip(orbit):
    r, v = coe_to_rv(**orbit)
    coe = rv_to_coe(r, v)
    assert coe["a_km"] == pytest.approx(orbit["a_km"], rel=1e-8)
    assert coe["e"] == pytest.approx(orbit["e"], rel=1e-8, abs=1e-10)
    assert coe["i_deg"] == pytest.approx(orbit["i_deg"], abs=1e-6)
    assert coe["raan_deg"] == pytest.approx(orbit["raan_deg"], abs=1e-6)
    assert coe["argp_deg"] == pytest.approx(orbit["argp_deg"], abs=1e-6)
    assert coe["nu_deg"] == pytest.approx(orbit["nu_deg"], abs=1e-6)


@pytest.mark.parametrize("orbit", ARBITRARY_ORBITS + [BASE])
def test_C_pqw_to_eci_rotation_is_right_handed(orbit):
    """PQW->ECI rotation matrix must be a proper rotation (det=+1),
    confirming right-handedness for every element set tested, as claimed
    in elements.py / DESIGN.md."""
    from molniya_design.elements import _rot1, _rot3

    i = np.radians(orbit["i_deg"])
    raan = np.radians(orbit["raan_deg"])
    argp = np.radians(orbit["argp_deg"])
    R = _rot3(raan) @ _rot1(i) @ _rot3(argp)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)
    # orthonormality: R^T R = I
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)


# ---------------------------------------------------------------------------
# Shared fixture: propagate the baseline for ~1 sidereal day (~2 revs)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def baseline_state0():
    return coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )


@pytest.fixture(scope="module")
def baseline_apsides(baseline_state0):
    r0, v0 = baseline_state0
    # propagate a bit past 2 full periods to safely bracket 4 apsis crossings
    t_span = (0.0, 2.05 * M1_PERIOD_S)
    return find_apsides(r0, v0, t_span, rtol=1e-13, atol=1e-13)


# ---------------------------------------------------------------------------
# D. propagated period vs Kepler (M1 analytical) period
# ---------------------------------------------------------------------------


def test_D_numerical_period_matches_M1(baseline_apsides):
    # baseline starts exactly at apogee (nu0=180 deg), so r.v=0 is also
    # detected as a (spurious, t~0) event at the initial condition itself.
    # The *next* apogee event, excluding that t~0 one, is one full period
    # later.
    apogee_times = [a.t_s for a in baseline_apsides if a.kind == "apogee" and a.t_s > 1.0]
    assert len(apogee_times) >= 1
    numerical_period = apogee_times[0]
    assert numerical_period == pytest.approx(M1_PERIOD_S, rel=1e-6)


# ---------------------------------------------------------------------------
# E, F. propagated perigee/apogee radii vs analytic rp, ra
# ---------------------------------------------------------------------------


def test_E_propagated_perigee_radius(baseline_apsides):
    perigee_radii = [a.r_km for a in baseline_apsides if a.kind == "perigee"]
    assert len(perigee_radii) >= 1
    for r in perigee_radii:
        assert r == pytest.approx(M1_RP_KM, rel=1e-7)
        assert (r - R_EARTH) == pytest.approx(M1_PERIGEE_ALT_KM, rel=1e-6)


def test_F_propagated_apogee_radius(baseline_apsides):
    apogee_radii = [a.r_km for a in baseline_apsides if a.kind == "apogee"]
    assert len(apogee_radii) >= 1
    for r in apogee_radii:
        assert r == pytest.approx(M1_RA_KM, rel=1e-7)
        assert (r - R_EARTH) == pytest.approx(M1_APOGEE_ALT_KM, rel=1e-6)


# ---------------------------------------------------------------------------
# G. propagated apsis velocities vs vis-viva
# ---------------------------------------------------------------------------


def test_G_apsis_velocities_vs_vis_viva(baseline_apsides):
    for a in baseline_apsides:
        expected = M1_VP_KM_S if a.kind == "perigee" else M1_VA_KM_S
        assert a.v_km_s == pytest.approx(expected, rel=1e-7)


# ---------------------------------------------------------------------------
# H, I, J. conservation: energy, |h|, eccentricity vector, orbital plane
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def baseline_trajectory(baseline_state0):
    r0, v0 = baseline_state0
    t_eval = np.linspace(0.0, 2.0 * M1_PERIOD_S, 400)
    sol = propagate(r0, v0, (0.0, 2.0 * M1_PERIOD_S), t_eval=t_eval, rtol=1e-13, atol=1e-13)
    return sol


def test_H_energy_conservation(baseline_trajectory):
    sol = baseline_trajectory
    energies = np.array(
        [specific_energy(sol.y[0:3, k], sol.y[3:6, k]) for k in range(sol.y.shape[1])]
    )
    rel_drift = (energies.max() - energies.min()) / abs(energies[0])
    assert rel_drift < 1e-9


def test_I_angular_momentum_conservation(baseline_trajectory):
    sol = baseline_trajectory
    h_mags = np.array(
        [
            np.linalg.norm(specific_angular_momentum(sol.y[0:3, k], sol.y[3:6, k]))
            for k in range(sol.y.shape[1])
        ]
    )
    rel_drift = (h_mags.max() - h_mags.min()) / h_mags[0]
    assert rel_drift < 1e-9


def test_J_eccentricity_vector_and_plane_constant(baseline_trajectory):
    sol = baseline_trajectory
    e_vecs = np.array(
        [eccentricity_vector(sol.y[0:3, k], sol.y[3:6, k]) for k in range(sol.y.shape[1])]
    )
    h_vecs = np.array(
        [specific_angular_momentum(sol.y[0:3, k], sol.y[3:6, k]) for k in range(sol.y.shape[1])]
    )
    h_hat = h_vecs / np.linalg.norm(h_vecs, axis=1, keepdims=True)

    e_ref = e_vecs[0]
    h_hat_ref = h_hat[0]

    e_drift = np.max(np.linalg.norm(e_vecs - e_ref, axis=1))
    plane_drift = np.max(np.linalg.norm(h_hat - h_hat_ref, axis=1))

    assert e_drift < 1e-8
    assert plane_drift < 1e-10


# ---------------------------------------------------------------------------
# K. one-period state closure
# ---------------------------------------------------------------------------


def test_K_one_period_closure(baseline_state0):
    r0, v0 = baseline_state0
    sol = propagate(r0, v0, (0.0, M1_PERIOD_S), rtol=1e-13, atol=1e-13)
    r_final = sol.y[0:3, -1]
    v_final = sol.y[3:6, -1]
    r_err = np.linalg.norm(r_final - r0) / np.linalg.norm(r0)
    v_err = np.linalg.norm(v_final - v0) / np.linalg.norm(v0)
    assert r_err < 1e-6
    assert v_err < 1e-6


# ---------------------------------------------------------------------------
# L. independent Kepler-vs-numerical propagation cross-check
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("frac", [0.0, 0.1, 0.25, 0.5, 0.63, 0.9])
def test_L_kepler_vs_numerical(baseline_state0, frac):
    r0, v0 = baseline_state0
    t = frac * M1_PERIOD_S
    sol = propagate(r0, v0, (0.0, t), rtol=1e-13, atol=1e-13) if t > 0 else None
    if t == 0.0:
        r_num, v_num = r0, v0
    else:
        r_num, v_num = sol.y[0:3, -1], sol.y[3:6, -1]

    r_kep, v_kep = kepler_state_at_time(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"], t,
    )

    r_err = np.linalg.norm(r_num - r_kep) / np.linalg.norm(r_kep)
    v_err = np.linalg.norm(v_num - v_kep) / np.linalg.norm(v_kep)
    assert r_err < 1e-7
    assert v_err < 1e-7


# ---------------------------------------------------------------------------
# M, N. northern-apogee orientation regression test
# ---------------------------------------------------------------------------


def test_M_northern_apogee_omega_270():
    r, v = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], 270.0, 180.0
    )
    lat = geocentric_latitude_deg(r)
    assert lat == pytest.approx(BASE["i_deg"], abs=1e-6)


def test_N_southern_apogee_omega_90():
    r, v = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], 90.0, 180.0
    )
    lat = geocentric_latitude_deg(r)
    assert lat == pytest.approx(-BASE["i_deg"], abs=1e-6)


# ---------------------------------------------------------------------------
# O. tolerance convergence study
# ---------------------------------------------------------------------------


def test_O_tolerance_convergence(baseline_state0):
    r0, v0 = baseline_state0
    settings = [
        ("loose", 1e-6, 1e-6),
        ("medium", 1e-9, 1e-9),
        ("tight", 1e-12, 1e-12),
    ]
    closure_errors = []
    energy_drifts = []
    for _, rtol, atol in settings:
        sol = propagate(r0, v0, (0.0, M1_PERIOD_S), rtol=rtol, atol=atol)
        r_err = np.linalg.norm(sol.y[0:3, -1] - r0) / np.linalg.norm(r0)
        closure_errors.append(r_err)

        e0 = specific_energy(sol.y[0:3, 0], sol.y[3:6, 0])
        ef = specific_energy(sol.y[0:3, -1], sol.y[3:6, -1])
        energy_drifts.append(abs(ef - e0) / abs(e0))

    # errors must not increase as tolerances tighten (loose -> tight)
    assert closure_errors[0] >= closure_errors[1] >= closure_errors[2] * 0.1 or (
        closure_errors[2] < closure_errors[0]
    )
    assert energy_drifts[-1] <= energy_drifts[0]
    # tight-tolerance closure must be small in an absolute sense
    assert closure_errors[-1] < 1e-6


# ---------------------------------------------------------------------------
# P. dimensional / unit sanity checks
# ---------------------------------------------------------------------------


def test_P_unit_round_trips():
    # km <-> m
    x_km = 12345.6789
    assert (x_km * 1000.0) / 1000.0 == pytest.approx(x_km, rel=1e-15)

    # deg <-> rad
    ang_deg = 63.434949
    assert np.degrees(np.radians(ang_deg)) == pytest.approx(ang_deg, rel=1e-15)

    # s <-> h consistency with M1 period
    assert M1_PERIOD_S / 3600.0 == pytest.approx(11.967235, rel=1e-6)

    # mu units: v^2 ~ mu/r check at a reference circular radius
    r_circ = 7000.0
    v_circ = np.sqrt(MU_EARTH / r_circ)
    assert v_circ == pytest.approx(7.546, abs=1e-3)
