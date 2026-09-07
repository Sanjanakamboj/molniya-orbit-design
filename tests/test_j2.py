"""Milestone 4 test suite: first-order secular J2 model, critical-
inclination verification, secular ground track, and an independent
Cartesian J2 cross-check.

Covers checklist items A-R from the M4 task description.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    J2_EARTH,
    M1_APOGEE_ALT_KM,
    M1_PERIGEE_ALT_KM,
    M1_PERIOD_S,
    M1_VA_KM_S,
    M1_VP_KM_S,
    MU_EARTH,
    R_EARTH,
)
from molniya_design.elements import coe_to_rv, rv_to_coe
from molniya_design.groundtrack import (
    compute_ground_track,
    compute_ground_track_j2,
    longitude_separation_deg,
)
from molniya_design.j2 import (
    SecularElements,
    compute_secular_rates,
    mean_anomaly_crossing_times,
    propagate_elements_j2_secular,
    reconstruct_eci_state,
)
from molniya_design.j2_cartesian import j2_acceleration, j2_eom
from molniya_design.propagation import find_apsides, propagate, true_to_mean_anomaly
from molniya_design.twobody import two_body_eom

BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))


@pytest.fixture(scope="module")
def elements0():
    return SecularElements(
        a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
        raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG,
    )


# ---------------------------------------------------------------------------
# A. baseline Omega_dot matches independent formula
# ---------------------------------------------------------------------------


def test_A_baseline_raan_dot_independent_formula():
    """Recompute RAAN_dot with a formula coded independently in this test
    (not by re-calling compute_secular_rates), as required by M4 §5."""
    a, e, i_deg = BASE["a_km"], BASE["e"], BASE["i_deg"]
    n_independent = (MU_EARTH / a**3) ** 0.5
    p_independent = a * (1 - e**2)
    i_rad = i_deg * np.pi / 180.0
    raan_dot_independent = -1.5 * J2_EARTH * n_independent * (R_EARTH / p_independent) ** 2 * np.cos(i_rad)
    raan_dot_deg_day_independent = raan_dot_independent * 180.0 / np.pi * 86400.0

    rates = compute_secular_rates(a, e, i_deg)
    assert rates.raan_dot_deg_day == pytest.approx(raan_dot_deg_day_independent, rel=1e-12)
    # M1 target
    assert rates.raan_dot_deg_day == pytest.approx(-0.145135, abs=1e-5)


# ---------------------------------------------------------------------------
# B. baseline omega_dot ~ 0
# ---------------------------------------------------------------------------


def test_B_baseline_argp_dot_near_zero():
    rates = compute_secular_rates(BASE["a_km"], BASE["e"], BASE["i_deg"])
    assert abs(rates.argp_dot_deg_day) < 1e-12


# ---------------------------------------------------------------------------
# C. critical relation 5cos^2(i)-1 ~ 0
# ---------------------------------------------------------------------------


def test_C_critical_relation():
    i_rad = np.radians(BASE["i_deg"])
    assert (5 * np.cos(i_rad) ** 2 - 1) == pytest.approx(0.0, abs=1e-14)


# ---------------------------------------------------------------------------
# D. omega_dot sign flips across critical inclination
# ---------------------------------------------------------------------------


def test_D_argp_dot_sign_flip():
    below = compute_secular_rates(BASE["a_km"], BASE["e"], 60.0).argp_dot_deg_day
    above = compute_secular_rates(BASE["a_km"], BASE["e"], 65.0).argp_dot_deg_day
    assert below > 0
    assert above < 0


# ---------------------------------------------------------------------------
# E. off-critical omega_dot magnitudes are nonzero
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("i_deg", [60.0, 63.0, 64.0, 65.0])
def test_E_off_critical_argp_dot_nonzero(i_deg):
    r = compute_secular_rates(BASE["a_km"], BASE["e"], i_deg)
    assert abs(r.argp_dot_deg_day) > 1e-4


# ---------------------------------------------------------------------------
# F. J2->0 secular model reproduces two-body mean-element evolution
# ---------------------------------------------------------------------------


def test_F_j2_zero_limit_rates(elements0):
    rates = compute_secular_rates(BASE["a_km"], BASE["e"], BASE["i_deg"], J2=0.0)
    n_expected = np.sqrt(MU_EARTH / BASE["a_km"] ** 3)
    assert rates.raan_dot_rad_s == 0.0
    assert rates.argp_dot_rad_s == 0.0
    assert rates.m_dot_rad_s == pytest.approx(n_expected, rel=1e-14)


def test_F_j2_zero_elements_match_two_body_kepler(elements0):
    from molniya_design.propagation import kepler_state_at_time

    t = 0.37 * M1_PERIOD_S
    elements_t = propagate_elements_j2_secular(elements0, t, J2=0.0)
    r_j2, v_j2 = reconstruct_eci_state(elements_t)
    r_kep, v_kep = kepler_state_at_time(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"], t,
    )
    assert np.linalg.norm(r_j2 - r_kep) / np.linalg.norm(r_kep) < 1e-10
    assert np.linalg.norm(v_j2 - v_kep) / np.linalg.norm(v_kep) < 1e-10


# ---------------------------------------------------------------------------
# G. a/e/i remain constant in secular model
# ---------------------------------------------------------------------------


def test_G_a_e_i_constant(elements0):
    for t in [0.0, 1e5, 1e6, 1.2e6]:
        e_t = propagate_elements_j2_secular(elements0, t)
        assert e_t.a_km == elements0.a_km
        assert e_t.e == elements0.e
        assert e_t.i_deg == elements0.i_deg


# ---------------------------------------------------------------------------
# H. Omega/omega/M evolve linearly at their rates
# ---------------------------------------------------------------------------


def test_H_linear_evolution(elements0):
    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg)
    t = 12345.678
    e_t = propagate_elements_j2_secular(elements0, t, wrap=False)
    expected_raan = elements0.raan_deg + rates.raan_dot_deg_day / 86400.0 * t
    expected_argp = elements0.argp_deg + rates.argp_dot_deg_day / 86400.0 * t
    expected_m = elements0.m_deg + rates.m_dot_deg_day / 86400.0 * t
    assert e_t.raan_deg == pytest.approx(expected_raan, rel=1e-10)
    assert e_t.argp_deg == pytest.approx(expected_argp, rel=1e-10)
    assert e_t.m_deg == pytest.approx(expected_m, rel=1e-10)


def test_H_no_mutation_of_input(elements0):
    before = (elements0.a_km, elements0.e, elements0.i_deg, elements0.raan_deg, elements0.argp_deg, elements0.m_deg)
    propagate_elements_j2_secular(elements0, 50000.0)
    after = (elements0.a_km, elements0.e, elements0.i_deg, elements0.raan_deg, elements0.argp_deg, elements0.m_deg)
    assert before == after


# ---------------------------------------------------------------------------
# I. Kepler reconstruction from secular elements is self-consistent
# ---------------------------------------------------------------------------


def test_I_reconstruction_radius_matches_kepler_formula(elements0):
    t = 5000.0
    e_t = propagate_elements_j2_secular(elements0, t)
    r, v = reconstruct_eci_state(e_t)
    r_mag = np.linalg.norm(r)

    M_rad = np.radians(e_t.m_deg)
    from molniya_design.propagation import mean_to_eccentric_anomaly

    E_rad = mean_to_eccentric_anomaly(M_rad, e_t.e)
    r_formula = e_t.a_km * (1.0 - e_t.e * np.cos(E_rad))
    assert r_mag == pytest.approx(r_formula, rel=1e-10)


# ---------------------------------------------------------------------------
# J. J2=0 reconstructed ECI state matches M2 propagation
# ---------------------------------------------------------------------------


def test_J_j2_zero_matches_m2_propagation(elements0):
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    t = 0.71 * M1_PERIOD_S
    sol = propagate(r0, v0, (0.0, t), rtol=1e-13, atol=1e-13)
    r_m2, v_m2 = sol.y[0:3, -1], sol.y[3:6, -1]

    e_t = propagate_elements_j2_secular(elements0, t, J2=0.0)
    r_j2, v_j2 = reconstruct_eci_state(e_t)

    assert np.linalg.norm(r_j2 - r_m2) / np.linalg.norm(r_m2) < 1e-9
    assert np.linalg.norm(v_j2 - v_m2) / np.linalg.norm(v_m2) < 1e-9


# ---------------------------------------------------------------------------
# K. angle wrapping behavior
# ---------------------------------------------------------------------------


def test_K_angle_wrapping(elements0):
    t = 200 * 86400.0  # long horizon, many revolutions -> unwrapped angle >> 360
    e_wrapped = propagate_elements_j2_secular(elements0, t, wrap=True)
    e_unwrapped = propagate_elements_j2_secular(elements0, t, wrap=False)

    assert 0.0 <= e_wrapped.raan_deg < 360.0
    assert 0.0 <= e_wrapped.argp_deg < 360.0
    assert 0.0 <= e_wrapped.m_deg < 360.0

    assert e_wrapped.raan_deg == pytest.approx(e_unwrapped.raan_deg % 360.0, abs=1e-6)
    assert e_wrapped.m_deg == pytest.approx(e_unwrapped.m_deg % 360.0, abs=1e-6)
    # unwrapped M should be very large (many revolutions)
    assert e_unwrapped.m_deg > 720.0


def test_K_array_time_input(elements0):
    ts = np.array([0.0, 1000.0, 5000.0, 10000.0])
    e_t = propagate_elements_j2_secular(elements0, ts)
    assert e_t.m_deg.shape == ts.shape
    assert e_t.raan_deg.shape == ts.shape


# ---------------------------------------------------------------------------
# L. apogee-time analytical calculation
# ---------------------------------------------------------------------------


def test_L_apogee_time_analytical(elements0):
    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg)
    times = mean_anomaly_crossing_times(elements0, 180.0, (0.0, 2.5 * M1_PERIOD_S))
    assert len(times) >= 2
    # first apogee at t=0 (m0=180 deg exactly)
    assert times[0] == pytest.approx(0.0, abs=1e-6)
    # next apogee: full 360 deg of mean anomaly later
    expected_next = 2.0 * np.pi / rates.m_dot_rad_s
    assert times[1] == pytest.approx(expected_next, rel=1e-10)

    # verify: at each returned time, M really is 180 deg (mod 360)
    for t in times:
        e_t = propagate_elements_j2_secular(elements0, t, wrap=True)
        assert e_t.m_deg == pytest.approx(180.0, abs=1e-6)


# ---------------------------------------------------------------------------
# M. two-body vs J2 ground-track divergence is nonzero
# ---------------------------------------------------------------------------


def test_M_ground_track_divergence_nonzero(elements0):
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    t_span = (0.0, 7 * 86400.0)
    gt_2body = compute_ground_track(r0, v0, t_span, dt_s=600.0)
    gt_j2 = compute_ground_track_j2(elements0, t_span, dt_s=600.0)

    lon_diff = np.array(
        [longitude_separation_deg(a, b) for a, b in zip(gt_2body.lon_deg, gt_j2.lon_deg)]
    )
    # at t=0 they coincide; by the end of 7 days they should have diverged measurably
    assert abs(lon_diff[0]) < 1e-6
    assert abs(lon_diff[-1]) > 0.5  # degrees


# ---------------------------------------------------------------------------
# N. baseline omega orientation remains effectively fixed over 14 days
# ---------------------------------------------------------------------------


def test_N_baseline_argp_fixed_over_14_days(elements0):
    e_14 = propagate_elements_j2_secular(elements0, 14 * 86400.0, wrap=False)
    assert e_14.argp_deg == pytest.approx(270.0, abs=1e-9)


# ---------------------------------------------------------------------------
# O. off-critical apogee orientation drifts measurably
# ---------------------------------------------------------------------------


def test_O_off_critical_argp_drifts(elements0):
    elements0_65 = SecularElements(
        a_km=elements0.a_km, e=elements0.e, i_deg=65.0,
        raan_deg=elements0.raan_deg, argp_deg=elements0.argp_deg, m_deg=elements0.m_deg,
    )
    e_14 = propagate_elements_j2_secular(elements0_65, 14 * 86400.0, wrap=False)
    drift = e_14.argp_deg - 270.0
    assert abs(drift) > 0.1  # degrees, clearly nonzero over 14 days
    assert drift < 0  # i > i_crit -> argp_dot < 0 -> omega decreases


# ---------------------------------------------------------------------------
# P. M1/M2/M3 regression tests still pass (spot check via production code)
# ---------------------------------------------------------------------------


def test_P_regression_apsis_still_matches_M1():
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    apsides = find_apsides(r0, v0, (0.0, 1.01 * M1_PERIOD_S), rtol=1e-13, atol=1e-13)
    perigee = next(a for a in apsides if a.kind == "perigee")
    apogee = next(a for a in apsides if a.kind == "apogee" and a.t_s > 1.0)
    assert (perigee.r_km - R_EARTH) == pytest.approx(M1_PERIGEE_ALT_KM, abs=1e-4)
    assert (apogee.r_km - R_EARTH) == pytest.approx(M1_APOGEE_ALT_KM, abs=1e-4)
    assert perigee.v_km_s == pytest.approx(M1_VP_KM_S, abs=1e-6)
    assert apogee.v_km_s == pytest.approx(M1_VA_KM_S, abs=1e-6)


# ---------------------------------------------------------------------------
# Q. units and deg/day conversion sanity
# ---------------------------------------------------------------------------


def test_Q_unit_conversion_sanity():
    rates = compute_secular_rates(BASE["a_km"], BASE["e"], BASE["i_deg"])
    assert rates.raan_dot_deg_day == pytest.approx(np.degrees(rates.raan_dot_rad_s) * 86400.0, rel=1e-14)
    assert rates.argp_dot_deg_day == pytest.approx(np.degrees(rates.argp_dot_rad_s) * 86400.0, rel=1e-14)
    assert rates.m_dot_deg_day == pytest.approx(np.degrees(rates.m_dot_rad_s) * 86400.0, rel=1e-14)
    # deg/day of unperturbed n should be ~721.97 (baseline ~2 rev/day)
    n_deg_day = np.degrees(rates.n_rad_s) * 86400.0
    assert n_deg_day == pytest.approx(721.971295, rel=1e-6)


# ---------------------------------------------------------------------------
# R. Cartesian J2 cross-check: fitted secular rates agree with theory
# ---------------------------------------------------------------------------


def test_R_j2_zero_cartesian_reduces_to_two_body():
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    y0 = np.concatenate([r0, v0])
    d_j2_zero = j2_eom(0.0, y0, J2=0.0)
    d_two_body = two_body_eom(0.0, y0)
    assert np.allclose(d_j2_zero, d_two_body, atol=1e-15)


def test_R_j2_acceleration_matches_two_body_plus_correction():
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    a_j2 = j2_acceleration(r0)
    a_two_body = -MU_EARTH * r0 / np.linalg.norm(r0) ** 3
    # J2 correction should be nonzero but much smaller than the two-body term
    correction = a_j2 - a_two_body
    assert np.linalg.norm(correction) > 0
    assert np.linalg.norm(correction) / np.linalg.norm(a_two_body) < 1e-2


def test_R_fitted_secular_rates_match_theory():
    """Independent cross-check: propagate the full osculating Cartesian
    J2 dynamics over 7 days, sample RAAN/argp at perigee crossings
    (removing most short-period oscillation), and linear-fit the secular
    trend. Compare against the first-order theory rates."""
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"],
    )
    y0 = np.concatenate([r0, v0])

    def rdotv(t, y):
        return np.dot(y[0:3], y[3:6])

    rdotv.direction = 0

    t_span = (0.0, 7 * 86400.0)
    sol = solve_ivp(
        j2_eom, t_span, y0, method="DOP853", rtol=1e-13, atol=1e-13, dense_output=True, events=rdotv,
    )
    event_times = sol.t_events[0]
    states = [sol.sol(t) for t in event_times]
    radii = [np.linalg.norm(s[0:3]) for s in states]
    r_mid = 0.5 * (min(radii) + max(radii))

    perigee_t, perigee_raan, perigee_argp = [], [], []
    for t, s, r in zip(event_times, states, radii):
        if r < r_mid:
            coe = rv_to_coe(s[0:3], s[3:6])
            perigee_t.append(t)
            perigee_raan.append(coe["raan_deg"])
            perigee_argp.append(coe["argp_deg"])

    perigee_t = np.array(perigee_t)
    perigee_raan = np.degrees(np.unwrap(np.radians(perigee_raan)))
    perigee_argp = np.degrees(np.unwrap(np.radians(perigee_argp)))
    assert len(perigee_t) >= 10

    A = np.vstack([perigee_t, np.ones_like(perigee_t)]).T
    raan_slope = np.linalg.lstsq(A, perigee_raan, rcond=None)[0][0] * 86400.0
    argp_slope = np.linalg.lstsq(A, perigee_argp, rcond=None)[0][0] * 86400.0

    theory = compute_secular_rates(BASE["a_km"], BASE["e"], BASE["i_deg"])
    # first-order theory vs a 7-day osculating fit: agreement within 0.1%
    # relative for RAAN_dot (dominant, well-resolved rate)
    assert raan_slope == pytest.approx(theory.raan_dot_deg_day, rel=1e-3)
    # argp_dot theory is ~0; the fitted value should be much smaller than
    # any off-critical rate (>1e-4 deg/day, see test_E), confirming near-zero
    assert abs(argp_slope) < 1e-3
