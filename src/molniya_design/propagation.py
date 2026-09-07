"""Two-body numerical propagation, apsis detection, and an independent
Kepler (eccentric/mean-anomaly) time-of-flight cross-check.

Two independent verification paths are provided deliberately:

1. **Numerical propagation** (``propagate``, ``find_apsides``) integrates
   the Cartesian two-body ODE from :mod:`molniya_design.twobody` with
   ``scipy.integrate.solve_ivp``.
2. **Analytical Kepler time-of-flight** (``kepler_state_at_time``) solves
   Kepler's equation directly from the classical elements and does **not**
   call ``solve_ivp`` at all — it is a genuinely independent check of the
   propagator, not a wrapper around it (M2 §6 requirement).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .constants import MU_EARTH
from .elements import coe_to_rv
from .twobody import two_body_eom

# ---------------------------------------------------------------------------
# 1. Numerical two-body propagation
# ---------------------------------------------------------------------------


def propagate(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    t_eval: np.ndarray | None = None,
    mu: float = MU_EARTH,
    method: str = "DOP853",
    rtol: float = 1e-12,
    atol: float = 1e-12,
    dense_output: bool = False,
):
    """Numerically propagate a Cartesian two-body state.

    Parameters
    ----------
    r0, v0 : initial position (km) and velocity (km/s), ECI.
    t_span : (t0, tf) in seconds.
    t_eval : optional array of output times, s.
    mu : gravitational parameter, km^3/s^2.
    method : scipy.integrate.solve_ivp integrator. DOP853 (explicit
        Runge-Kutta order 8(5,3)) is used by default for high accuracy on
        this smooth, non-stiff two-body field.
    rtol, atol : integrator tolerances, explicit and configurable (used in
        the §9 convergence study with looser/tighter values).
    dense_output : if True, returns a continuous solution (``sol.sol``)
        usable for event-independent interpolation.

    Returns
    -------
    scipy.integrate.OdeResult (the raw ``solve_ivp`` return value).
    """
    y0 = np.concatenate([np.asarray(r0, dtype=float), np.asarray(v0, dtype=float)])
    sol = solve_ivp(
        two_body_eom,
        t_span,
        y0,
        args=(mu,),
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        dense_output=dense_output,
    )
    if not sol.success:
        raise RuntimeError(f"Propagation failed: {sol.message}")
    return sol


@dataclass
class Apsis:
    t_s: float
    kind: str  # "perigee" or "apogee"
    r_km: float
    v_km_s: float
    state: np.ndarray  # [rx, ry, rz, vx, vy, vz]


def find_apsides(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    mu: float = MU_EARTH,
    method: str = "DOP853",
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[Apsis]:
    """Detect apsides (radial-velocity zero crossings) along a propagated
    trajectory, without assuming or re-using the input elements.

    Method: locate zero crossings of ``r . v`` (radial-velocity sign
    change) via ``solve_ivp`` events, then classify each detected apsis as
    perigee or apogee purely from the *relative* radii of the detected
    crossings (smaller = perigee, larger = apogee) — no analytic rp/ra is
    referenced in the classification itself, only in the caller's
    cross-check.
    """

    def rdotv(t, y, mu=mu):
        r = y[0:3]
        v = y[3:6]
        return np.dot(r, v)

    rdotv.direction = 0

    y0 = np.concatenate([np.asarray(r0, dtype=float), np.asarray(v0, dtype=float)])
    sol = solve_ivp(
        two_body_eom,
        t_span,
        y0,
        args=(mu,),
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=True,
        events=rdotv,
    )
    if not sol.success:
        raise RuntimeError(f"Propagation failed: {sol.message}")

    event_times = sol.t_events[0]
    if len(event_times) == 0:
        return []

    states = [sol.sol(t) for t in event_times]
    radii = [np.linalg.norm(s[0:3]) for s in states]
    r_mid = 0.5 * (min(radii) + max(radii))

    apsides = []
    for t, s, r in zip(event_times, states, radii):
        kind = "perigee" if r < r_mid else "apogee"
        v = np.linalg.norm(s[3:6])
        apsides.append(Apsis(t_s=t, kind=kind, r_km=r, v_km_s=v, state=s))
    return apsides


# ---------------------------------------------------------------------------
# 2. Independent analytical Kepler time-of-flight check
# ---------------------------------------------------------------------------


def true_to_mean_anomaly(nu_rad: float, e: float) -> float:
    """True anomaly -> mean anomaly via eccentric anomaly, radians."""
    E = 2.0 * np.arctan2(
        np.sqrt(1 - e) * np.sin(nu_rad / 2.0), np.sqrt(1 + e) * np.cos(nu_rad / 2.0)
    )
    M = E - e * np.sin(E)
    return M % (2.0 * np.pi)


def mean_to_eccentric_anomaly(
    M_rad: float, e: float, tol: float = 1e-14, maxiter: int = 100
) -> float:
    """Solve Kepler's equation M = E - e sin E for E via Newton iteration."""
    M_rad = M_rad % (2.0 * np.pi)
    E = M_rad if e < 0.8 else np.pi  # standard starting-value heuristic
    for _ in range(maxiter):
        f = E - e * np.sin(E) - M_rad
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    else:
        raise RuntimeError("Kepler's equation did not converge")
    return E % (2.0 * np.pi)


def mean_to_eccentric_anomaly_array(
    M_rad: np.ndarray, e: float, tol: float = 1e-14, maxiter: int = 100
) -> np.ndarray:
    """Vectorized Newton solve of Kepler's equation for an array of mean
    anomalies (M5 performance primitive — regional-grid coverage needs
    thousands of samples). Algebraically identical Newton iteration to
    the scalar :func:`mean_to_eccentric_anomaly`; kept as a separate
    function rather than generalizing the scalar one, so the original
    M2 scalar contract (and its tests) is left completely untouched.
    """
    M_rad = np.asarray(M_rad, dtype=float) % (2.0 * np.pi)
    E = np.where(e < 0.8, M_rad, np.pi * np.ones_like(M_rad))
    for _ in range(maxiter):
        f = E - e * np.sin(E) - M_rad
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E = E + dE
        if np.max(np.abs(dE)) < tol:
            break
    else:
        raise RuntimeError("Kepler's equation did not converge (array form)")
    return E % (2.0 * np.pi)


def eccentric_to_true_anomaly(E_rad: float, e: float) -> float:
    nu = 2.0 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E_rad / 2.0), np.sqrt(1 - e) * np.cos(E_rad / 2.0)
    )
    return nu % (2.0 * np.pi)


def kepler_state_at_time(
    a_km: float,
    e: float,
    i_deg: float,
    raan_deg: float,
    argp_deg: float,
    nu0_deg: float,
    t_s: float,
    mu: float = MU_EARTH,
) -> tuple[np.ndarray, np.ndarray]:
    """Independent analytical propagation: elements at epoch (nu0 at t=0)
    -> state at time t, purely via Kepler's equation and the perifocal
    rotation in :func:`molniya_design.elements.coe_to_rv`. Does **not**
    call :func:`propagate` / ``solve_ivp`` — an independent check path.
    """
    n = np.sqrt(mu / a_km**3)
    M0 = true_to_mean_anomaly(np.radians(nu0_deg), e)
    M_t = (M0 + n * t_s) % (2.0 * np.pi)
    E_t = mean_to_eccentric_anomaly(M_t, e)
    nu_t = eccentric_to_true_anomaly(E_t, e)
    return coe_to_rv(a_km, e, i_deg, raan_deg, argp_deg, np.degrees(nu_t), mu=mu)


# ---------------------------------------------------------------------------
# 3. Small geometry helper (latitude only — no ECEF/longitude; that's M3)
# ---------------------------------------------------------------------------


def geocentric_latitude_deg(r_eci: np.ndarray) -> float:
    """Geocentric latitude of a position vector, degrees.

    Latitude (unlike longitude) does not depend on Earth's rotation about
    the shared Z axis, so this is identical whether ``r_eci`` is expressed
    in ECI or ECEF — no Earth-rotation/GMST transform (M3 scope) is
    needed to compute it.
    """
    r_eci = np.asarray(r_eci, dtype=float)
    r_mag = np.linalg.norm(r_eci)
    return np.degrees(np.arcsin(r_eci[2] / r_mag))
