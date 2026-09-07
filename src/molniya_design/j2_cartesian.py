"""Independent Cartesian J2 acceleration and osculating propagator.

**Supporting cross-check only** (M4 §11): used to independently verify
the first-order secular rates in :mod:`molniya_design.j2` over a
short-to-medium horizon by numerically integrating the full osculating
Cartesian equations of motion and fitting a secular trend to the
resulting osculating RAAN/argp time series. This module is explicitly
**not** promoted to a production propagator for M4/M5 — the secular
mean-element model in `j2.py` remains primary.

Acceleration model
-------------------
Two-body plus the standard J2 perturbation term, in the form given in the
M4 task description:

    ax = -mu*x/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 1)]
    ay = -mu*y/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 1)]
    az = -mu*z/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 3)]

This is algebraically identical to the more commonly quoted factored form
    a_J2 = -(3/2) J2 (mu/r^2) (Re/r)^2 [ (1-5(z/r)^2) x/r, (1-5(z/r)^2) y/r,
                                          (3-5(z/r)^2) z/r ]
(expand both and the x/y/z terms match term-by-term). Using the task's
form directly here for traceability to the M4 instructions.

Setting J2=0 makes this reduce exactly to the M2 two-body EOM
(:func:`molniya_design.twobody.two_body_eom`) — tested explicitly.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .constants import J2_EARTH, MU_EARTH, R_EARTH


def j2_acceleration(
    r: np.ndarray, mu: float = MU_EARTH, J2: float = J2_EARTH, Re: float = R_EARTH
) -> np.ndarray:
    """Two-body + J2 acceleration (km/s^2) at ECI position r (km)."""
    x, y, z = r
    r_mag = np.linalg.norm(r)
    z2_r2 = (z / r_mag) ** 2
    common = 1.0 - 1.5 * J2 * (Re / r_mag) ** 2 * (5.0 * z2_r2 - 1.0)
    common_z = 1.0 - 1.5 * J2 * (Re / r_mag) ** 2 * (5.0 * z2_r2 - 3.0)

    ax = -mu * x / r_mag**3 * common
    ay = -mu * y / r_mag**3 * common
    az = -mu * z / r_mag**3 * common_z
    return np.array([ax, ay, az])


def j2_eom(
    t: float, y: np.ndarray, mu: float = MU_EARTH, J2: float = J2_EARTH, Re: float = R_EARTH
) -> np.ndarray:
    """y' = f(t, y) for state y=[rx,ry,rz,vx,vy,vz], two-body + J2."""
    r = y[0:3]
    v = y[3:6]
    a = j2_acceleration(r, mu, J2, Re)
    return np.concatenate([v, a])


def propagate_j2_cartesian(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    t_eval: np.ndarray | None = None,
    mu: float = MU_EARTH,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    method: str = "DOP853",
    rtol: float = 1e-12,
    atol: float = 1e-12,
    dense_output: bool = False,
):
    """Numerically propagate the osculating two-body + J2 Cartesian state.
    Thin `solve_ivp` wrapper, mirroring `propagation.propagate`."""
    y0 = np.concatenate([np.asarray(r0, dtype=float), np.asarray(v0, dtype=float)])
    sol = solve_ivp(
        j2_eom, t_span, y0, args=(mu, J2, Re), method=method,
        t_eval=t_eval, rtol=rtol, atol=atol, dense_output=dense_output,
    )
    if not sol.success:
        raise RuntimeError(f"J2 Cartesian propagation failed: {sol.message}")
    return sol
