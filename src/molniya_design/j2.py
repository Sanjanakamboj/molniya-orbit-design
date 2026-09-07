"""First-order secular J2 mean-element model.

Scope (M4): a clean **mean-element secular** model — not a full
osculating-element (Brouwer/Lyddane) theory and not a high-fidelity
numerical propagator. An optional short-horizon Cartesian J2 cross-check
lives in :mod:`molniya_design.j2_cartesian` and is used only to
independently verify the secular rates below, per the M4 scope guard.

Secular rate formulas (standard first-order, e.g. Vallado *Fundamentals
of Astrodynamics and Applications*)
-----------------------------------------------------------------------
For a Keplerian mean-element orbit with

    n = sqrt(mu / a^3)              (unperturbed mean motion, rad/s)
    p = a (1 - e^2)                 (semilatus rectum, km)

the first-order secular rates used throughout this module are:

    RAAN_dot = -(3/2) * J2 * n * (Re/p)^2 * cos(i)

    argp_dot =  (3/4) * J2 * n * (Re/p)^2 * (5 cos^2(i) - 1)

    M_dot    = n + (3/4) * J2 * n * (Re/p)^2 * sqrt(1-e^2) * (3 cos^2(i) - 1)

**M-dot convention (explicit):** the mean anomaly does *not* advance at
the unperturbed rate `n` alone. The secular J2 perturbation also shifts
the mean anomaly's own secular rate (sometimes folded into a "drag-like"
period change); the correction term above is the standard first-order
result for the secular rate of mean anomaly due to J2 (identical in
structure to the RAAN/argp rate derivations, sharing the same
`J2 * n * (Re/p)^2` prefactor). This is tested independently in
`tests/test_j2.py` (checklist item F: the J2->0 limit must recover
`M_dot = n` exactly, i.e. plain two-body mean-anomaly advance).

All three rates share the common factor `J2 * n * (Re/p)^2`, which is
computed once per :func:`compute_secular_rates` call.

Units: km, s, radians internally; degrees only at the public dataclass
boundary (matching the convention established in `elements.py`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import J2_EARTH, MU_EARTH, R_EARTH
from .elements import coe_to_rv
from .propagation import (
    eccentric_to_true_anomaly,
    mean_to_eccentric_anomaly,
)


@dataclass(frozen=True)
class SecularRates:
    """First-order secular rates, all in rad/s (see module docstring)."""

    n_rad_s: float
    raan_dot_rad_s: float
    argp_dot_rad_s: float
    m_dot_rad_s: float

    @property
    def raan_dot_deg_day(self) -> float:
        return np.degrees(self.raan_dot_rad_s) * 86400.0

    @property
    def argp_dot_deg_day(self) -> float:
        return np.degrees(self.argp_dot_rad_s) * 86400.0

    @property
    def m_dot_deg_day(self) -> float:
        return np.degrees(self.m_dot_rad_s) * 86400.0


def compute_secular_rates(
    a_km: float,
    e: float,
    i_deg: float,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    mu: float = MU_EARTH,
) -> SecularRates:
    """First-order secular RAAN/argp/mean-anomaly rates for a Keplerian
    mean-element orbit (a, e, i held fixed). See module docstring for the
    exact formulas and M-dot convention."""
    n = np.sqrt(mu / a_km**3)
    p = a_km * (1.0 - e**2)
    i = np.radians(i_deg)
    cos_i = np.cos(i)
    cos2_i = cos_i**2

    common = J2 * n * (Re / p) ** 2

    raan_dot = -1.5 * common * cos_i
    argp_dot = 0.75 * common * (5.0 * cos2_i - 1.0)
    m_dot = n + 0.75 * common * np.sqrt(1.0 - e**2) * (3.0 * cos2_i - 1.0)

    return SecularRates(
        n_rad_s=n, raan_dot_rad_s=raan_dot, argp_dot_rad_s=argp_dot, m_dot_rad_s=m_dot
    )


@dataclass(frozen=True)
class SecularElements:
    """Mean classical elements at a given time. a, e, i are held fixed in
    the M4 secular model; raan_deg/argp_deg/m_deg are the (possibly
    unwrapped) mean angles. Immutable — propagation returns a new
    instance, never mutates its input (M4 §6 requirement)."""

    a_km: float
    e: float
    i_deg: float
    raan_deg: float
    argp_deg: float
    m_deg: float


def propagate_elements_j2_secular(
    elements0: SecularElements,
    t_s,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    mu: float = MU_EARTH,
    wrap: bool = True,
) -> SecularElements:
    """Propagate mean elements under the first-order secular J2 model.

    a, e, i are constant. RAAN, argp, and M evolve linearly at the rates
    from :func:`compute_secular_rates`, computed once from ``elements0``
    (a, e, i do not change, so the rates are constant over the whole
    propagation — this is the defining property of the *secular* model).

    ``t_s`` may be a scalar or a numpy array (both supported via
    broadcasting on the returned angle fields); the returned
    ``SecularElements`` fields will then themselves be scalars or arrays.

    Internally the angles are advanced **unwrapped** (raan0 + rate*t, with
    no modulo applied at any intermediate step) and are wrapped to
    [0, 360) degrees only at the very end if ``wrap=True`` (the default).
    Pass ``wrap=False`` to get the continuous (unwrapped) angle, useful
    for e.g. plotting long-horizon drift without a wraparound
    discontinuity, or for exactly recovering the linear secular trend in
    a regression/tests.

    Does not mutate ``elements0``.
    """
    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg, J2, Re, mu)

    t_s = np.asarray(t_s, dtype=float)

    raan = np.radians(elements0.raan_deg) + rates.raan_dot_rad_s * t_s
    argp = np.radians(elements0.argp_deg) + rates.argp_dot_rad_s * t_s
    m = np.radians(elements0.m_deg) + rates.m_dot_rad_s * t_s

    raan_deg = np.degrees(raan)
    argp_deg = np.degrees(argp)
    m_deg = np.degrees(m)

    if wrap:
        raan_deg = raan_deg % 360.0
        argp_deg = argp_deg % 360.0
        m_deg = m_deg % 360.0

    # convert back to plain float when t_s was a 0-d/scalar input, for a
    # clean, non-surprising scalar-in/scalar-out API
    if t_s.ndim == 0:
        raan_deg = float(raan_deg)
        argp_deg = float(argp_deg)
        m_deg = float(m_deg)

    return SecularElements(
        a_km=elements0.a_km,
        e=elements0.e,
        i_deg=elements0.i_deg,
        raan_deg=raan_deg,
        argp_deg=argp_deg,
        m_deg=m_deg,
    )


def reconstruct_eci_state(
    elements: SecularElements, mu: float = MU_EARTH
) -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct an ECI (r, v) state from mean elements at a single
    instant (scalar fields only) by solving Kepler's equation for the
    eccentric/true anomaly and reusing :func:`elements.coe_to_rv` — no
    duplicated orbital-geometry code, per M4 §7."""
    M_rad = np.radians(elements.m_deg) % (2.0 * np.pi)
    E_rad = mean_to_eccentric_anomaly(M_rad, elements.e)
    nu_rad = eccentric_to_true_anomaly(E_rad, elements.e)
    return coe_to_rv(
        elements.a_km, elements.e, elements.i_deg, elements.raan_deg,
        elements.argp_deg, np.degrees(nu_rad), mu=mu,
    )


def mean_anomaly_crossing_times(
    elements0: SecularElements,
    target_m_deg: float,
    t_span: tuple[float, float],
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    mu: float = MU_EARTH,
) -> np.ndarray:
    """Analytical times within ``t_span`` at which the secularly
    propagated mean anomaly equals ``target_m_deg`` (mod 360), e.g.
    ``target_m_deg=180`` for apogee or ``0`` for perigee.

    Solved directly from the linear M(t) = M0 + M_dot*t secular relation
    (no event detection / numerical search needed, since M(t) is exactly
    linear in this model) — the robust, preferred method per M4 §9.
    """
    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg, J2, Re, mu)
    m_dot = rates.m_dot_rad_s
    if m_dot <= 0:
        raise ValueError("mean anomaly rate must be positive")

    m0 = np.radians(elements0.m_deg)
    target = np.radians(target_m_deg)

    # t_k = (target - m0 + 2*pi*k) / m_dot for integer k, restricted to t_span
    t_start, t_end = t_span
    k_min = int(np.floor((m_dot * t_start - target + m0) / (2.0 * np.pi))) - 1
    k_max = int(np.ceil((m_dot * t_end - target + m0) / (2.0 * np.pi))) + 1

    times = []
    for k in range(k_min, k_max + 1):
        t_k = (target - m0 + 2.0 * np.pi * k) / m_dot
        if t_start - 1e-6 <= t_k <= t_end + 1e-6:
            times.append(t_k)
    return np.array(sorted(times))
