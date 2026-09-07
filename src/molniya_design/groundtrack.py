"""Two-body ECI propagation -> ECEF -> geocentric ground track.

No J2. This module only chains together M2's two-body propagator
(:mod:`molniya_design.propagation`) with M3's frame transformation
(:mod:`molniya_design.frames`) and the M2 apsis detector — it adds no new
dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MU_EARTH
from .frames import ecef_to_geocentric_latlon, eci_to_ecef
from .propagation import Apsis, find_apsides, propagate


@dataclass
class GroundTrack:
    t_s: np.ndarray
    r_eci: np.ndarray  # (3, N) km
    r_ecef: np.ndarray  # (3, N) km
    lat_deg: np.ndarray  # (N,)
    lon_deg: np.ndarray  # (N,)


def compute_ground_track(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    dt_s: float,
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> GroundTrack:
    """Propagate (two-body, no J2) and compute the Earth-fixed ground
    track at uniform time steps ``dt_s`` over ``t_span``."""
    n = int(round((t_span[1] - t_span[0]) / dt_s)) + 1
    t_eval = np.linspace(t_span[0], t_span[1], n)
    sol = propagate(r0, v0, t_span, t_eval=t_eval, mu=mu, rtol=rtol, atol=atol)

    r_eci = sol.y[0:3, :]
    r_ecef = np.zeros_like(r_eci)
    lat = np.zeros(sol.t.shape)
    lon = np.zeros(sol.t.shape)
    for k, t in enumerate(sol.t):
        r_ecef[:, k] = eci_to_ecef(r_eci[:, k], t, theta_g0_rad)
        lat[k], lon[k] = ecef_to_geocentric_latlon(r_ecef[:, k])

    return GroundTrack(t_s=sol.t, r_eci=r_eci, r_ecef=r_ecef, lat_deg=lat, lon_deg=lon)


@dataclass
class ApsisGroundPoint:
    t_s: float
    kind: str
    lat_deg: float
    lon_deg: float
    r_km: float


def apsis_ground_points(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> list[ApsisGroundPoint]:
    """Detect apsides (M2 method) and report their Earth-fixed lat/lon."""
    apsides: list[Apsis] = find_apsides(r0, v0, t_span, mu=mu, rtol=rtol, atol=atol)
    points = []
    for a in apsides:
        r_ecef = eci_to_ecef(a.state[0:3], a.t_s, theta_g0_rad)
        lat, lon = ecef_to_geocentric_latlon(r_ecef)
        points.append(ApsisGroundPoint(t_s=a.t_s, kind=a.kind, lat_deg=lat, lon_deg=lon, r_km=a.r_km))
    return points


def longitude_separation_deg(lon1_deg: float, lon2_deg: float) -> float:
    """Signed shortest angular separation lon2 - lon1, wrapped to
    (-180, 180] deg, avoiding the +180/-180 discontinuity when comparing
    successive ground-track passes."""
    d = (lon2_deg - lon1_deg + 180.0) % 360.0 - 180.0
    if d == -180.0:
        d = 180.0
    return d
