"""ECI <-> ECEF transformation, and spherical-Earth geocentric lat/lon.

Earth rotation model (M3 scope, no J2, no real UTC epoch)
-----------------------------------------------------------
A **constant-rate** Earth rotation is used:

    theta_G(t) = theta_G0 + omega_E * t

- ``theta_G0 = 0`` at ``t = 0`` — this is an **engineering reference
  epoch**, not a real UTC/GMST epoch. It simply says "ECEF and ECI axes
  are coincident at the start of the propagated timeline." Tying this to
  a real calendar date/GMST is out of scope for this project (DESIGN.md
  §0.4 already flags this as a genuine limitation shared by every
  Earth-fixed-longitude result in this milestone).
- ``omega_E = 2*pi / SIDEREAL_DAY_S`` — the same sidereal rotation rate
  that was already baked into the M1 half-sidereal-day period design
  choice, so the ground-track repeat analysis in this module is
  self-consistent with the M1/M2 baseline by construction.

Rotation sign convention
-------------------------
The ECEF frame rotates **eastward** (counterclockwise viewed from +Z/
north) relative to ECI, at rate ``omega_E``, matching Earth's actual
prograde rotation. Equivalently: a point that is fixed in the ECEF frame
has ECI coordinates equal to its ECEF coordinates rotated by **+theta_G**
about Z (active rotation, using the same ``_rot3`` convention as
``elements.py``):

    r_ECI  = R3(+theta_G) @ r_ECEF
    r_ECEF = R3(-theta_G) @ r_ECI      (since R3 is orthogonal, R3(-x) = R3(x)^T)

This sign is documented here explicitly, tested at theta=0 (identity),
tested for round-trip consistency, and tested against known 90-degree
axis rotations, so longitude signs cannot silently flip in later
milestones.

Geocentric latitude/longitude (spherical Earth only)
-------------------------------------------------------
    lat = asin(z / r)            in [-90, +90] deg
    lon = atan2(y, x)             wrapped to [-180, +180) deg

This is **geocentric** latitude on a **spherical** Earth model — explicitly
*not* WGS-84 geodetic latitude (which would additionally require Earth's
oblateness/flattening; out of scope here and in the M1 limitations list).
"""

from __future__ import annotations

import numpy as np

from .constants import R_EARTH, SIDEREAL_DAY_S

EARTH_ROTATION_RATE_RAD_S = 2.0 * np.pi / SIDEREAL_DAY_S
"""omega_E, rad/s — Earth's rotation rate, derived from the same sidereal
day used throughout M1 (SIDEREAL_DAY_S), not a separately chosen value."""


def _rot3(theta_rad: float) -> np.ndarray:
    """Right-handed active rotation of a vector by +theta about Z.

    Duplicated (not imported) from elements.py deliberately: frames.py is
    the ECI<->ECEF boundary and should not depend on elements.py's PQW/COE
    machinery, keeping the two rotation conventions independently
    inspectable even though they use the same mathematical form.
    """
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def theta_g(t_s: float, theta_g0_rad: float = 0.0) -> float:
    """Greenwich-reference rotation angle at time t (s), radians, wrapped
    to [0, 2*pi)."""
    return (theta_g0_rad + EARTH_ROTATION_RATE_RAD_S * t_s) % (2.0 * np.pi)


def eci_to_ecef(r_eci: np.ndarray, t_s: float, theta_g0_rad: float = 0.0) -> np.ndarray:
    """Rotate an ECI position vector (km) into ECEF at time t (s)."""
    theta = theta_g(t_s, theta_g0_rad)
    R = _rot3(-theta)
    return R @ np.asarray(r_eci, dtype=float)


def ecef_to_eci(r_ecef: np.ndarray, t_s: float, theta_g0_rad: float = 0.0) -> np.ndarray:
    """Rotate an ECEF position vector (km) into ECI at time t (s). Inverse
    of :func:`eci_to_ecef`."""
    theta = theta_g(t_s, theta_g0_rad)
    R = _rot3(theta)
    return R @ np.asarray(r_ecef, dtype=float)


def ecef_to_geocentric_latlon(r_ecef: np.ndarray) -> tuple[float, float]:
    """ECEF position (km) -> (geocentric latitude, longitude), both
    degrees. Longitude wrapped to [-180, 180)."""
    r_ecef = np.asarray(r_ecef, dtype=float)
    r = np.linalg.norm(r_ecef)
    lat = np.degrees(np.arcsin(r_ecef[2] / r))
    lon = np.degrees(np.arctan2(r_ecef[1], r_ecef[0]))
    lon = ((lon + 180.0) % 360.0) - 180.0
    return lat, lon


def eci_array_to_ecef_array(
    r_eci_array: np.ndarray, t_s_array: np.ndarray, theta_g0_rad: float = 0.0
) -> np.ndarray:
    """Vectorized :func:`eci_to_ecef`: rotate an (N,3) array of ECI
    positions to ECEF at the corresponding (N,) array of times — an M5
    performance primitive (regional-grid coverage needs the satellite
    trajectory rotated at thousands of samples). Builds one batched
    rotation-matrix array via `np.matmul`/`einsum` rather than looping in
    Python. Algebraically identical to calling :func:`eci_to_ecef` once
    per sample (cross-checked in tests)."""
    r_eci_array = np.asarray(r_eci_array, dtype=float)  # (N,3)
    t_s_array = np.asarray(t_s_array, dtype=float)  # (N,)
    theta = (theta_g0_rad + EARTH_ROTATION_RATE_RAD_S * t_s_array) % (2.0 * np.pi)

    c, s = np.cos(-theta), np.sin(-theta)
    zeros, ones = np.zeros_like(theta), np.ones_like(theta)
    R = np.stack(
        [np.stack([c, -s, zeros], axis=-1),
         np.stack([s, c, zeros], axis=-1),
         np.stack([zeros, zeros, ones], axis=-1)],
        axis=-2,
    )  # (N,3,3), each R3(-theta_k)
    return np.einsum("nij,nj->ni", R, r_eci_array)


def geocentric_latlon_to_ecef(
    lat_deg: float, lon_deg: float, r_km: float = R_EARTH
) -> np.ndarray:
    """Spherical (lat, lon, radius) -> ECEF position vector (km). Inverse
    of :func:`ecef_to_geocentric_latlon` (given the same radius)."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return r_km * np.array(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )
