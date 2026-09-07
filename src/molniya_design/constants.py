"""Physical constants and the M1 baseline element set.

All values below reproduce the M1 hand-calculation checkpoint table in
DESIGN.md §9. Units are documented explicitly on every constant since this
module is the single source of truth for every later milestone.

Unit convention for this whole package: **kilometers, seconds, radians**
internally in all numerical code. Degrees are used only at the user-facing
edges (element definitions, documentation, printed reports) and are
converted to radians immediately on entry.

M2 precision note (transparent correction; see DESIGN.md §M2 "Genuine
discrepancy found and fixed")
--------------------------------------------------------------------------
M1's DESIGN.md table reports a, e, ra, etc. rounded to 4-7 significant
figures for human-readable hand-calculation purposes (e.g. e = 0.737286).
Multiplying the *rounded* a and e to recover ra = a(1+e) does not
self-consistently reproduce the rounded ra to the ~1e-7 relative tolerance
needed for automated regression testing (rounding e alone introduces an
~0.01 km error in ra when multiplied through a ~26562 km semi-major axis).

This is a display-precision artifact, not a physics error: the *design
choices* (target period = half sidereal day, perigee altitude = 600 km)
are unchanged from M1. To eliminate this whole class of truncation drift,
this module now *re-derives* a, e, and every dependent quantity directly
from those two design choices at import time, in full float64 precision,
rather than storing hand-rounded decimals. The DESIGN.md M1 table (rounded
for readability) is left untouched as the historical record; the values
below agree with it to the precision DESIGN.md itself reports.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Earth physical constants (WGS-84 / standard values) — identical to DESIGN.md
# ---------------------------------------------------------------------------

MU_EARTH = 398600.4418
"""Earth gravitational parameter GM, km^3/s^2."""

R_EARTH = 6378.137
"""Earth equatorial radius, km."""

J2_EARTH = 1.08262668e-3
"""Earth J2 (dimensionless second zonal harmonic). Not used until M4."""

SIDEREAL_DAY_S = 86164.0905
"""Length of one sidereal day, s."""

# ---------------------------------------------------------------------------
# M1 baseline Molniya orbit: re-derived at full precision from the two
# design choices made in DESIGN.md §0 (target period, perigee altitude).
# ---------------------------------------------------------------------------

_T_TARGET_S = SIDEREAL_DAY_S / 2.0
_PERIGEE_ALT_KM = 600.0
_I_CRIT_DEG = np.degrees(np.arccos(np.sqrt(1.0 / 5.0)))  # critical inclination, DESIGN.md §2

_A_KM = (MU_EARTH * (_T_TARGET_S / (2.0 * np.pi)) ** 2) ** (1.0 / 3.0)
_RP_KM = R_EARTH + _PERIGEE_ALT_KM
_E = 1.0 - _RP_KM / _A_KM
_RA_KM = _A_KM * (1.0 + _E)

BASELINE_ELEMENTS_DEG = {
    "a_km": _A_KM,             # semi-major axis, km
    "e": _E,                   # eccentricity, dimensionless
    "i_deg": _I_CRIT_DEG,      # inclination, deg (critical inclination)
    "raan_deg": 0.0,           # RAAN, deg (conventional epoch, DESIGN.md §0.4)
    "argp_deg": 270.0,         # argument of perigee, deg
    "nu0_deg": 180.0,          # initial true anomaly (apogee), deg
}

# M1 analytical regression targets, now self-consistently re-derived at
# full precision (previously hand-rounded; see module docstring). Compare
# against DESIGN.md §9 to the precision DESIGN.md reports (4-7 sig figs);
# agreement there is exact within that reporting precision.
M1_RP_KM = _RP_KM
M1_RA_KM = _RA_KM
M1_PERIGEE_ALT_KM = _PERIGEE_ALT_KM
M1_APOGEE_ALT_KM = _RA_KM - R_EARTH
M1_PERIOD_S = 2.0 * np.pi * np.sqrt(_A_KM**3 / MU_EARTH)
M1_PERIOD_H = M1_PERIOD_S / 3600.0
M1_MEAN_MOTION_RAD_S = np.sqrt(MU_EARTH / _A_KM**3)
M1_SEMILATUS_RECTUM_KM = _A_KM * (1.0 - _E**2)
M1_VP_KM_S = np.sqrt(MU_EARTH * (2.0 / _RP_KM - 1.0 / _A_KM))
M1_VA_KM_S = np.sqrt(MU_EARTH * (2.0 / _RA_KM - 1.0 / _A_KM))
M1_SPECIFIC_ENERGY_KM2_S2 = -MU_EARTH / (2.0 * _A_KM)
M1_SPECIFIC_ANG_MOMENTUM_KM2_S = np.sqrt(MU_EARTH * M1_SEMILATUS_RECTUM_KM)
