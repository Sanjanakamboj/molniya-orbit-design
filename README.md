# Molniya Orbit Design

High-latitude communications coverage via a classical Molniya-style,
critical-inclination, half-sidereal-day orbit. A defensible baseline
orbital design, ground-track and high-latitude dwell/coverage analysis,
sensitivity study, and independent numerical verification — built and
verified milestone by milestone.

## Status

**Milestones 1–4 of 6 — complete.** M1 is the analytical design and
derivation; M2 implements and numerically verifies classical-element ↔
Cartesian-state conversion and two-body propagation; M3 adds the
ECI→ECEF transform, the Earth-fixed two-body ground track, and basic
spherical-Earth topocentric access geometry; M4 adds a first-order
secular J2 mean-element model, verifies the critical-inclination result
numerically, and quantifies J2-driven ground-track drift over 1–14 days.
**No coverage/revisit study, no constellation design, no link budget**
yet — and **M3's access-geometry numbers (80.30% access fraction etc.)
were computed with the two-body-only model; they do not yet include J2
drift**, which M5 will need to account for. See [`DESIGN.md`](DESIGN.md)
for the full derivation and M2/M3/M4 numerical results.

## Milestone roadmap

| Milestone | Scope | Status |
|---|---|---|
| M1 | Analytical design + verification plan | ✅ complete |
| M2 | Element/state conversion + two-body propagation + apsis/period verification | ✅ complete |
| M3 | ECI/ECEF transformation + Earth-fixed ground track + access geometry | ✅ complete |
| M4 | First-order J2 secular propagation + critical-inclination verification | ✅ complete |
| M5 | High-latitude dwell/coverage/revisit + parameter sensitivity | not started |
| M6 | Independent validation + portfolio polish + CI/reproducibility audit | not started |

## M4 verification highlights

- Baseline critical inclination (63.434949°): `argp_dot` = 7.2×10⁻¹⁷
  deg/day (zero to double-precision floor); `RAAN_dot` = **-0.145135
  deg/day**, matching the M1 prediction and an independently coded
  cross-check formula
- `argp_dot` sign flips exactly at critical inclination: +0.0406 deg/day
  at i=60°, **0** at i=63.4349°, -0.0174 deg/day at i=65° — see
  [`figures/m4_j2_rate_sensitivity.png`](figures/m4_j2_rate_sensitivity.png)
- ω(t) stays exactly at 270.000000° over 14 days at the baseline
  inclination, vs. a measurable -0.243° drift at an off-critical i=65°
  case — see [`figures/m4_element_drift.png`](figures/m4_element_drift.png)
- Two-body vs. J2 ground-track divergence grows from -0.167° (1 day) to
  -2.335° (14 days), fully explained by a RAAN-regression + M-dot-timing
  decomposition — see
  [`figures/m4_ground_track_two_body_vs_j2.png`](figures/m4_ground_track_two_body_vs_j2.png)
- Independent Cartesian J2 propagation (osculating, 7 days, perigee-
  sampled fit) confirms the secular RAAN_dot to 0.022% relative error
- J2=0 limit reproduces the M2 two-body dynamics exactly, through two
  independent code paths

Full numbers, the RAAN/M-dot drift decomposition, and one documented
figure-layout fix are in
[`DESIGN.md` — Milestone 4](DESIGN.md#milestone-4--first-order-j2-secular-propagation--critical-inclination-verification).

## M3 verification highlights

- Ground track (two-body, no J2) over one sidereal day; apogee latitude
  **+63.434949°** at every apogee, confirming the M1/M2 northern-apogee
  orientation survives the Earth-fixed transform unchanged
- Successive-apogee longitude separation: exactly **-180.00000000°**
  (since `omega_E * T = pi` by construction) — the classic two-lobe
  Molniya ground-track pattern, quantified numerically
- One-sidereal-day ground-track repeat error: **1.17×10⁻¹⁰ deg**
- Representative site (65°N, 40°E, 10° min elevation) one-day access:
  **80.30% access fraction**, longest pass **8.74 h**, max gap **2.38 h**,
  peak elevation **65.48°**
- Independent triangle-geometry elevation cross-check vs. the ENU-based
  method: max error **1.6×10⁻¹³ deg** (double-precision floor)
- Timestep convergence (120/60/15 s): every access metric shrinks
  monotonically as sampling is refined — see
  [`figures/m3_ground_track.png`](figures/m3_ground_track.png),
  [`figures/m3_access_vs_time.png`](figures/m3_access_vs_time.png),
  [`figures/m3_range_vs_elevation.png`](figures/m3_range_vs_elevation.png)

Full numbers, the pass-count windowing caveat, and one documented
figure-rendering fix (not a numerical bug) are in
[`DESIGN.md` — Milestone 3](DESIGN.md#milestone-3--ecivecef-transformation--ground-track--basic-access-geometry).

## M2 verification highlights

- Numerical (DOP853, rtol=atol=1e-13) period vs. M1 analytical period:
  relative error **1.3×10⁻¹³**
- One-period state closure: relative position error **3.4×10⁻¹³**
- Specific-energy / angular-momentum conservation drift over 2 periods:
  **~1×10⁻¹²** (double-precision floor)
- Independent Kepler-equation time-of-flight cross-check (not a
  `solve_ivp` self-check): max relative error **6.2×10⁻¹³**
- Northern-apogee orientation (ω=270°) numerically confirmed at
  **+63.4349°** latitude; ω=90° control case confirmed **-63.4349°**
  (southern) — see [`figures/m2_orbit_geometry.png`](figures/m2_orbit_geometry.png)
  and [`figures/m2_conservation_and_apsides.png`](figures/m2_conservation_and_apsides.png)

Full numbers, a tolerance-convergence study, and one documented
precision-rounding fix (not a physics bug) are in
[`DESIGN.md` — Milestone 2](DESIGN.md#milestone-2--elementstate-conversion--two-body-propagation-verification).

## Baseline design (M1 summary)

| Element | Value |
|---|---|
| Semi-major axis a | 26561.762 km |
| Eccentricity e | 0.737286 |
| Inclination i | 63.4349° (critical inclination) |
| Argument of perigee ω | 270° |
| Period T | 11.967 h (≈ half sidereal day) |
| Perigee / apogee altitude | 600 km / 39767 km |

Full derivation, justification, hand-calculation checkpoint table, and
verification plan: [`DESIGN.md`](DESIGN.md).

## Repository layout

```
DESIGN.md              analytical design & verification plan
src/molniya_design/    Python package (constants, elements, twobody,
                        propagation, frames, groundtrack, access,
                        j2, j2_cartesian)
tests/                 pytest suite (98 tests: M1 + M2 + M3 + M4)
scripts/               verification-report / figure-generation scripts
figures/               generated figures (M2 orbit/conservation; M3
                        ground track / access / range-vs-elevation;
                        M4 J2 rate sensitivity / ground-track comparison
                        / element drift)
results/               generated numeric verification reports
```

## Development

```bash
pip install -e ".[dev]"
pytest -W error
```

## Limitations (M1)

Point-mass spacecraft, spherical-Earth coverage geometry, two-body baseline
with first-order secular J2 planned for M4 only, no drag/SRP/lunisolar
perturbations, no stationkeeping, no RF link budget, no claim of continuous
single-satellite coverage. Full list in [`DESIGN.md`](DESIGN.md#11-limitations).
