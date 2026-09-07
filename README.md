# Molniya Orbit Design

High-latitude communications coverage via a classical Molniya-style,
critical-inclination, half-sidereal-day orbit. A defensible baseline
orbital design, ground-track and high-latitude dwell/coverage analysis,
sensitivity study, and independent numerical verification — built and
verified milestone by milestone.

## Status

**Milestones 1–2 of 6 — complete.** M1 is the analytical design and
derivation; M2 implements and numerically verifies classical-element ↔
Cartesian-state conversion and two-body propagation against the M1
baseline (period, apsis radii/velocities, conservation, an independent
Kepler time-of-flight cross-check, and the northern-apogee orientation
claim). No Earth-fixed ground track, J2, or coverage code exists yet. See
[`DESIGN.md`](DESIGN.md) for the full derivation, baseline element set,
M2 numerical results, and the verification plan later milestones are held
to.

## Milestone roadmap

| Milestone | Scope | Status |
|---|---|---|
| M1 | Analytical design + verification plan | ✅ complete |
| M2 | Element/state conversion + two-body propagation + apsis/period verification | ✅ complete |
| M3 | ECI/ECEF transformation + Earth-fixed ground track + access geometry | not started |
| M4 | First-order J2 secular propagation + critical-inclination verification | not started |
| M5 | High-latitude dwell/coverage/revisit + parameter sensitivity | not started |
| M6 | Independent validation + portfolio polish + CI/reproducibility audit | not started |

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
src/molniya_design/    Python package (placeholder in M1)
tests/                 pytest suite
scripts/               analysis/plotting scripts (added from M2 onward)
figures/               generated figures (added from M3 onward)
results/               generated numeric results (added from M2 onward)
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
