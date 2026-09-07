# Molniya Orbit Design

High-latitude communications coverage via a classical Molniya-style,
critical-inclination, half-sidereal-day orbit. A defensible baseline
orbital design, ground-track and high-latitude dwell/coverage analysis,
sensitivity study, and independent numerical verification — built and
verified milestone by milestone.

## Status

**Milestone 1 of 6 — complete.** Analytical design and derivation only; no
propagation or coverage code has been implemented yet. See
[`DESIGN.md`](DESIGN.md) for the full derivation, baseline element set, and
verification plan that later milestones are held to.

## Milestone roadmap

| Milestone | Scope |
|---|---|
| M1 | Analytical design + verification plan |
| M2 | Element/state conversion + two-body propagation + apsis/period verification |
| M3 | ECI/ECEF transformation + Earth-fixed ground track + access geometry |
| M4 | First-order J2 secular propagation + critical-inclination verification |
| M5 | High-latitude dwell/coverage/revisit + parameter sensitivity |
| M6 | Independent validation + portfolio polish + CI/reproducibility audit |

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
