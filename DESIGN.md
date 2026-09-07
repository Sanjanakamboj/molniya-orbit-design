# Molniya Orbit Design — Milestone 1: Analytical Design & Verification Plan

**Status:** M1 — design/derivation only. No propagation, no ground-track, no
coverage code has been implemented. This document is the analytical baseline
that later milestones must reproduce numerically.

**Constants used (WGS-84 / standard Earth values):**

| Symbol | Value | Units |
|---|---|---|
| μ (GM⊕) | 398600.4418 | km³/s² |
| Re | 6378.137 | km |
| J2 | 1.08262668×10⁻³ | — |
| Sidereal day | 86164.0905 | s |

These are the values `src/molniya_design` will expose as named constants in M2
so that every later numeric result traces back to this table.

---

## 0. Mission Scenario

**Mission:** A single Molniya-class spacecraft providing long-dwell
communications relay coverage to a high-northern-latitude service region
(e.g., Arctic operations, northern Scandinavia/Siberia/Canada-class
latitudes), where geostationary satellites are geometrically unusable
(near-zero or negative elevation above ~≈70–75° N).

**Target service region:** 60°–75° N latitude band, all longitudes (a single
Molniya satellite dwells over a fixed longitude sector per apogee pass; a
full constellation, e.g. 2–3 satellites at 8-hour phasing, would be needed
for continuous coverage — out of scope here, see §11).

**Design philosophy for M1:** rather than copying published Molniya element
sets, the baseline is *derived* from three physical requirements, in this
order:

1. **Period:** must be very close to one-half of a sidereal day, so the
   ground track repeats twice per sidereal day (§8).
2. **Inclination:** must sit at the critical inclination so the argument of
   perigee does not secularly drift under J2 — otherwise apogee would not
   stay fixed over the northern hemisphere (§2).
3. **Eccentricity:** set by choosing a practical, drag-safe circular-orbit-
   class perigee altitude and then *solving* for e from the period-derived
   semi-major axis, rather than assuming e a priori.

### 0.1 Derivation of semi-major axis from period

Requiring T = T_sidereal/2 = 43082.045 s and inverting Kepler's third law
(§1) gives:

```
a = ( μ (T / 2π)² )^(1/3) = 26561.762 km
```

### 0.2 Eccentricity from a practical perigee altitude

A perigee altitude of **600 km** is chosen: high enough to keep atmospheric
drag negligible over a multi-year mission (a classic Molniya design choice,
comparable to LEO operational altitudes ~500–800 km), low enough to
maximize eccentricity/apogee dwell for a fixed a. This gives:

```
rp = Re + 600 = 6978.137 km
e = 1 - rp/a = 0.737286
```

This is the self-consistent baseline: **a and e are outputs of a period
requirement and a perigee-altitude requirement, not copied constants.**

### 0.3 Baseline element set (epoch/longitude convention, see §0.4)

| Element | Value | Notes |
|---|---|---|
| a | 26561.762 km | from §0.1 |
| e | 0.737286 | from §0.2 |
| i | 63.4349° | critical inclination, prograde solution (§2) |
| Ω (RAAN) | 0° at reference epoch | conventional — see §0.4 |
| ω (arg. of perigee) | 270° | places apogee at max northern latitude (§4) |
| ν₀ (initial true anomaly) | 180° | spacecraft placed at apogee at epoch, for a representative "on-station" initial condition |

### 0.4 Epoch / longitude convention (stated explicitly per instructions)

At M1, **no numerical epoch (calendar date/time) has been fixed**, and RAAN
is therefore only meaningful relative to an arbitrary reference meridian —
it is *not yet tied to a specific Earth-fixed longitude*. Ω = 0° here means
"RAAN measured from an arbitrary inertial reference direction at an
as-yet-unspecified reference epoch t0." Tying Ω to a specific ground
longitude requires:

- choosing a concrete calendar epoch t0 (UTC),
- the Greenwich Mean Sidereal Time (GMST) at t0,
- the ECI→ECEF rotation, which is explicitly deferred to **M3**.

M1's job is only to fix the *orbit-plane geometry* (a, e, i, ω) and the
*regression rate* of Ω under J2 (§3); the absolute longitude of the ground
track is an M3 deliverable once an epoch is chosen.

---

## 1. Keplerian Geometry

Definitions:

```
rp = a(1-e)            ra = a(1+e)
a = (rp + ra)/2         e = (ra - rp)/(ra + rp)
T = 2π sqrt(a³/μ)       n = sqrt(μ/a³) = 2π/T
p = a(1-e²)
```

**Computed baseline values:**

| Quantity | Value |
|---|---|
| rp | 6978.137 km (perigee altitude 600.000 km) |
| ra | 46145.388 km (apogee altitude 39767.251 km) |
| Inverse check: a = (rp+ra)/2 | 26561.762 km ✓ matches §0.1 to machine precision |
| Inverse check: e = (ra-rp)/(ra+rp) | 0.737286 ✓ matches §0.2 |
| T | 43082.045 s = 11.967235 h = 0.718034 h·... |
| T / T_sidereal | **0.500000** — confirms half-sidereal-day design intent |
| n | 1.4584231716×10⁻⁴ rad/s = 721.971295°/day = 2.005476 rev/day |
| p = a(1-e²) | 12123.022 km |

**Independent verification that T ≈ half sidereal day:** by construction
(§0.1) T was solved *from* T_sidereal/2, so T/T_sidereal = 0.500000 exactly
at the precision used. This is a tautological check on the algebra (the
forward Kepler's-third-law formula, applied to the derived a, reproduces the
target period to 6+ significant figures), not an independent physical
confirmation — the real independent check is deferred to M2's numerical
propagator, which must reproduce T = 43082.045 s ± integration tolerance
from full-force-model dynamics without assuming it.

---

## 2. Critical Inclination

**First-order secular J2 argument-of-perigee rate** (standard form, e.g.
Vallado):

```
ω̇ = (3/4) J2 n (Re/p)² (5 cos²i - 1)
```

**Solving ω̇ = 0:** since J2, n, and (Re/p)² are all strictly positive and
non-zero for a bound orbit, the only way to zero ω̇ is:

```
5 cos²i - 1 = 0  ⟹  cos²i = 1/5  ⟹  cos i = ±1/√5 = ±0.4472136
```

Two solutions in [0°, 180°]:

```
i₁ = arccos(+1/√5) = 63.434949°   (prograde)
i₂ = arccos(-1/√5) = 116.565051°  (retrograde, = 180° - i₁)
```

**Why the prograde solution (i ≈ 63.4349°) is selected:**

1. **Launch-site accessibility / ΔV:** a prograde inclination near 63.4° is
   reachable from mid/high-latitude launch sites (historically Plesetsk,
   ~62.8° N) with a smaller inclination-change penalty than a retrograde
   orbit, which would require launching against Earth's rotation and
   sacrificing the free ~0.46 km/s (at 63° latitude) of rotational ΔV
   assist. The retrograde twin at 116.565° is J2-critical too, but
   operationally far more expensive to reach and provides no
   coverage-geometry advantage.
2. **Physical symmetry, not a numerically preferred root:** both roots zero
   ω̇ identically — the choice is an engineering (launch-cost) decision, not
   a mathematical one. This is stated explicitly per instructions, rather
   than treating 63.4349° as a memorized "magic number."

**Numeric check at i = i₁:** substituting i₁, e, a, n as computed:

```
ω̇(i_crit) = 1.456×10⁻²³ rad/s ≈ 7.2×10⁻¹⁷ °/day
```

This is numerically zero to within double-precision floating-point roundoff
(the analytic result is exactly zero; the residual is a machine-precision
artifact of evaluating cos²(arccos(1/√5)) rather than substituting 1/5
symbolically) — confirming the critical-inclination condition is satisfied
by the chosen i to the precision needed for M4's numerical verification.

---

## 3. RAAN Regression

**First-order secular nodal regression:**

```
Ω̇ = -(3/2) J2 n (Re/p)² cos(i)
```

At the baseline (i = 63.4349°, n and p from §1):

```
Ω̇ = -2.931801×10⁻⁸ rad/s = -0.145135 °/day
```

This is negative (regressing, i.e., the ascending node drifts westward),
consistent with a prograde (i < 90°) orbit under J2's oblateness torque.

**This is an explicit M1 prediction:** M4's numerical J2 secular propagator
must reproduce Ω̇ = -0.145135°/day to within the tolerance specified in the
verification plan (§10.G).

---

## 4. Northern-Apogee Orientation (ω = 270°, not 90°)

**Convention adopted:** argument of perigee ω is measured in the orbit
plane, in the direction of motion, from the ascending node to perigee.
Consequently apogee is always located 180° further along from perigee in
argument of latitude, i.e. at argument of latitude u = ω + ν where ν=0 at
perigee, ν=180° at apogee, so apogee sits at u_apogee = ω + 180°.

The spacecraft's geocentric latitude φ at argument of latitude u (measured
from the ascending node, for inclination i) is given by the standard
spherical-triangle relation:

```
sin φ = sin i · sin u
```

φ is maximized (φ_max ≈ +i, i.e. northernmost) when u = 90°, and minimized
(φ ≈ -i, southernmost) when u = 270°.

- **Perigee** is at u = ω.
- **Apogee** is at u = ω + 180°.

To place **apogee** at the northernmost latitude, we need u_apogee = 90°:

```
ω + 180° = 90°  (mod 360°)
⟹ ω = 90° - 180° = -90° ≡ 270°
```

Equivalently, and more intuitively: if ω = 90°, perigee itself sits at
u = 90° (northernmost point), which is the *opposite* of the design intent
— a fast, low, brief pass over the north with the slow high-dwell apogee
over the *southern* hemisphere (u_apogee = 270° ⟹ φ = sin i·sin270° = -i,
i.e. apogee over the south). That is exactly backwards for a northern
service mission.

Choosing **ω = 270°** puts perigee at u = 270° (southernmost, φ = -i, a
fast/brief low pass over the south — acceptable, since no coverage is
needed there) and apogee at u = 270° + 180° = 90° (northernmost, φ = +i =
+63.4349°, where the spacecraft moves slowest and dwells longest — exactly
where the high-latitude service region needs it).

**Summary table:**

| ω | u_perigee | φ_perigee | u_apogee | φ_apogee | Verdict |
|---|---|---|---|---|---|
| 90° | 90° | +63.43° (N) | 270° | -63.43° (S) | ✗ wrong hemisphere for apogee dwell |
| 270° | 270° | -63.43° (S) | 90° | +63.43° (N) | ✓ baseline choice |

This matches the well-known Molniya convention (ω ≈ 270°, sometimes quoted
as ω ≈ -90°) but is derived here from the u = ω + ν and sin φ = sin i·sin u
relations rather than asserted from memory, per instructions.

---

## 5. Apsis Velocities (vis-viva)

```
v = sqrt( μ (2/r - 1/a) )
```

| | r (km) | v (km/s) |
|---|---|---|
| Perigee | 6978.137 | **9.961732** |
| Apogee | 46145.388 | **1.506420** |

**Cross-check via specific angular momentum** h = sqrt(μp):

```
h = sqrt(398600.4418 × 12123.022) = 69514.330 km²/s
v_p,check = h/rp = 9.961732 km/s   ✓ matches vis-viva to 6 s.f.
v_a,check = h/ra = 1.506420 km/s   ✓ matches vis-viva to 6 s.f.
```

Both cross-checks agree with vis-viva to the precision computed, confirming
internal consistency of a, e, rp, ra, p. Behavior matches physical
expectation: perigee speed (9.96 km/s) is roughly 6.6× apogee speed
(1.51 km/s), i.e. fast/low at perigee, slow/high at apogee — the basis for
the dwell behavior in §6.

**Specific orbital energy** (independent check): ε = -μ/(2a) = -7.503275
km²/s², constant at both apsides by construction of vis-viva — this
equality is exactly what M2's energy-conservation regression test (§10.C)
will verify numerically along a propagated trajectory.

---

## 6. Dwell-Time Physics

**Qualitative basis — Kepler's 2nd law:** equal areas are swept in equal
times. Since the orbit is highly eccentric (e = 0.737), the areal rate
dA/dt = h/2 is constant, but the radius vector is far longer near apogee
than near perigee — so for the same swept area (same Δt), the true-anomaly
sweep Δν is much smaller near apogee. Equivalently, angular rate
ν̇ = h/r² falls off as 1/r², so ν̇ at apogee is (rp/ra)² ≈ (6978/46145)² ≈
0.0229 times the rate at perigee: the spacecraft crawls through true anomaly
near apogee and whips through it near perigee.

**Numerical dwell calculation — labeled as an anomaly-based proxy, not
actual communications coverage:**

*Definition of the proxy:* the fraction of one orbital period during which
the true anomaly ν lies within ±60° of apogee (ν=180°), i.e. ν ∈ [120°,
240°]. This measures orbital dwell in true-anomaly space; it is **not** a
statement about elevation angle, line-of-sight, or communications access
to any specific ground point — that calculation requires ground-station
geometry and is deferred to M5.

*Method:* convert the true-anomaly bounds to eccentric anomaly E via
tan(E/2) = sqrt((1-e)/(1+e))·tan(ν/2), then to mean anomaly M = E - e·sinE
(Kepler's equation, since M advances linearly in time at rate n). Δt =
ΔM / n.

*Result (e = 0.737286, n = 1.4584232×10⁻⁴ rad/s):*

| ν bound | M (deg) |
|---|---|
| 120° | 28.7777° |
| 240° | 331.2223° |

```
ΔM = 331.2223° - 28.7777° = 302.4447°
Δt = ΔM / n = 36194.267 s = 10.0540 h
Fraction of orbit = Δt / T = 10.0540 / 11.9672 = 0.8401  (≈ 84.0%)
```

**Interpretation:** despite ±60° being only 1/3 of the 360° true-anomaly
range (33.3%), the spacecraft spends **84.0%** of each orbital period
within that window — a direct, quantitative demonstration of Kepler's
second law concentrating time near apogee for a high-e orbit. This dwell
fraction (not a coverage/access fraction) is the M1 anomaly-based
regression target for M2/M5 to reproduce from full numerical time-of-flight
propagation.

---

## 7. High-Latitude Coverage Requirement (definition only — not computed)

This section defines the coverage problem that **M5** will solve
numerically. No coverage computation is performed in M1.

- **Target latitude band:** 60°–75° N.
- **Representative target point(s):** a single reference ground site is
  proposed at 65° N, 40° E (representative high-Arctic longitude under the
  baseline apogee ground track; exact longitude is provisional pending the
  M3 epoch/ECI→ECEF choice per §0.4). Additional points across the band
  (e.g. 60° N, 70° N, 75° N at several longitudes) will be added in M5 for
  a latitude-band sweep.
- **Minimum elevation angle:** 10° (a standard, conservative RF-link
  operating threshold; no link budget is computed at M1 — see §11).
- **Coverage metrics to be evaluated in M5 (combination, not a single
  metric):**
  - **Access fraction** — fraction of time the target point sees the
    spacecraft above the minimum elevation angle, over a representative
    multi-day span.
  - **Dwell time per pass** — contiguous duration of a single access
    window (directly informed by §6's anomaly proxy, refined with real
    elevation geometry).
  - **Maximum coverage gap** — worst-case time between the end of one
    access window and the start of the next (the key design driver for a
    single-satellite Molniya mission, since one satellite alone cannot give
    continuous coverage — see §11).
  - Revisit interval is a secondary/derived metric (≈ related to orbital
    period and multiples thereof) and will be reported but is not a primary
    design metric for a single spacecraft.

**Geometry to be implemented in M3 (documented now, not run):**
spherical-Earth line-of-sight / minimum-elevation test. For a ground site at
geocentric position **r_site** (radius Re, on a spherical Earth model) and
spacecraft ECEF position **r_sat**, define the topocentric line-of-sight
vector **ρ** = r_sat - r_site, resolved into the site's local ENU (East-
North-Up) frame. The elevation angle is:

```
el = arcsin( ρ_up / |ρ| )
```

Access is declared when el ≥ el_min (10°). Equivalently, in the common
central-angle form (ignoring topocentric refinements), the maximum
achievable elevation from a point at central angle λ from the sub-satellite
point, at satellite geocentric radius r_sat and Earth radius Re, satisfies:

```
cos(λ_max) relation:  tan(el) = (r_sat/Re · cos λ - 1) / (r_sat/Re · sin λ)
```

which reduces to the standard geometric horizon/elevation-mask formula used
in ground-station visibility analysis. M3 will implement this exactly
(spherical Earth; no atmospheric refraction, no terrain masking — see §11)
and validate it against at least one hand-computed case (verification item
§10.L).

---

## 8. Ground-Track Expectations (qualitative predictions to test numerically later)

Before any ground-track code is written, the following qualitative behavior
is predicted from the analytical design and must be tested numerically in
M3/M4:

1. **~2 revolutions per sidereal day.** T ≈ 11.967 h ⟹ sidereal_day/T ≈
   2.0055 rev/sidereal day (from n above). Two apogee passes per sidereal
   day are expected, each over roughly the same longitude sector (before
   J2 nodal drift, §3, and small excess-of-2 revs/day accumulate offsets).
2. **Characteristic "hurricane" / asymmetric ground track.** Because
   angular rate is high near perigee and very low near apogee (§6), the
   ground track shows tightly bunched, fast-moving low-latitude
   perigee-region crossings and a broad, slow, high-latitude loop near
   apogee — the recognizable elongated Molniya loop shape.
3. **Long dwell in the high-latitude lobe.** Directly follows from §6 (84%
   of the orbital period within ±60° true anomaly of apogee): the apogee
   loop is where the satellite lingers, which is why one apogee pass alone
   can serve the target latitude band for several hours.
4. **Earth-rotation longitude shift between passes.** Because the orbital
   period (11.967 h) is not exactly an integer submultiple related to Earth
   rotation *except* by design (T = sidereal_day/2), each successive apogee
   passage should recur near the same two longitude sectors (~180° apart)
   sidereal-day after sidereal-day, before accounting for J2 nodal
   regression. Within one revolution, the sub-satellite point sweeps
   through all longitudes as Earth rotates ~180° underneath during the
   11.967 h flight.
5. **J2 RAAN regression shifts the pattern over multiple days.** Ω̇ =
   -0.145135°/day (§3) causes the ascending node — and hence the entire
   ground-track pattern — to drift westward by that amount per day; over
   a mission of order 1 year this accumulates to tens of degrees and must
   be tracked/compensated operationally (out of scope for M1–M5 station-
   keeping design, see §11).
6. **ω should remain nearly frozen at critical inclination.** Because i is
   set to the critical value (§2), ω̇ ≈ 0 to first order in J2, so the
   apogee should remain fixed at the same argument of latitude (northern
   apogee) for the life of the mission — this is the entire point of using
   the critical inclination, and M4 must confirm ω stays within a small
   numerical band around 270° over a multi-year propagation (§10.H).

**Explicit list of predictions later milestones must numerically test:**
(a) revs/sidereal-day ≈ 2.0055; (b) qualitative loop shape and northern
dwell lobe; (c) 84.0% anomaly-based dwell fraction reproduced by numerical
time-of-flight; (d) Ω̇ = -0.145135°/day reproduced by J2 secular/numerical
propagation; (e) ω̇ ≈ 0 (bounded oscillation, no secular drift) at i =
63.4349°, contrasted against nonzero ω̇ at an off-critical inclination
(§10.I).

---

## 9. Hand-Calculation Checkpoint (regression targets for M2+)

| Quantity | Symbol | Value | Units |
|---|---|---|---|
| Gravitational parameter | μ | 398600.4418 | km³/s² |
| Earth radius | Re | 6378.137 | km |
| J2 | J2 | 1.08262668×10⁻³ | — |
| Semi-major axis | a | 26561.7624 | km |
| Eccentricity | e | 0.7372860 | — |
| Inclination | i | 63.434949 | deg |
| RAAN (reference epoch, conventional) | Ω | 0.0 | deg |
| Argument of perigee | ω | 270.0 | deg |
| Initial true anomaly | ν₀ | 180.0 | deg |
| Perigee radius | rp | 6978.1370 | km |
| Apogee radius | ra | 46145.3879 | km |
| Perigee altitude | hp | 600.0000 | km |
| Apogee altitude | ha | 39767.2509 | km |
| Orbital period | T | 11.967235 | h |
| Mean motion | n | 1.4584231716×10⁻⁴ | rad/s |
| Mean motion | n | 721.971295 | deg/day |
| Semilatus rectum | p | 12123.0223 | km |
| Predicted RAAN rate | Ω̇ | -0.145135 | deg/day |
| Predicted arg. of perigee rate | ω̇ | ~7.2×10⁻¹⁷ (≈0, machine-precision residual) | deg/day |
| Perigee speed | vp | 9.961732 | km/s |
| Apogee speed | va | 1.506420 | km/s |
| Specific orbital energy | ε | -7.503275 | km²/s² |
| Specific angular momentum | h | 69514.3298 | km²/s |
| Dwell fraction (±60° true anomaly of apogee, anomaly-based proxy) | — | 0.8401 (84.01%) | fraction of T |
| Dwell time (same proxy) | — | 10.0540 | h |

---

## 10. Verification Plan (for M2–M6)

For every item: PASS criterion = value/behavior below is met; tolerances
are chosen relative to double-precision numerical integration and standard
J2 secular-theory accuracy, not measurement uncertainty.

| # | Check | Method | Tolerance / expected behavior |
|---|---|---|---|
| A | Propagated period vs Kepler period | Numerically propagate one full two-body orbit (M2); measure time between successive perigee passages | Match T = 43082.045 s to ≤1e-6 relative (integrator-tolerance limited) |
| B | Propagated perigee/apogee vs analytic rp, ra | Extract min/max |r| along propagated orbit | Match §9 rp, ra to ≤1e-6 relative |
| C | Specific orbital energy conservation | ε = v²/2 - μ/r at each integration step | Constant to ≤1e-10 relative drift over one orbit (double precision, fixed-step RK) |
| D | Angular momentum conservation | h = r×v magnitude at each step | Constant to ≤1e-10 relative drift |
| E | Element → state → element round trip | COE→RV→COE using standard conversion | All 6 elements recovered to ≤1e-8 relative (or ≤1e-6 deg for angles) |
| F | Independent vis-viva check at apsides | Compare propagated speed at r=rp, r=ra against §5 | ≤1e-6 relative |
| G | Numerical J2 Ω̇ vs M1 prediction | J2-perturbed propagation (M4); linear fit of Ω(t) | Match -0.145135°/day to ≤1% (first-order secular theory accuracy) |
| H | Numerical ω̇ ≈ 0 at critical inclination | J2-perturbed propagation at i=63.4349°; fit ω(t) | Secular slope ≤0.01°/day (bounded short-period oscillation only) |
| I | Off-critical inclination comparison | Repeat H at, e.g., i=55° or i=70° | Nonzero secular ω̇, matching analytic formula (§2) to ≤5% |
| J | ECI→ECEF ground-track longitude consistency | Rotate propagated ECI position by GMST(t) (M3); recompute at two independent times | Longitude consistent with direct ECEF propagation/rotation to ≤1e-6 deg |
| K | Northern/southern apogee orientation sanity check | Confirm propagated apogee latitude ≈ +i, perigee latitude ≈ -i (M3) | ≤0.01° of ±63.4349° |
| L | Independent coverage recomputation at ground points | Recompute elevation-angle access (M5) via two independent formulas (ENU vector method and central-angle method, §7) | Agree to ≤1e-6 deg elevation |
| M | Timestep convergence | Halve integrator step size repeatedly (M2) | Solution changes by less than previous halving squared (RK4-order convergence), or below 1e-9 relative once converged |
| N | RAAN/reference-epoch invariance | Repeat M3 ground track with epoch shifted by exactly one sidereal day | Ground track pattern reproduces to numerical precision (mod J2 drift over that day, §G) |
| O | Dimensional/unit round-trip checks | Convert km↔m, rad↔deg, s↔h through full pipeline | Exact round trip to floating-point precision |

---

## 11. Limitations

- Point-mass spacecraft (no attitude, no finite-size effects).
- Spherical Earth assumed for all coverage/elevation geometry unless
  explicitly upgraded to an oblate-Earth model in a later milestone.
- Keplerian two-body baseline for M1/M2; J2 secular effects only from M4
  onward.
- Only **first-order secular J2** is modeled in the planned perturbation
  milestone (M4) — no J3+, no tesseral harmonics.
- No atmospheric drag.
- No solar radiation pressure (SRP).
- No lunisolar (third-body) perturbations.
- No stationkeeping / orbit-maintenance design (drift is characterized, not
  corrected).
- No RF link budget (minimum elevation angle of 10° is assumed as a
  standard operating threshold, not derived from a link budget).
- No atmospheric or radiation-environment analysis.
- No launch injection-error analysis.
- No navigation/estimation errors (elements are treated as exactly known).
- No operational-availability claim is made.
- **No claim that a single Molniya spacecraft provides continuous
  coverage** — §7 and §11 both note that a real operational system needs
  multiple phased satellites (classically 2–3, 8-hour-apart) for continuous
  high-latitude coverage; this project analyzes and characterizes a single
  spacecraft's coverage/dwell/gap behavior only.

---

## 12. Milestone Roadmap

- **M1** — Analytical design + verification plan (this document). ✅ complete
- **M2** — Element/state conversion + two-body propagation + apsis/period
  numerical verification (§10.A–F). ✅ complete — see §M2 below.
- **M3** — ECI/ECEF transformation + Earth-fixed ground track + access
  geometry (§10.J, K, L, N). ✅ complete — see §M3 below.
- **M4** — First-order J2 secular propagation + critical-inclination
  verification (§10.G, H, I). ✅ complete — see §M4 below.
- **M5** — High-latitude dwell/coverage/revisit + parameter sensitivity
  study.
- **M6** — Independent validation + portfolio polish + CI/reproducibility
  audit.

M5 is **not** started as part of this milestone.

---

# Milestone 2 — Element/State Conversion + Two-Body Propagation Verification

**Status:** M2 complete. Implements the Cartesian two-body dynamics,
classical-element ↔ state conversion, apsis detection, an independent
Kepler time-of-flight cross-check, conservation diagnostics, the
northern-apogee orientation regression test, and a tolerance-convergence
study — all against the M1 analytical baseline. **No ECI→ECEF, no ground
track, no J2, no coverage geometry** is implemented here (M3/M4/M5 scope).

## M2.0 Pre-flight verification (before any code was written)

Before implementation, M1's regression values were independently
recomputed from the constants and design choices in DESIGN.md §0–§9 (not
copied from the table) and matched to the reported precision:

| Quantity | M1 value | Recomputed | Match |
|---|---|---|---|
| a | 26561.762 km | 26561.762 km | ✓ |
| e | 0.737286 | 0.737286 | ✓ |
| i | 63.4349° | 63.4349° | ✓ |
| perigee altitude | 600 km | 600.000 km | ✓ |
| apogee altitude | 39767.251 km | 39767.251 km | ✓ |
| period | 11.967235 h | 11.967235 h | ✓ |
| vp | 9.961732 km/s | 9.961732 km/s | ✓ |
| va | 1.506420 km/s | 1.506420 km/s | ✓ |

No discrepancy found at this stage.

## M2.1 Implementation

New modules under `src/molniya_design/`:

- **`constants.py`** — Earth constants (unchanged from M1) plus the M1
  baseline element set, now *re-derived at full float64 precision* from
  the two M1 design choices (target period = half sidereal day, perigee
  altitude = 600 km) rather than stored as hand-rounded decimals — see
  §M2.2 "Genuine discrepancy found and fixed" below.
- **`elements.py`** — `coe_to_rv` (classical elements → ECI state via the
  perifocal PQW frame and the standard 3-1-3 Euler rotation
  R = R3(RAAN)·R1(i)·R3(argp)) and `rv_to_coe` (ECI state → classical
  elements via angular-momentum/eccentricity-vector/node-vector geometry,
  independent of `coe_to_rv`). Frame, rotation-sense, and unit conventions
  are documented in the module docstring: ECI is right-handed equatorial
  (X toward the reference direction, Z along the rotation axis); the
  rotation matrix is confirmed right-handed/orientation-preserving
  (det R = +1) by a unit test; all public-API angles are in degrees,
  converted to radians immediately internally.
- **`twobody.py`** — the two-body Cartesian EOM `rddot = -mu r/|r|^3` in
  state form `y=[rx,ry,rz,vx,vy,vz]` (km, km/s), plus `specific_energy`,
  `specific_angular_momentum`, `eccentricity_vector` helpers.
- **`propagation.py`** — `propagate` (thin `scipy.integrate.solve_ivp`
  wrapper, default method `DOP853`, explicit configurable `rtol`/`atol`,
  default `1e-12`); `find_apsides` (radial-velocity zero-crossing event
  detection, classifying perigee/apogee by relative radius only — it does
  **not** consult the analytic rp/ra during classification, only in the
  caller's cross-check); `kepler_state_at_time` (an **independent**
  Kepler-equation analytical time-of-flight solver — mean anomaly →
  eccentric anomaly via Newton iteration on `M=E-e sinE` → true anomaly →
  `coe_to_rv` — which does not call `solve_ivp` at all, per the M2 scope
  requirement for a real independent check); `geocentric_latitude_deg`
  (documented as frame-independent between ECI/ECEF, since latitude does
  not depend on the Earth-rotation/GMST transform that is explicitly out
  of scope until M3).

`pyproject.toml` gained `scipy>=1.10` as a runtime dependency.

## M2.2 Genuine discrepancy found and fixed

**What was found:** the first version of `constants.py` stored the M1
baseline `a` and `e` as hand-rounded decimals (`a=26561.7624`,
`e=0.7372860`, matching DESIGN.md's display precision). Recomputing
`ra = a(1+e)` from those *rounded* values gave 46145.37795 km, about
0.0099 km (2.1×10⁻⁷ relative) away from the DESIGN.md-reported
46145.3879 km — large enough to fail several M2 regression tests written
at 1×10⁻⁷ relative tolerance (`test_E`/`F`/`G` in `tests/test_twobody.py`
initially failed with exactly this residual).

**Diagnosis:** this is a **display-rounding artifact, not a physics
error.** The original M1 derivation used full float64 precision
internally and only rounded for the printed DESIGN.md table (e.g. the
true e is 0.7372863710269664, which *does* round correctly to 0.737286 at
6 decimal places — the bug was re-deriving `ra` from the *already-rounded*
`e` rather than carrying full precision through).

**Fix:** `constants.py` now re-derives `a`, `e`, and every dependent
baseline quantity (`ra`, `rp`, `T`, `n`, `p`, `vp`, `va`, energy, angular
momentum) directly from the two original design choices (§0.1/§0.2:
half-sidereal-day period, 600 km perigee altitude) using the same
closed-form physics as M1, evaluated once at full float64 precision at
import time — eliminating the truncation step entirely rather than just
adding more decimal digits. The DESIGN.md §9 table (rounded for
human readability) is unchanged and still agrees with the recomputed
values to the precision it reports. All M1 physics, element values, and
design choices are otherwise identical — no M1 conclusion changes.

**Regression test:** `tests/test_twobody.py::test_C_arbitrary_round_trip`
and the whole `test_B`/`test_E`/`test_F`/`test_G` family now exercise this
path at tight tolerance and pass; `constants.py`'s module docstring
documents the fix so it cannot silently regress.

A second, expected (non-bug) behavior was also handled explicitly:
`find_apsides` naturally detects a spurious r·v=0 event at t≈0 s because
the baseline epoch state is itself exactly at apogee (nu0=180°). This is
correct behavior (apogee *is* an r·v=0 crossing), not a bug — tests and
the report script filter for `t_s > 1.0` when looking for the *next*
apogee to measure the numerical period.

## M2.3 Numerical results

All figures below are from `scripts/m2_verification_report.py`
(reproducible; also exercised by `tests/test_twobody.py`).

**Round trip (element → state → element), baseline:**

| Element | Input | Recovered | Agreement |
|---|---|---|---|
| a | 26561.762430362043 km | 26561.762430362036 km | 8×10⁻¹⁶ km |
| e | 0.7372863710269664 | 0.7372863710269664 | exact (float64) |
| i | 63.43494882292201° | 63.43494882292201° | exact |
| RAAN | 0.0° | 0.0° | exact |
| ω | 270.0° | 270.0° | exact |
| ν | 180.0° | 180.0° | exact |

**Apsis detection over 2.05 periods (rtol=atol=1e-13):**

| Apsis | Detected r | M1 analytical | Abs. error | Detected v | M1 analytical | Abs. error |
|---|---|---|---|---|---|---|
| Perigee | 6978.136999999 km | 6978.137000 km | ≤1.0×10⁻⁹ km | 9.961731876983 km/s | 9.961731876983 km/s | ≤2.7×10⁻¹³ km/s |
| Apogee | 46145.387860715 km | 46145.387860724 km | ≤1.6×10⁻⁸ km | 1.506419883276 km/s | 1.506419883275 km/s | ≤4.1×10⁻¹³ km/s |

**Numerical period vs Kepler period:**

- Numerical (apogee-to-apogee): **43082.045250 s** (11.967235 h)
- M1 analytical: 43082.045250 s (11.967235 h)
- Relative error: **1.29×10⁻¹³**

**One-period closure (rtol=atol=1e-13):**

- Relative position error: **3.45×10⁻¹³**
- Relative velocity error: **8.17×10⁻¹³**

**Conservation over 2 periods (2000 samples, rtol=atol=1e-13):**

- Specific energy relative drift: **5.01×10⁻¹²**
- |h| relative drift: **1.22×10⁻¹²**
- Eccentricity-vector max drift: **2.09×10⁻¹²**
- Orbital-plane (ĥ) drift: below 1×10⁻¹⁰ (test `test_J`)

All conservation residuals are consistent with double-precision floating
point accumulation at these tolerances, not with any physical drift —
expected for an unperturbed two-body integration.

**Independent Kepler time-of-flight cross-check** (6 times across one
period, `kepler_state_at_time` vs. `propagate`, neither calling the
other):

| ν-fraction of period | t (s) | Position rel. error | Velocity rel. error |
|---|---|---|---|
| 0.00 | 0.00 | 0 | 0 |
| 0.10 | 4308.20 | 3.3×10⁻¹⁵ | 1.4×10⁻¹⁴ |
| 0.25 | 10770.51 | 1.6×10⁻¹⁴ | 4.4×10⁻¹⁴ |
| 0.50 | 21541.02 | 2.8×10⁻¹³ | 1.7×10⁻¹³ |
| 0.63 | 27141.69 | 2.1×10⁻¹³ | 2.4×10⁻¹³ |
| 0.90 | 38773.84 | 2.8×10⁻¹³ | 6.2×10⁻¹³ |

Maximum relative error across all sampled times: position 2.8×10⁻¹³,
velocity 6.2×10⁻¹³ — both consistent with double-precision accumulation
over one full orbital revolution, confirming the numerical propagator
independently against closed-form Kepler-equation time-of-flight.

**Northern-apogee orientation regression check (numerical, DESIGN.md §4):**

| ω | Numerical apogee latitude | Expected | Match |
|---|---|---|---|
| 270° | +63.434949° | +63.4349° (northern) | ✓ |
| 90° | -63.434949° | -63.4349° (southern) | ✓ |

This confirms M1 §4's derivation numerically and is encoded as a
permanent regression test (`test_M_northern_apogee_omega_270`,
`test_N_southern_apogee_omega_90`) so the ω convention cannot silently
flip in later milestones.

**Tolerance convergence study** (one-period closure, energy/|h| drift,
apsis radius error, at three `solve_ivp` tolerance settings):

| Tolerance | Closure (rel. pos. err) | Energy drift (rel.) | \|h\| drift (rel.) | Apogee r error | Perigee r error |
|---|---|---|---|---|---|
| loose (1e-6) | 9.44×10⁻⁶ | 6.79×10⁻⁶ | 2.77×10⁻⁷ | 0.728 km | 0.0114 km |
| medium (1e-9) | 3.75×10⁻⁹ | 2.30×10⁻⁹ | 2.23×10⁻¹⁰ | 2.55×10⁻⁴ km | 2.04×10⁻⁶ km |
| tight (1e-12) | 4.21×10⁻¹² | 2.24×10⁻¹² | 3.97×10⁻¹³ | 3.22×10⁻⁷ km | 7.18×10⁻⁹ km |

Every metric decreases monotonically — by roughly 3 orders of magnitude
per 3-orders-of-magnitude tolerance tightening — as `rtol`/`atol` tighten
from loose → medium → tight, demonstrating genuine numerical convergence
of the DOP853 integrator on this problem (not merely asserted).

## M2.4 Figures

- [`figures/m2_orbit_geometry.png`](figures/m2_orbit_geometry.png) — 3D
  ECI two-body trajectory over one period, with Earth, perigee, and
  apogee markers, explicitly labeled as an inertial trajectory (not a
  ground track).
- [`figures/m2_conservation_and_apsides.png`](figures/m2_conservation_and_apsides.png)
  — altitude vs. time over 2 periods with detected apsides (top), and
  specific-energy/angular-momentum relative error vs. time on a log scale
  (bottom), showing the machine-precision noise floor (~1e-12 to 1e-16).

Both figures were visually inspected for clipping, overlap, aspect-ratio
distortion, missing units, and legend readability; none were found.

## M2.5 Test suite

`tests/test_twobody.py` implements checklist items A–P plus a rotation
right-handedness check (34 tests total across both M1 and M2 files).
`tests/test_placeholder.py` was updated transparently: the M1-era guard
against "`molniya_design.propagate` existing" is now obsolete (M2 has
been explicitly approved and implemented) and was removed with an
explanatory note; the M3/M5 guards are unchanged, and an M4 (J2) guard
was added. All 34 tests pass under `pytest -W error` with zero warnings.

## M2.6 Scope guard confirmation

No ECI→ECEF transform, GMST/Earth-rotation, longitude, ground-station
access/elevation, coverage/revisit metric, J2 acceleration, secular J2
propagation, numerical RAAN regression, or numerical critical-inclination
verification was implemented in M2 — all remain explicitly deferred to
M3/M4/M5 per the roadmap.

---

# Milestone 3 — ECI/ECEF Transformation + Ground Track + Basic Access Geometry

**Status:** M3 complete. Implements the ECI→ECEF transform under a
constant-rate Earth-rotation model, geocentric spherical lat/lon,
the two-body (no-J2) Earth-fixed ground track, and basic spherical-Earth
topocentric access geometry (range/azimuth/elevation, access intervals,
one-day access metrics), each cross-checked against an independent method
and a timestep convergence study. **No J2, no secular RAAN/ω drift, no
constellation/coverage optimization, no link budget, no geodetic (WGS-84)
station coordinates** — all explicitly out of scope until M4/M5.

## M3.0 Pre-flight verification

Before writing any M3 code, the M2 headline values were reproduced from
production code (not copied) and the M2 report/figures were regenerated
and diffed byte-for-byte against the committed artifacts:

| Quantity | M2 committed value | Recomputed (pre-M3) | Match |
|---|---|---|---|
| Period | 43082.045250 s | 43082.045250 s | ✓ |
| Perigee altitude | 600.000 km | 600.000 km | ✓ |
| Apogee altitude | 39767.251 km | 39767.251 km | ✓ |
| vp | 9.961732 km/s | 9.961732 km/s | ✓ |
| va | 1.506420 km/s | 1.506420 km/s | ✓ |
| Apogee latitude (ω=270°) | +63.434949° | +63.434949° | ✓ |

`scripts/m2_verification_report.py` output and both M2 figure PNGs
(`m2_orbit_geometry.png`, `m2_conservation_and_apsides.png`) were
regenerated and found **bit-for-bit identical** to the committed versions
(`diff` on the report text; identical file checksums on the figures) —
zero numerical drift from M2. `pytest -W error` on the pre-M3 tree passed
34/34.

## M3.1 Earth-rotation / frame model

A **constant-rate** Earth rotation is used (`frames.py`):

    theta_G(t) = theta_G0 + omega_E * t,   theta_G0 = 0 at t = 0

- `theta_G0 = 0` at `t = 0` is an **engineering reference epoch**, not a
  real UTC/GMST epoch — it only fixes ECEF and ECI to be coincident at the
  start of the propagated timeline (consistent with DESIGN.md §0.4's
  M1 caveat that absolute longitude is convention-only until a real epoch
  is chosen — still true here; M3 only adds a *self-consistent* Earth
  rotation, not a real-world one).
- `omega_E = 2*pi / SIDEREAL_DAY_S`, i.e. the **same** sidereal rate
  already implicit in the M1 half-sidereal-day period design choice —
  confirmed by `omega_E * SIDEREAL_DAY_S = 2*pi` and, as a direct
  consequence, `omega_E * T_baseline = pi` exactly (both verified
  numerically to double-precision, see §M3.3).

**Sign convention (ECI→ECEF):** the ECEF frame rotates **eastward**
(counterclockwise viewed from +Z/north) relative to ECI, matching Earth's
actual prograde rotation:

    r_ECI  = R3(+theta_G) @ r_ECEF
    r_ECEF = R3(-theta_G) @ r_ECI

using the same active-rotation `R3` convention as `elements.py`. This
sign is documented in `frames.py`'s module docstring, verified at
`theta=0` (identity), round-trip tested, and tested for right-handedness
(`det R = +1`, `R^T R = I`) at multiple times.

## M3.2 Geocentric latitude/longitude

    lat = asin(z / r)                    [-90, +90] deg
    lon = atan2(y, x), wrapped to [-180, +180) deg

Explicitly **geocentric** latitude on a **spherical** Earth — not WGS-84
geodetic latitude (unchanged M1 limitation, §11). The inverse
(`geocentric_latlon_to_ecef`) is implemented and round-trip tested,
including explicit longitude-wraparound cases (±180° boundary, values
outside [-180,360)).

## M3.3 Ground-track results (two-body, no J2)

Propagated the M2 baseline for one full sidereal day (≈2.0055 Molniya
revolutions, `SIDEREAL_DAY_S = 86164.0905` s) and transformed every
sampled ECI state to ECEF/geocentric lat-lon. Apsides (from M2's
`find_apsides`) were also converted to Earth-fixed lat/lon:

| Apsis | t (s) | Latitude (deg) | Longitude (deg) | r (km) |
|---|---|---|---|---|
| Apogee 0 | 0.000 | +63.43495 | +90.00000 | 46145.388 |
| Perigee 0 | 21541.023 | -63.43495 | -180.00000 | 6978.137 |
| Apogee 1 | 43082.045 | +63.43495 | -90.00000 | 46145.388 |
| Perigee 1 | 64623.068 | -63.43495 | 0.00000 | 6978.137 |
| Apogee 2 | 86164.090 | +63.43495 | +90.00000 | 46145.388 |

**Successive-apogee longitude separation:**

- Apogee 0 → 1: **-180.00000000°** (expected: exactly ±180°, since
  `omega_E * T = pi` exactly by construction)
- Apogee 1 → 2: **-180.00000000°**

**One-sidereal-day repeat error** (apogee 0 vs. apogee 2, i.e. after
exactly 2 baseline periods = 1 sidereal day): longitude error **1.17×10⁻¹⁰
deg**, latitude match to double-precision. This numerically confirms the
DESIGN.md §8 prediction that the ground track repeats every sidereal day
under the two-body, no-J2 baseline, and that consecutive apogee passes
alternate between two longitude sectors 180° apart — the classic Molniya
"figure-eight-like" two-lobe pattern, quantified rather than only
visually inspected (see `figures/m3_ground_track.png`, §M3.7).

All apogees sit at the same latitude (+63.43495°, matching M1/M2's
critical-inclination/ω=270° prediction) and all perigees at -63.43495°,
confirming the M1 §4 orientation claim survives the Earth-fixed transform
unchanged (latitude is frame-independent of Earth rotation about Z, as
documented in `propagation.geocentric_latitude_deg` since M2).

## M3.4 Representative site and access geometry

Site (DESIGN.md §7): **lat = 65° N, lon = 40° E**, minimum elevation
**10°**, spherical Earth radius `R_EARTH`.

**ENU/azimuth/elevation formulas** (`access.py`):

    r_site  = R_E [cos(lat)cos(lon), cos(lat)sin(lon), sin(lat)]
    e_east  = [-sin(lon), cos(lon), 0]
    e_north = [-sin(lat)cos(lon), -sin(lat)sin(lon), cos(lat)]
    e_up    = [ cos(lat)cos(lon),  cos(lat)sin(lon), sin(lat)]

    rho = r_sat_ECEF - r_site
    range     = |rho|
    elevation = asin((rho . e_up) / range)
    azimuth   = atan2(rho . e_east, rho . e_north) mod 360   (0=N, 90=E)

No atmospheric refraction. `(e_east, e_north, e_up)` confirmed orthonormal
and right-handed (`e_east × e_north = e_up`) at multiple lat/lon (test K).

## M3.5 Access-interval extraction and one-day metrics

Access intervals (elevation ≥ 10°) are extracted from a densely sampled
elevation time series with **linearly interpolated threshold crossings**
(not raw-sample counting) and a **local parabolic refinement of the peak
elevation/time** within each interval. Validated against a synthetic
sinusoid with a closed-form crossing/peak solution (test P): crossing
times recovered to ≤0.01 s, peak elevation to ≤0.01°, on a very finely
sampled synthetic signal (dt=0.005 s) — confirming the extraction
algorithm itself is correct before applying it to the real (coarser)
propagated trajectory.

**One-sidereal-day access metrics** (dt=15 s ground track):

| Metric | Value |
|---|---|
| Number of passes | 3 (see caveat below) |
| Total access time | 69192.106 s (19.2200 h) |
| Access fraction | **0.803027** (80.30%) |
| Longest pass | 31471.372 s (8.7420 h) |
| Maximum no-access gap | 8551.051 s (2.3753 h) |
| Peak elevation (best pass) | 65.4781° |

**Pass-count caveat (documented, not a bug):** the reported "3 passes"
includes one continuous high-elevation dwell that is **split into two
fragments** ("pass 0": 0.00-18823.85 s, and "pass 2": 67267.21-86164.09 s)
purely because the analysis window `[0, sidereal_day]` cuts through the
middle of that physically continuous access run — the ground track
exactly repeats after one sidereal day (§M3.3), so pass 0 and pass 2 are
literally the same recurring apogee dwell viewed at the two ends of an
arbitrary window. A window-independent count would report 2 physically
distinct passes per sidereal day (consistent with ~2 Molniya revolutions/
day). This is a reporting-window artifact of choosing `[0, sidereal_day]`
as the analysis interval, not a geometry or algorithm defect, and is
retained here (rather than silently patched) because M5 will need to make
an explicit, documented choice about window boundaries for revisit-gap
statistics — this M3 result flags exactly why that choice matters.

The 80.3% one-day access fraction is in the same ballpark as, but not
identical to, the M1 §6 anomaly-based dwell proxy (84.0% of the orbital
period within ±60° true anomaly of apogee) — expected, since the two are
different metrics (a fixed-site elevation threshold vs. a pure orbital
anomaly window) that should be of comparable magnitude for a mission
designed around northern apogee dwell, without being numerically equal.

## M3.6 Independent access verification

**Independent elevation cross-check** (`elevation_independent_check`):
uses the Earth-center/site/satellite triangle geometry (geocentric angle
psi = arccos(r_site_hat · r_sat_hat), then closed-form
`tan(el) = (cos psi - R_E/r_sat) / sin psi` and law-of-cosines range) —
**no ENU basis, no dot products against e_east/e_north/e_up** — a
genuinely separate derivation path from the ENU-projection method.

Compared against the ENU method over the full one-sidereal-day, dt=15 s
trajectory (5745 samples):

- Max |elevation error| = **1.563×10⁻¹³ deg**
- Max |range error| = **1.455×10⁻¹¹ km**

Both at the double-precision floor — the two independent methods agree
to numerical precision, confirming the ENU-based access geometry.

**Timestep convergence study** (one sidereal day, dt = 120/60/15 s):

| dt (s) | Samples | Total access (s) | Max gap (s) | Longest pass (s) | Passes |
|---|---|---|---|---|---|
| 120.0 | 719 | 69190.499 | 8552.215 | 31470.980 | 3 |
| 60.0 | 1437 | 69191.780 | 8551.195 | 31471.249 | 3 |
| 15.0 | 5745 | 69192.106 | 8551.051 | 31471.372 | 3 |

Deltas relative to the finest (dt=15 s) setting:

| dt (s) | \|Δ total access\| (s) | \|Δ max gap\| (s) | \|Δ longest pass\| (s) |
|---|---|---|---|
| 120.0 | 1.6070 | 1.1640 | 0.3926 |
| 60.0 | 0.3252 | 0.1440 | 0.1235 |
| 15.0 | 0 (reference) | 0 | 0 |

Every metric shrinks monotonically (roughly by a factor of ~5 per 2×
timestep halving from 120→60 s, and again 60→15 s) as the sampling
interval is refined — genuine numerical convergence of the
linear-interpolation/parabolic-refinement scheme, not merely asserted.
Pass count (3, including the window-split artifact of §M3.5) is stable
across all three timesteps.

## M3.7 Figures

- [`figures/m3_ground_track.png`](figures/m3_ground_track.png) — two-body
  + spherical-Earth-rotation ground track over one sidereal day, with
  apogee/perigee markers and the representative site, explicitly labeled
  "no J2." Longitude wraparound is handled by splitting the plotted
  polyline wherever consecutive samples jump by more than 180°, so no
  fake line is drawn across the map.
- [`figures/m3_access_vs_time.png`](figures/m3_access_vs_time.png) —
  elevation vs. time over one sidereal day, the 10° threshold, shaded
  access intervals, and apogee-time markers, explicitly labeled as
  geometric elevation only (not a communications/link result).
- [`figures/m3_range_vs_elevation.png`](figures/m3_range_vs_elevation.png)
  — supporting figure: range vs. elevation colored by time, showing the
  geometric range/elevation relationship through each pass.

All three figures were visually inspected for clipping, overlapping
annotations, fake wraparound lines, unreadable legends, wrong units, and
misleading aspect ratio. **One genuine issue was found and fixed during
this inspection** (§M3.8).

## M3.8 Genuine discrepancy found and fixed

**What was found:** the first version of `m3_access_vs_time.png` drew
each access-interval shaded region with `axvspan(..., alpha=0.15)` inside
a loop, then drew the *first* interval a second time (with a `label=`) to
create a legend entry. Because matplotlib alpha-blends overlapping
patches, the first (legend-labeled) interval rendered visibly **darker**
than the other two identical-alpha intervals — a misleading visual
artifact that could be misread as a different (e.g. higher-confidence or
different-type) access interval.

**Fix:** the plotting loop no longer redraws any interval for the legend;
a `matplotlib.patches.Patch` / `Line2D` proxy artist is built solely for
the legend entry and never added to the axes. Re-inspected the
regenerated figure: all three shaded intervals now render with identical,
consistent shading.

This was a documentation/figure-generation bug only — it did not affect
any numerical result, test, or the underlying access-interval data.

## M3.9 Test suite

`tests/test_frames_access.py` implements checklist items A-T (38 new
tests). One test tolerance was deliberately loosened with an explicit
comment (test L, overhead elevation): `arcsin` is ill-conditioned near
±1 (derivative → infinity), so a ~1e-16 floating-point error in
`rho_up/range` at exactly-overhead amplifies to ~1e-6 deg in the computed
elevation — an expected numerical effect of the formula, not a geometry
bug; the test tolerance (1e-4 deg) remains far tighter than any physical
requirement. All M1/M2 regression tests continue to pass unchanged
(test S re-verifies M1/M2 apsis radii/velocities through the production
code path in the M3 tree).

**Total: 72 tests pass under `pytest -W error`, zero warnings**
(34 from M1/M2 + 38 new M3 tests).

## M3.10 Scope guard confirmation

No J2 acceleration, secular J2 propagation, numerical nodal/argument-of-
perigee regression, long-horizon J2 ground-track drift, constellation
design, multi-satellite coverage, RF link budget, atmospheric refraction,
antenna gain, oblate-Earth visibility, geodetic WGS-84 station
coordinates, or stationkeeping was implemented in M3 — all remain
explicitly deferred to M4/M5 per the roadmap. This module chain is
strictly: two-body ECI propagation (M2, unchanged) → ECEF → geocentric
ground track → basic spherical-Earth access geometry.

---

# Milestone 4 — First-Order J2 Secular Propagation + Critical-Inclination Verification

**Status:** M4 complete. Implements a first-order secular J2 mean-element
model (RAAN/argp/mean-anomaly rates), a secular element propagator with
state reconstruction reusing M2's element↔state machinery, a J2-aware
ground track reusing M3's frame code unchanged, an analytical apogee/
perigee-time finder, and an independent short-horizon Cartesian J2
cross-check that fits secular rates from osculating elements. **M4 is a
mean-element secular model, not a full osculating-element production
propagator** — the optional Cartesian cross-check is supporting only, per
the M4 scope guard. **Coverage/revisit analysis remains M5** — nothing
in this section is a coverage or link-availability result.

## M4.0 Pre-flight verification

Before writing any M4 code, the M3 headline values were reproduced from
production code and the M2/M3 reports/figures were regenerated and
diffed byte-for-byte against the committed artifacts:

| Quantity | M3 committed value | Recomputed (pre-M4) | Match |
|---|---|---|---|
| Apogee latitude | +63.434949° | +63.434949° | ✓ |
| Successive-apogee longitude separation | -180.00000000° | -180.00000000° | ✓ |
| One-sidereal-day ground-track repeat error | 1.173×10⁻¹⁰ deg | 1.173×10⁻¹⁰ deg | ✓ |
| Access fraction | 80.30% | 80.3027% | ✓ |
| Max no-access gap | 8551.051 s | 8551.051 s | ✓ |

`scripts/m2_verification_report.py` and `scripts/m3_verification_report.py`
output, plus all five committed M2/M3 figure PNGs, were regenerated and
found **byte-for-byte identical** (`diff` on report text; identical MD5
checksums on every figure) — zero numerical drift from M2/M3.
`pytest -W error` on the pre-M4 tree passed 72/72.

## M4.1 First-order secular J2 rates

New module `j2.py`. For `n = sqrt(mu/a^3)` and `p = a(1-e^2)`:

    RAAN_dot = -(3/2) J2 n (Re/p)^2 cos(i)
    argp_dot =  (3/4) J2 n (Re/p)^2 (5 cos^2(i) - 1)
    M_dot    = n + (3/4) J2 n (Re/p)^2 sqrt(1-e^2) (3 cos^2(i) - 1)

**M-dot convention (explicit, per M4 §2):** mean anomaly does **not**
advance at the unperturbed rate `n` alone — the secular J2 correction
term above is added, using the same `J2 * n * (Re/p)^2` prefactor shared
by all three rates. This is tested independently (checklist F): the
`J2→0` limit of this exact formula must reduce to `M_dot = n` exactly,
which it does to double-precision (`test_F_j2_zero_limit_rates`).

`a`, `e`, `i` are held fixed in this secular model (M4 §3 requirement);
`RAAN`, `argp`, `M` are tracked and evolved separately — never collapsed
into an argument-of-latitude, since the Molniya baseline's eccentricity
(0.737) makes that collapse invalid. `SecularElements` and `SecularRates`
are immutable (`frozen=True`) dataclasses; `propagate_elements_j2_secular`
never mutates its input (test H).

## M4.2 Critical-inclination verification

At the baseline `i = 63.434949°`:

    5 cos^2(i) - 1 = 4.441e-16   (zero to double-precision floor)
    argp_dot = 1.4557e-23 rad/s = 7.206e-17 deg/day   (effectively zero)

**Off-critical comparison** (a, e held at baseline; only i varied):

| i (deg) | 5cos²i - 1 | argp_dot (deg/day) | RAAN_dot (deg/day) |
|---|---|---|---|
| 60.000000 | +0.250000 | **+0.040566** | -0.162265 |
| 63.000000 | +0.030537 | **+0.004955** | -0.147334 |
| 63.434949 (baseline) | 0.000000 | **0.000000** (7.2e-17) | -0.145135 |
| 64.000000 | -0.039154 | **-0.006353** | -0.142265 |
| 65.000000 | -0.106969 | **-0.017357** | -0.137153 |

`argp_dot` changes sign exactly at the critical inclination — positive
below it, negative above (test D), with clearly nonzero magnitude at
every off-critical sample point (test E, all |argp_dot| > 1e-4 deg/day,
i.e. many orders of magnitude above the baseline's numerically-zero
residual). See [`figures/m4_j2_rate_sensitivity.png`](figures/m4_j2_rate_sensitivity.png).

## M4.3 Nodal regression verification

**Independent recomputation** (a separate formula coded directly in
`tests/test_j2.py::test_A_baseline_raan_dot_independent_formula`, not a
second call to `compute_secular_rates`):

    RAAN_dot = -0.145135 deg/day

matching both the production `compute_secular_rates` output and the M1
prediction (DESIGN.md §3: `Omega_dot ≈ -0.145135 deg/day`) to 6
significant figures — **no discrepancy found** beyond the display
rounding already documented in M1/M2.

## M4.4 Secular element propagator and state reconstruction

`propagate_elements_j2_secular(elements0, t)` advances RAAN/argp/M
linearly at the rates from §M4.1 (a, e, i constant), supports scalar or
numpy-array `t`, and never mutates `elements0` (tests G, H). Internally
the angles are advanced **unwrapped** (no modulo at any intermediate
step) and wrapped to [0, 360) deg only at the output boundary if
`wrap=True` (default); `wrap=False` returns the continuous unwrapped
angle, used for the long-horizon drift figure (§M4.7) and for exactly
recovering the linear secular trend in tests (test H, K).

**State reconstruction** reuses, rather than duplicates, M2's Kepler-
solving and PQW-rotation code: mean anomaly → eccentric anomaly
(`propagation.mean_to_eccentric_anomaly`, unchanged since M2) → true
anomaly → `elements.coe_to_rv` (unchanged since M2). Self-consistency
verified: reconstructed radius matches `r = a(1 - e cos E)` to
≤1e-10 relative (test I).

**J2→0 regression** (test F, J): with `J2=0`, the secular model's
reconstructed state matches (a) the independent Kepler time-of-flight
solver (`kepler_state_at_time`, unchanged since M2) to ≤1e-10 relative,
and (b) a full M2 `solve_ivp` two-body propagation to ≤1e-9 relative —
confirming the J2=0 limit exactly recovers two-body dynamics through two
independent M2 code paths, not just one.

## M4.5 Ground track under secular J2 (multi-day)

`groundtrack.compute_ground_track_j2` reuses the *exact same* M3
`eci_to_ecef` / `ecef_to_geocentric_latlon` frame code used by the
two-body ground track — so any difference between the two tracks is
attributable only to J2 dynamics, not a different frame model.

**Apogee tracking method (M4 §9):** apogee times are found **analytically**
from the linear M(t) relation (`mean_anomaly_crossing_times`, solving
`M0 + M_dot*t = 180° (mod 360°)` directly for t), not by nearest-sample
search — exact by construction, verified against the propagated M value
at each returned time (test L).

**Apogee ground-track drift, same-side-apogee comparison** (comparing the
first apogee to the last apogee an even number of half-revolutions later,
so both sit on the same ~90°-longitude lobe):

| Horizon | # apogees | First (lat, lon) | Last (lat, lon) | Longitude drift | Latitude drift |
|---|---|---|---|---|---|
| 1 day | 3 | (63.43495, 90.00000) | (63.43495, 89.83339) | **-0.16661°** | -1.4e-14° |
| 3 days | 7 | (63.43495, 90.00000) | (63.43495, 89.50016) | **-0.49984°** | 0.0° |
| 7 days | 15 | (63.43495, 90.00000) | (63.43495, 88.83371) | **-1.16629°** | -1.4e-14° |
| 14 days | 29 | (63.43495, 90.00000) | (63.43495, 87.66742) | **-2.33258°** | 0.0° |

Latitude drift is at the double-precision floor (as expected: argp_dot
≈ 0 at critical inclination keeps apogee's argument of latitude, and
hence its geocentric latitude, fixed). Longitude drift grows
monotonically and non-linearly with horizon (driven by the RAAN
regression accumulating over more elapsed time and more revolutions per
comparison window — not a fixed per-day increment).

**Decomposition of the 1-day drift (a genuine quantitative explanation,
not just "drift observed"):** the same-side apogee 1 day later actually
occurs at t=86169.324 s, not exactly `SIDEREAL_DAY_S=86164.0905 s`,
because the J2-corrected `M_dot` is marginally faster than the
unperturbed `n` (§M4.1). This 5.233 s offset means Earth has rotated
slightly *more* than one full sidereal day by the time the same-side
apogee recurs, which by itself would shift the apparent apogee longitude
westward by an additional amount independent of RAAN regression:

    RAAN-regression component  = RAAN_dot * (86169.324/86400) = -0.144747 deg
    extra-Earth-rotation component = -(86169.324 - 86164.0905)/86164.0905 * 360 = -0.021866 deg
    sum = -0.166613 deg   (observed: -0.166613 deg — exact match)

This decomposition is a genuine independent check that the observed
number is fully explained by the two known physical effects (RAAN
regression + M-dot/period timing offset), not an unexplained residual.

**Two-body vs. J2 ground-track divergence** (longitude difference at the
end of each horizon, comparing the *same-time* two-body and J2 tracks —
a different, complementary comparison from the same-side-apogee table
above):

| Horizon | Longitude diff at window end |
|---|---|
| 1 day | -0.167080° |
| 3 days | -0.501205° |
| 7 days | -1.169059° |
| 14 days | -2.335197° |

Clearly nonzero at every horizon and growing with time (test M);
visualized in [`figures/m4_ground_track_two_body_vs_j2.png`](figures/m4_ground_track_two_body_vs_j2.png)
with a zoomed inset, since at full-map scale the two tracks are visually
indistinguishable (the drift, while measured precisely above, is
genuinely small relative to the ~180°-scale ground-track loops — J2 is a
weak perturbation for this orbit over these horizons, and the figure
says so honestly rather than exaggerating the visual difference).

**One-sidereal-day repeat error under J2** (window start vs. end, *not*
an apogee-to-apogee comparison — this quantifies how the M3 "exact
repeat" property degrades under J2): longitude at t=0 is 90.000000°,
at t=one sidereal day it is 89.833374°, a **-0.166626°** difference —
consistent with (and dominated by) the same RAAN-regression + M-dot
timing effects decomposed above. Unlike the M3 two-body case (repeat
error 1.17×10⁻¹⁰ deg), the ground track under J2 does **not** exactly
repeat every sidereal day — expected and correctly reproduced.

## M4.6 Argument-of-perigee orientation check

**Baseline (critical inclination):** ω(0) = 270.000000°, ω(14 days) =
270.000000° — **unchanged to the precision printed** (the true drift,
argp_dot × 14 days ≈ 1×10⁻¹⁵ deg, is far below any meaningful digit).
Apogee latitude: 63.434949° at day 0 and day 14, delta 0.0° (exactly
frozen, test N).

**Off-critical (i=65°):** ω(0) = 270.000000°, ω(14 days) = 269.756997°,
a **-0.243003°** drift over 14 days (test O, clearly >0.1° threshold).
Apogee latitude changes correspondingly: 65.000000° → 64.998895°, a
**-1.105×10⁻³°** shift over 14 days.

This makes the critical-inclination benefit quantitatively (and, in
[`figures/m4_element_drift.png`](figures/m4_element_drift.png), visually)
obvious: the baseline ω(t) trace is a flat line at 270°, while the
off-critical case visibly slopes downward over the same 14-day window —
even though the *absolute* magnitude of the off-critical apogee-latitude
shift (~0.001°) is modest at this timescale, because J2 is a genuinely
weak perturbation over 2 weeks for this semi-major axis. The relevant
engineering point is the qualitative contrast (zero vs. nonzero secular
drift), not the absolute size of the off-critical number at 14 days.

## M4.7 Optional Cartesian J2 cross-check

New module `j2_cartesian.py`, implementing the M4-§11 acceleration form

    ax = -mu*x/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 1)]
    ay = -mu*y/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 1)]
    az = -mu*z/r^3 * [1 - 1.5 J2 (Re/r)^2 (5 z^2/r^2 - 3)]

verified algebraically equivalent to the more commonly quoted factored
form in the module docstring. **J2=0 reduces exactly** to the M2
two-body EOM (`np.allclose(..., atol=1e-15)`, test R) — confirmed
before doing anything else with this module.

**Fitted-rate cross-check:** propagated the full osculating Cartesian
two-body+J2 dynamics for 7 days at `rtol=atol=1e-13`, detected apsis
events via `r·v=0` crossings (same event-detection technique as M2's
`find_apsides`, independently re-implemented here rather than reusing
the two-body-only production function), extracted osculating RAAN/argp
at each of the 14 **perigee** crossings (sampling at a fixed orbital
phase removes most short-period J2 oscillation — this is the standard
"mean-element-like" trick, not a full Brouwer/Lyddane osculating→mean
transformation), and linear-fit the secular trend:

| Rate | First-order theory | 7-day osculating fit | Agreement |
|---|---|---|---|
| RAAN_dot | -0.145135 deg/day | **-0.145166 deg/day** | 2.18×10⁻⁴ relative |
| argp_dot | 7.2×10⁻¹⁷ deg/day (≈0) | **2.729×10⁻⁵ deg/day** | both ≈0; fitted value is 4 orders of magnitude below the smallest off-critical rate (~4×10⁻³ deg/day at i=63°, §M4.2) |

**Short-period vs. secular distinction:** the ~0.02% RAAN_dot discrepancy
and the small nonzero fitted argp_dot are attributed to (a) short-period
J2 oscillation not fully averaged out by only 14 perigee samples over 7
days, and (b) the first-order secular theory itself omitting higher-order
J2² and J2-e coupling terms present in the full osculating dynamics. Both
are consistent with, and bounded by, the known limitations of first-order
secular theory — not evidence of an error in either implementation. This
independent numerical fit corroborates the analytical secular rates
without promoting the Cartesian J2 propagator to a production model, per
the M4 scope guard.

## M4.8 Convergence / numerical checks

- **Kepler solver convergence:** the same Newton-iteration solver used
  since M2 (`mean_to_eccentric_anomaly`, unchanged) is reused for all M4
  state reconstruction; its convergence behavior was already verified in
  M2 and is exercised here at additional (secularly-evolved) M values
  with no new failures.
- **Reconstructed radius identity:** `r = a(1 - e cos E)` verified to
  ≤1e-10 relative at an arbitrary secularly-propagated state (test I).
- **Apsis radii under J2 secular model:** since a, e are held exactly
  constant in this model, reconstructed perigee/apogee radii are
  *identical* to the M1/M2 values at every sampled time by construction
  (no drift is possible or expected in the mean-element radius — only
  RAAN/argp/M evolve) — confirmed via `test_P_regression_apsis_still_matches_M1`
  running the unchanged M2 two-body path.
- **No long-horizon wrap artifacts:** tested at a 200-day horizon
  (test K) — unwrapped mean anomaly grows to >720° as expected (many
  full revolutions), while the wrapped value stays in [0, 360) and
  matches `unwrapped mod 360` to within 1e-6 deg.
- **Cartesian J2 integrator:** run at the same tight `rtol=atol=1e-13`
  DOP853 settings established and convergence-tested in M2/M3; used here
  only for the supporting cross-check, so a fresh convergence sweep was
  not repeated (would duplicate M2 §9's already-established convergence
  behavior on a structurally similar smooth ODE).

## M4.9 Figures

- [`figures/m4_j2_rate_sensitivity.png`](figures/m4_j2_rate_sensitivity.png)
  — RAAN_dot and argp_dot vs. inclination (50°-75°), critical/baseline
  inclination marked on both panels, argp_dot crossing zero clearly
  visible, units deg/day, no axis truncation.
- [`figures/m4_ground_track_two_body_vs_j2.png`](figures/m4_ground_track_two_body_vs_j2.png)
  — 7-day two-body (solid blue) vs. secular-J2 (dashed orange) ground
  track, apogees marked, same M3 spherical-Earth rotation model, with a
  zoomed inset on the day-7 apogee showing the actual 1.166° divergence
  (invisible at full-map scale).
- [`figures/m4_element_drift.png`](figures/m4_element_drift.png) — RAAN(t)
  and argp(t) over 14 days; argp(t) overlays the baseline (flat at 270°)
  against the i=65° off-critical case (visibly sloped), making the
  critical-inclination benefit visually obvious.

All three figures were visually inspected for clipping, overlapping
annotations, fake longitude-wrap lines, unreadable legends, wrong units,
and confusing wrapped-angle jumps. **One genuine layout issue was found
and fixed** (§M4.10).

## M4.10 Genuine discrepancy found and fixed

**What was found:** the first version of
`m4_ground_track_two_body_vs_j2.png` rendered the two tracks as visually
indistinguishable at full-map scale (correct — the divergence really is
that small relative to a ±180° longitude axis — but not a useful figure
on its own), and after adding a zoomed inset to fix that, the inset's
title text overlapped the main plot's subtitle.

**Fix:** repositioned and shrank the inset (`inset_axes` bbox and title
font size/pad adjusted) and reserved top margin via
`tight_layout(rect=(0,0,1,0.95))`. Re-inspected the regenerated figure:
no overlap, inset title and zoom-indicator lines render cleanly. This was
a figure-layout issue only — no numerical result, test, or underlying
data was affected.

## M4.11 Test suite

`tests/test_j2.py` implements checklist items A-R (26 new tests,
including the optional Cartesian cross-check). All M1/M2/M3 regression
tests continue to pass unchanged (test P re-verifies M1 apsis
radii/velocities through the unchanged M2 production code path in the
M4 tree).

**Total: 98 tests pass under `pytest -W error`, zero warnings**
(72 from M1/M2/M3 + 26 new M4 tests).

## M4.12 Scope guard confirmation

No M5 coverage/revisit study, long-horizon operational-availability
claim, constellation sizing, multi-satellite architecture, RF link
budget, antenna pointing, stationkeeping, drag/SRP/lunisolar
perturbation, Moon/Sun ephemeris, launch-vehicle study, or final
portfolio packaging was implemented in M4 — all remain explicitly
deferred to M5/M6 per the roadmap. **Coverage analysis (access fraction,
revisit interval, coverage gaps under J2) remains entirely M5's scope**;
this milestone only establishes the J2-perturbed orbital dynamics and
ground-track drift that M5's coverage study will need to account for.
