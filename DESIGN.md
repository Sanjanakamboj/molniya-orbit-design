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

- **M1** — Analytical design + verification plan (this document). ✅ current
- **M2** — Element/state conversion + two-body propagation + apsis/period
  numerical verification (§10.A–F).
- **M3** — ECI/ECEF transformation + Earth-fixed ground track + access
  geometry (§10.J, K, L, N).
- **M4** — First-order J2 secular propagation + critical-inclination
  verification (§10.G, H, I).
- **M5** — High-latitude dwell/coverage/revisit + parameter sensitivity
  study.
- **M6** — Independent validation + portfolio polish + CI/reproducibility
  audit.

M2 is **not** started as part of this milestone.
