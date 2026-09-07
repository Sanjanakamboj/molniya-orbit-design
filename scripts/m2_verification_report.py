"""M2 verification report generator.

Runs the propagation/element-conversion checks and prints a numeric
summary used to populate the M2 section of DESIGN.md. Also produces the
two M2 diagnostic figures in figures/.

Run:
    python scripts/m2_verification_report.py
"""

from __future__ import annotations

import numpy as np

from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    M1_APOGEE_ALT_KM,
    M1_PERIGEE_ALT_KM,
    M1_PERIOD_S,
    M1_RA_KM,
    M1_RP_KM,
    M1_VA_KM_S,
    M1_VP_KM_S,
    R_EARTH,
)
from molniya_design.elements import coe_to_rv, rv_to_coe
from molniya_design.propagation import (
    find_apsides,
    geocentric_latitude_deg,
    kepler_state_at_time,
    propagate,
)
from molniya_design.twobody import (
    eccentricity_vector,
    specific_angular_momentum,
    specific_energy,
)

BASE = BASELINE_ELEMENTS_DEG


def main():
    print("=" * 78)
    print("M2 VERIFICATION REPORT")
    print("=" * 78)

    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )
    print(f"\nBaseline state at epoch (apogee, nu0=180 deg):")
    print(f"  r0 = {r0} km, |r0| = {np.linalg.norm(r0):.6f} km")
    print(f"  v0 = {v0} km/s, |v0| = {np.linalg.norm(v0):.6f} km/s")
    print(f"  M1 ra = {M1_RA_KM:.6f} km, M1 va = {M1_VA_KM_S:.6f} km/s")

    # --- round trip ---
    coe = rv_to_coe(r0, v0)
    print(f"\nRound-trip element recovery (element->state->element):")
    for k in ["a_km", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg"]:
        print(f"  {k}: input={BASE.get(k, BASE.get('nu0_deg') if k=='nu_deg' else None)}  recovered={coe[k]}")

    # --- apsides over ~2 periods ---
    t_span = (0.0, 2.05 * M1_PERIOD_S)
    apsides = find_apsides(r0, v0, t_span, rtol=1e-13, atol=1e-13)
    print(f"\nDetected {len(apsides)} apsides over 2.05 periods:")
    for a in apsides:
        print(f"  t={a.t_s:12.4f} s  {a.kind:8s} r={a.r_km:.6f} km  v={a.v_km_s:.6f} km/s")

    apogee_times = [a.t_s for a in apsides if a.kind == "apogee" and a.t_s > 1.0]
    numerical_period = apogee_times[0]
    print(f"\nNumerical period (apogee-to-apogee) = {numerical_period:.6f} s "
          f"({numerical_period/3600:.6f} h)")
    print(f"M1 analytical period                = {M1_PERIOD_S:.6f} s ({M1_PERIOD_S/3600:.6f} h)")
    print(f"Relative error                       = {abs(numerical_period-M1_PERIOD_S)/M1_PERIOD_S:.3e}")

    perigee_r = [a.r_km for a in apsides if a.kind == "perigee"]
    apogee_r = [a.r_km for a in apsides if a.kind == "apogee"]
    print(f"\nDetected perigee radii: {perigee_r}")
    print(f"  M1 rp = {M1_RP_KM:.6f} km; max abs error = {max(abs(r-M1_RP_KM) for r in perigee_r):.3e} km")
    print(f"Detected apogee radii: {apogee_r}")
    print(f"  M1 ra = {M1_RA_KM:.6f} km; max abs error = {max(abs(r-M1_RA_KM) for r in apogee_r):.3e} km")

    perigee_v = [a.v_km_s for a in apsides if a.kind == "perigee"]
    apogee_v = [a.v_km_s for a in apsides if a.kind == "apogee"]
    print(f"\nDetected perigee speeds: {perigee_v}")
    print(f"  M1 vp = {M1_VP_KM_S:.6f} km/s; max abs error = {max(abs(v-M1_VP_KM_S) for v in perigee_v):.3e} km/s")
    print(f"Detected apogee speeds: {apogee_v}")
    print(f"  M1 va = {M1_VA_KM_S:.6f} km/s; max abs error = {max(abs(v-M1_VA_KM_S) for v in apogee_v):.3e} km/s")

    # --- one-period closure ---
    sol = propagate(r0, v0, (0.0, M1_PERIOD_S), rtol=1e-13, atol=1e-13)
    r_final, v_final = sol.y[0:3, -1], sol.y[3:6, -1]
    r_err = np.linalg.norm(r_final - r0) / np.linalg.norm(r0)
    v_err = np.linalg.norm(v_final - v0) / np.linalg.norm(v0)
    print(f"\nOne-period closure (tight tol rtol=atol=1e-13):")
    print(f"  relative position error = {r_err:.3e}")
    print(f"  relative velocity error = {v_err:.3e}")

    # --- conservation over 2 periods ---
    t_eval = np.linspace(0.0, 2.0 * M1_PERIOD_S, 2000)
    sol2 = propagate(r0, v0, (0.0, 2.0 * M1_PERIOD_S), t_eval=t_eval, rtol=1e-13, atol=1e-13)
    energies = np.array([specific_energy(sol2.y[0:3, k], sol2.y[3:6, k]) for k in range(sol2.y.shape[1])])
    h_mags = np.array([
        np.linalg.norm(specific_angular_momentum(sol2.y[0:3, k], sol2.y[3:6, k]))
        for k in range(sol2.y.shape[1])
    ])
    e_vecs = np.array([eccentricity_vector(sol2.y[0:3, k], sol2.y[3:6, k]) for k in range(sol2.y.shape[1])])
    e_drift = np.max(np.linalg.norm(e_vecs - e_vecs[0], axis=1))

    print(f"\nConservation over 2 periods ({len(t_eval)} samples, rtol=atol=1e-13):")
    print(f"  specific energy: min={energies.min():.9f}, max={energies.max():.9f}, "
          f"rel drift={(energies.max()-energies.min())/abs(energies[0]):.3e}")
    print(f"  |h|: min={h_mags.min():.6f}, max={h_mags.max():.6f}, "
          f"rel drift={(h_mags.max()-h_mags.min())/h_mags[0]:.3e}")
    print(f"  eccentricity vector max drift = {e_drift:.3e}")

    # --- independent Kepler cross-check ---
    print(f"\nIndependent Kepler time-of-flight cross-check:")
    max_r_err, max_v_err = 0.0, 0.0
    for frac in [0.0, 0.1, 0.25, 0.5, 0.63, 0.9]:
        t = frac * M1_PERIOD_S
        if t == 0.0:
            r_num, v_num = r0, v0
        else:
            s = propagate(r0, v0, (0.0, t), rtol=1e-13, atol=1e-13)
            r_num, v_num = s.y[0:3, -1], s.y[3:6, -1]
        r_kep, v_kep = kepler_state_at_time(
            BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
            BASE["argp_deg"], BASE["nu0_deg"], t,
        )
        r_err = np.linalg.norm(r_num - r_kep) / np.linalg.norm(r_kep)
        v_err = np.linalg.norm(v_num - v_kep) / np.linalg.norm(v_kep)
        max_r_err = max(max_r_err, r_err)
        max_v_err = max(max_v_err, v_err)
        print(f"  frac={frac:.2f}  t={t:10.2f} s  r_rel_err={r_err:.3e}  v_rel_err={v_err:.3e}")
    print(f"  max relative position error = {max_r_err:.3e}")
    print(f"  max relative velocity error = {max_v_err:.3e}")

    # --- northern/southern apogee orientation ---
    r270, _ = coe_to_rv(BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], 270.0, 180.0)
    r90, _ = coe_to_rv(BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], 90.0, 180.0)
    lat270 = geocentric_latitude_deg(r270)
    lat90 = geocentric_latitude_deg(r90)
    print(f"\nNorthern-apogee orientation check:")
    print(f"  omega=270 deg -> apogee latitude = {lat270:.6f} deg (expect +{BASE['i_deg']:.4f})")
    print(f"  omega=90  deg -> apogee latitude = {lat90:.6f} deg (expect -{BASE['i_deg']:.4f})")

    # --- convergence study ---
    print(f"\nTolerance convergence study (one-period closure + energy drift):")
    settings = [("loose", 1e-6, 1e-6), ("medium", 1e-9, 1e-9), ("tight", 1e-12, 1e-12)]
    for name, rtol, atol in settings:
        s = propagate(r0, v0, (0.0, M1_PERIOD_S), rtol=rtol, atol=atol)
        r_err = np.linalg.norm(s.y[0:3, -1] - r0) / np.linalg.norm(r0)
        e0 = specific_energy(s.y[0:3, 0], s.y[3:6, 0])
        ef = specific_energy(s.y[0:3, -1], s.y[3:6, -1])
        e_drift = abs(ef - e0) / abs(e0)
        h0 = np.linalg.norm(specific_angular_momentum(s.y[0:3, 0], s.y[3:6, 0]))
        hf = np.linalg.norm(specific_angular_momentum(s.y[0:3, -1], s.y[3:6, -1]))
        h_drift = abs(hf - h0) / h0

        # apsis errors at these tolerances
        aps = find_apsides(r0, v0, t_span, rtol=rtol, atol=atol)
        ap_r = [a.r_km for a in aps if a.kind == "apogee"]
        pe_r = [a.r_km for a in aps if a.kind == "perigee"]
        ap_err = max(abs(r - M1_RA_KM) for r in ap_r)
        pe_err = max(abs(r - M1_RP_KM) for r in pe_r)

        print(f"  {name:6s} (rtol=atol={rtol:.0e}): closure_r_err={r_err:.3e}  "
              f"energy_drift={e_drift:.3e}  h_drift={h_drift:.3e}  "
              f"apogee_r_err={ap_err:.3e} km  perigee_r_err={pe_err:.3e} km")

    print("\n" + "=" * 78)
    print("END M2 VERIFICATION REPORT")
    print("=" * 78)


if __name__ == "__main__":
    main()
