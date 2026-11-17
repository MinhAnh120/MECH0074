"""
run_stage4.py - ONE-COMMAND RUNNER for MECH0074 CW1 Stage 4 (coupled response)
==============================================================================
Usage (from the folder containing stage4.py):
    python run_stage4.py                              # uses the built-in MOC solver
    python run_stage4.py --stage1 path/to/MECH0074    # uses Minh-Anh's Stage-1 solver (recommended)
    python run_stage4.py --out results_stage4

EDIT ONLY THE TWO BLOCKS BELOW (GEOMETRY and CASES). Everything else is automatic.
Outputs (in --out): comparison table CSV, per-case time histories, figures, verification log.
"""
import argparse, csv, json, os
from dataclasses import replace
import numpy as np
import stage4 as S

# ----------------------------------------------------------------------------
# BLOCK 1 - GEOMETRY and agreed inputs.  Replace the placeholders with Stage 3's real loop.
# ----------------------------------------------------------------------------
BASE = S.Par(
    c=1100.0,                  # GROUP DECISION: bulk sound speed
    f_darcy=0.01010783190,     # Stage 1 nominal Colebrook friction factor
    # --- loop layout: PLACEHOLDER until Stage 3 supplies the real geometry ---
    H=10.0, W=10.0,            # loop leg lengths (m)
    La=30.0, Lb=30.0,          # straight runs before / after the loop (m)
    guides=(),                 # lateral restraints, e.g. ((45.0, 1),) = y-guide at s=45 m
    # --- modelling choices to justify in the talk ---
    zeta=0.01,                 # damping ratio (sensitivity study runs 0.5-5 %)
    flex_bends=False,          # True = Karman bend flexibility
    sif=False,                 # True = stress intensification on bend stress
    T_end=4.0, dt=5e-4,
)

# ----------------------------------------------------------------------------
# BLOCK 2 - The designs to compare. First one must be the baseline.
# ----------------------------------------------------------------------------
CASES = [
    ("baseline",            BASE),
    ("slow closure 0.3 s",  replace(BASE, t_close=0.3)),
    ("larger radius 2.5 m", replace(BASE, R=2.5)),
    ("larger loop 14x14 m", replace(BASE, H=14.0, W=14.0)),
]

PHASE_LIMIT_PA = 95873.89895   # Stage 1 methane saturation screen at 111 K (absolute Pa)


def hydraulics(p, stage1_repo):
    """Pressure/velocity histories at the bends for this case (Stage-1 solver if available)."""
    geo = S.build_loop(p)
    if stage1_repo:
        import sys, importlib
        # the Stage-1 solver takes closure time from its own shared_case, so make that match this case
        sys.path.insert(0, stage1_repo)
        from dataclasses import replace as rep
        import complete_stage1, fluid_assumptions, stage1
        base = complete_stage1.shared_case("bulk", 0)
        f = fluid_assumptions.friction_estimate(base.density, base.velocity, base.diameter)["darcy_factor"]
        case = rep(base, friction=f, length=geo["L_hyd"], closure_time=p.t_close)
        sts = sorted({round(float(x), 9) for b in geo["bends"] for x in (b["s_start"], b["s_end"])})
        N = int(round(geo["L_hyd"] / 0.3125))
        r = stage1.simulate(case, N, p.T_end, stations=sts)
        return S.HydAtStations(r["time"], np.array(sts), r["pressure"], r["flow"] / case.area), r
    return S.moc(p, geo["L_hyd"]), None


def first_crossing(hyd, limit):
    """First time any exported station drops below the phase limit (None if never)."""
    below = (hyd.P <= limit).any(axis=1)
    return float(hyd.t[np.argmax(below)]) if below.any() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage1", default=None, help="path to the cloned Stage-1 repo")
    ap.add_argument("--out", default="results_stage4")
    ap.add_argument("--skip-sensitivity", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = {}

    # ---- 1. hand calculations to quote on the slides ----
    log["hand_calcs"] = {k: float(v) for k, v in S.joukowsky_numbers(BASE).items()}
    kf, sif = S.bend_flex_factor(BASE)
    log["bend_flexibility_factor"], log["SIF"] = float(kf), float(sif)
    log["restrained_thermal_stress_MPa"] = BASE.E * BASE.alpha * BASE.dTc / 1e6
    print("Hand calcs:", json.dumps(log["hand_calcs"], indent=1))

    # ---- 2. verification checks (put these in the 'verification' slide) ----
    log["cantilever_check"] = {}
    for L in (4.0, 10.0, 20.0):
        fe, th = S.cantilever_check(BASE, L)
        log["cantilever_check"][str(L)] = dict(FE=fe.tolist(), EulerBernoulli=th.tolist(),
                                               max_err_pct=float(np.max(abs(fe - th) / th) * 100))
    log["uniform_bend_check"] = {k: float(v) for k, v in
                                 S.check_uniform_bend(BASE, S.build_loop(BASE)["bends"][0]).items()}
    print("Cantilever max error %:", {k: round(v["max_err_pct"], 6) for k, v in log["cantilever_check"].items()})
    print("Uniform-bend check:", log["uniform_bend_check"])

    # ---- 3. run every case ----
    results, rows = [], []
    for label, p in CASES:
        hyd, raw = hydraulics(p, a.stage1)
        res = S.run_case(p, label, hyd=hyd, verbose=True)
        # Stage-1 checks EVERY node for the crossing; bend stations alone would detect it too late
        tcross = (raw["first_phase_crossing_s"] if raw is not None and raw.get("first_phase_crossing_s") is not None
                  else (None if raw is not None else first_crossing(hyd, PHASE_LIMIT_PA)))
        mag = np.hypot(res["u_dyn"][0::3], res["u_dyn"][1::3])
        if tcross is not None:                     # peak response BEFORE the phase limit is first crossed
            mag_valid = mag[:, res["t"] < tcross].max() * 1e3 if (res["t"] < tcross).any() else 0.0
        else:
            mag_valid = mag.max() * 1e3
        res["phase_cross_s"], res["surge_valid_mm"] = tcross, float(mag_valid)
        results.append(res)
        st = res["stress"]
        rows.append([label, p.t_close, p.R, p.H, p.W, res["therm_move_mm"], res["surge_mm"], mag_valid,
                     res["react_peak_kN"], res["react_static_kN"], st["vm_peak"] / 1e6,
                     tcross if tcross is not None else "none",
                     hyd.P.min() / 1e6])
        np.savetxt(os.path.join(a.out, f"history_{label.replace(' ', '_')}.csv"),
                   np.column_stack([res["t"], mag[res["k_node"]] * 1e3]), delimiter=",",
                   header="t_s,surge_mm_at_most_excited_node", comments="")
    head = ["case", "t_close_s", "R_m", "H_m", "W_m", "thermal_move_mm", "peak_surge_mm",
            "peak_surge_before_phase_crossing_mm", "anchor_peak_kN", "anchor_static_kN",
            "VM_peak_MPa", "first_phase_crossing_s", "min_pressure_MPa_at_bends"]
    with open(os.path.join(a.out, "comparison.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(head); w.writerows(rows)
    print("\n" + " | ".join(head[:1] + head[5:11]))
    for r in rows:
        print(f"{r[0]:22s} {r[5]:8.1f} {r[6]:8.1f} {r[7]:8.1f} {r[8]:8.0f} {r[9]:8.0f} {r[10]:8.1f}")

    # ---- 4. cross-check the two time integrators on the baseline ----
    hyd0, _ = hydraulics(CASES[0][1], a.stage1)
    rd = S.run_case(CASES[0][1], "baseline-direct", method="direct", hyd=hyd0, verbose=False)
    log["modal_vs_direct_surge_mm"] = [results[0]["surge_mm"], rd["surge_mm"]]
    print(f"\nModal vs direct Newmark: {results[0]['surge_mm']:.3f} vs {rd['surge_mm']:.3f} mm")

    # ---- 5. figures ----
    S.plot_pressure(results[0], os.path.join(a.out, "fig_pressure.png"))
    S.plot_modes(results[0], os.path.join(a.out, "fig_modes.png"))
    S.plot_spectra(results[0], os.path.join(a.out, "fig_spectra.png"))
    S.plot_response(results, os.path.join(a.out, "fig_response.png"))

    # ---- 6. convergence / sensitivity studies ----
    if not a.skip_sensitivity:
        log["sensitivity"] = {}
        for wh in ("damping", "hydraulic_dx", "structural_mesh", "time_step"):
            log["sensitivity"][wh] = S.sensitivity(BASE, wh)

    with open(os.path.join(a.out, "verification_log.json"), "w") as fh:
        json.dump(log, fh, indent=2, default=float)
    print(f"\nDone. Everything is in ./{a.out}/")


if __name__ == "__main__":
    main()
