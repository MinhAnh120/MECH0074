"""Run the three labelled Stage 1 teaching cases and write reproducible outputs."""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
from stage1 import Case, estimate_period, simulate, spectrum, theory
from resolution import DEFAULT_REACHES, TOLERANCE_PERCENT, refinement_change, resolution_study, save_resolution, validate_reaches


def metrics(case, result):
    target = theory(case)
    t, pressure = result["time"], result["pressure"]
    delta = pressure-pressure[0]
    first = t < target["round_trip_s"]
    rise = float(delta[first, -1].max())
    arrivals = []
    for i, x in enumerate(result["stations"]):
        found = np.flatnonzero(abs(delta[:, i]) > target["surge_pa"]*1e-6)
        observed = float(t[found[0]]) if found.size else None
        expected = (case.length-float(x))/case.speed
        arrivals.append({"station_m": float(x), "expected_s": expected,
                         "first_detected_s": observed,
                         "error_s": None if observed is None else observed-expected})
    period = estimate_period(t, pressure[:, -1]-case.pressure_abs, case.closure_time)
    return {"dx_m": result["dx"], "dt_s": result["dt"],
            "first_valve_rise_pa": rise,
            "surge_error_percent": 100*(rise/target["surge_pa"]-1),
            "period_s": period,
            "period_error_percent": None if period is None else 100*(period/target["period_s"]-1),
            "arrivals": arrivals, "minimum_absolute_pa": result["minimum_absolute_pa"],
            "first_nonpositive_s": result["first_nonpositive_s"],
            "first_phase_crossing_s": result["first_phase_crossing_s"],
            "phase_status": result["phase_status"],
            "interpolated_stations": result["interpolated_stations"]}


def save_case(directory, name, case, short, frequency, amplitude):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, time, pressure = short["stations"], short["time"], short["pressure"]
    columns = ["time_s"]+[f"p_abs_x{v:g}_m_pa" for v in x]+[f"Q_x{v:g}_m_m3_s" for v in x]
    np.savetxt(directory/f"{name}_histories.csv",
               np.column_stack((time, pressure, short["flow"])), delimiter=",",
               header=",".join(columns), comments="", fmt="%.12g")
    np.savetxt(directory/f"{name}_spectrum.csv", np.column_stack((frequency, amplitude)),
               delimiter=",", header="frequency_hz,pressure_amplitude_pa", comments="")
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), constrained_layout=True)
    for i, station in enumerate(x):
        label = f"x = {station:g} m"+(" (valve)" if station == case.length else "")
        axes[0].plot(time, pressure[:, i]/1e6, label=label)
        axes[1].plot(time, (pressure[:, i]-pressure[0, i])/1e6, label=label)
    axes[0].set(xlabel="Time [s]", ylabel="Absolute pressure [MPa]", xlim=(0, 1))
    axes[1].set(xlabel="Time [s]", ylabel="Pressure change [MPa]", xlim=(0, 1))
    if case.phase_limit is not None:
        axes[0].axhline(case.phase_limit/1e6, color="black", ls=":",
                        label="Entered phase threshold (see metadata)")
    axes[0].legend(ncol=3, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    axes[2].plot(frequency, amplitude/1e6)
    axes[2].axvline(case.speed/(4*case.length), color="tab:red", ls="--", label="a/(4L)")
    axes[2].set(title="Valve-pressure spectrum from the longer record",
                xlabel="Frequency [Hz]", ylabel="Valve amplitude [MPa]", xlim=(0, 20))
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.3)
    validity = {"LNG_phase_limit_not_supplied": "LNG phase limit unknown; physical validity unverified",
                "nonpositive_absolute_pressure": "Negative absolute pressure; formal verification only",
                "phase_limit_crossed": "Entered phase limit crossed; formal verification only",
                "above_entered_phase_limit_only": "Above entered phase limit; other model limits still apply"}
    title = name.replace("_", " ").capitalize()
    fig.suptitle(f"{title}: a = {case.speed:.3f} m/s; Darcy f = {case.friction:g}\n"
                 f"{validity[short['phase_status']]}", fontsize=10)
    fig.savefig(directory/f"{name}.png", dpi=160)
    plt.close(fig)


def run(output, reaches, periods, friction, phase_limit):
    import matplotlib

    reaches = validate_reaches(reaches)
    if periods < 5:
        raise ValueError("Use at least five periods for the longer FFT record.")
    output.mkdir(parents=True, exist_ok=True)
    report = {"python": platform.python_version(), "numpy": np.__version__,
              "matplotlib": matplotlib.__version__, "sources": {}, "cases": {}}
    for filename in ("stage1.py", "run_stage1.py", "resolution.py"):
        report["sources"][filename] = hashlib.sha256(Path(__file__).with_name(filename).read_bytes()).hexdigest()
    cases = {"elastic_ideal": Case(phase_limit=phase_limit),
             "bulk_ideal": Case(speed_mode="bulk", phase_limit=phase_limit),
             "elastic_friction_sensitivity": Case(friction=friction, phase_limit=phase_limit)}
    for name, case in cases.items():
        short, rows, passing = resolution_study(case, reaches, metrics)
        finest = reaches[-1]
        # A fixed 50-period default window has the same duration on every grid.
        start = math.ceil(case.closure_time/short["dt"])+4*finest
        samples = 4*finest*periods
        long = simulate(case, finest, (start+samples)*short["dt"])
        section = slice(start, start+samples)
        frequency, amplitude, peak = spectrum(long["time"][section], long["pressure"][section, -1])
        item = {"parameters_SI": asdict(case), "analytical": theory(case), "grids": rows,
                "reaches": list(reaches), "production_reaches": finest,
                "tolerance_percent": TOLERANCE_PERCENT,
                "smallest_tested_passing_reaches": passing,
                "long_record": {"samples": samples, "window": "Hann; mean removed",
                                "start_s": float(long["time"][start]),
                                "DFT_duration_s": samples*long["dt"],
                                "bin_spacing_hz": float(frequency[1]), "peak_hz": peak,
                                "minimum_absolute_pa": long["minimum_absolute_pa"],
                                "phase_status": long["phase_status"]}}
        report["cases"][name] = item
        save_case(output, name, case, short, frequency, amplitude)
        print(f"{name}: surge={rows[-1]['first_valve_rise_pa']/1e6:.6f} MPa, "
              f"FFT={peak:.6f} Hz; {long['phase_status']}")
    (output/"verification.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    save_resolution(output, report["cases"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/latest"))
    parser.add_argument("--reaches", type=int, nargs="+", default=list(DEFAULT_REACHES))
    parser.add_argument("--periods", type=int, default=50)
    parser.add_argument("--friction", type=float, default=0.020,
                        help="Assumed Darcy factor for the sensitivity case only.")
    parser.add_argument("--phase-limit-pa", type=float, default=None,
                        help="Optional justified LNG absolute phase pressure; no default is invented.")
    args = parser.parse_args()
    run(args.out, args.reaches, args.periods, args.friction, args.phase_limit_pa)
