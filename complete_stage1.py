"""Run the researched Stage 1 verification, sensitivity and benchmark handoff."""

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import params as p
from fluid_assumptions import PHASE_BASIS, friction_estimate, methane_saturation_pa
from run_stage1 import metrics, save_case
from resolution import DEFAULT_REACHES, TOLERANCE_PERCENT, resolution_study, save_resolution, validate_reaches
from stage1 import Case, simulate, spectrum, theory


def shared_case(mode, friction):
    return Case(length=p.L_BENCHMARK, diameter=p.D_IN, outer_diameter=p.D_OUT,
                density=p.RHO_LNG, bulk_speed=p.C_BULK, young=p.E_STEEL,
                velocity=p.U0, pressure_abs=p.P0_ABS, external_abs=p.P_EXT_ABS,
                closure_time=p.T_CLOSE, speed_mode=mode, friction=friction,
                phase_limit=methane_saturation_pa(p.T_OPERATE))


def provenance(code_version):
    files = ("stage1.py", "run_stage1.py", "complete_stage1.py", "params.py",
             "fluid_assumptions.py", "export_bends.py", "resolution.py")
    return {"code_version": code_version, "python": platform.python_version(),
            "numpy": np.__version__, "source_sha256": {
                name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                for name in files}}


def write_csv(path, headings, rows):
    with path.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headings)
        writer.writerows(rows)


def comparison(case, short, measured, fft_peak):
    """Measure arrivals independently of their theoretical target times."""
    target, rows = theory(case), []
    def add(quantity, unit, expected, observed):
        rows.append([quantity, unit, expected, observed,
                     None if observed is None else observed-expected])
    add("first_valve_rise", "Pa", target["surge_pa"], measured["first_valve_rise_pa"])
    for arrival in measured["arrivals"]:
        add(f"first_arrival_x{arrival['station_m']:g}", "s",
            arrival["expected_s"], arrival["first_detected_s"])
    boundary = simulate(case, round(case.length/short["dx"]), 1, stations=[0, case.length])
    changed = np.flatnonzero(abs(boundary["flow"][:, 0]-boundary["flow"][0, 0])
                            > case.area*case.velocity*1e-6)
    add("reservoir_one_way_arrival_using_flow", "s", target["one_way_s"],
        float(boundary["time"][changed[0]]) if changed.size else None)
    # A frictionless closed-valve plateau ends at the first reflected wave.
    if case.friction == 0:
        falls = np.flatnonzero((np.diff(short["pressure"][:, -1]) < -target["surge_pa"]*1e-6)
                              & (short["time"][1:] > case.closure_time))
        add("reservoir_reflection_return_at_valve", "s", target["round_trip_s"],
            float(short["time"][falls[0]+1]) if falls.size else None)
    add("pressure_period", "s", target["period_s"], measured["period_s"])
    add("long_record_fft_peak", "Hz", target["fundamental_hz"], fft_peak)
    return rows


def run(output, code_version="uncommitted; use source_sha256", reaches=DEFAULT_REACHES):
    reaches = validate_reaches(reaches)
    finest = reaches[-1]
    output.mkdir(parents=True, exist_ok=True)
    estimate = friction_estimate(p.RHO_LNG, p.U0, p.D_IN)
    report = {**provenance(code_version), "phase_basis": PHASE_BASIS,
              "actual_LNG_composition_known": False,
              "phase_limit_pa": methane_saturation_pa(p.T_OPERATE),
              "friction_assumptions": estimate, "cases": {}}
    comparisons, sensitivities = [], []
    factors = {"ideal": 0.0, "nominal": estimate["darcy_factor"],
               "lower": estimate["lower"], "upper": estimate["upper"], "teaching_f002": 0.020}
    for mode in ("bulk", "elastic"):
        for label, factor in factors.items():
            name, case = f"{mode}_{label}", shared_case(mode, factor)
            short, rows, passing = resolution_study(case, reaches, metrics)
            start = math.ceil(case.closure_time/short["dt"])+4*finest
            samples = 4*finest*50
            long = simulate(case, finest, (start+samples)*short["dt"])
            section = slice(start, start+samples)
            frequency, amplitude, peak = spectrum(long["time"][section], long["pressure"][section, -1])
            item = {"parameters_SI": asdict(case), "analytical": theory(case), "grids": rows,
                    "reaches": list(reaches), "production_reaches": finest,
                    "tolerance_percent": TOLERANCE_PERCENT,
                    "smallest_tested_passing_reaches": passing,
                    "interpretation": "Formal single-phase verification; surrogate phase limit crossed"
                    if long["first_phase_crossing_s"] is not None else "Above surrogate threshold only; actual LNG unverified",
                    "long_record": {"start_s": float(long["time"][start]), "samples": samples,
                                    "DFT_duration_s": samples*short["dt"], "window": "Hann; mean removed",
                                    "bin_spacing_hz": float(frequency[1]), "peak_hz": peak,
                                    "first_phase_crossing_s": long["first_phase_crossing_s"],
                                    "minimum_absolute_pa": long["minimum_absolute_pa"],
                                    "phase_status": long["phase_status"]}}
            report["cases"][name] = item
            save_case(output, name, case, short, frequency, amplitude)
            if label in ("ideal", "nominal"):
                comparisons.extend([[name]+r for r in comparison(case, short, rows[-1], peak)])
                np.savetxt(output/f"{name}_long_valve.csv",
                           np.column_stack((long["time"], long["pressure"][:, -1], long["flow"][:, -1])),
                           delimiter=",", fmt="%.10g", header="time_s,p_abs_valve_pa,Q_valve_m3_s", comments="")
            sensitivities.append([name, factor, item["analytical"]["speed_m_s"],
                                  rows[-1]["first_valve_rise_pa"], short["minimum_absolute_pa"],
                                  rows[-1]["period_s"], peak, short["first_phase_crossing_s"],
                                  short["phase_status"]])
            print(f"{name}: f_D={factor:.8f}; first rise={rows[-1]['first_valve_rise_pa']:.3f} Pa; {item['interpretation']}")
    save_resolution(output, report["cases"])
    write_csv(output/"theory_comparison.csv", ["case", "quantity", "unit", "theory", "computed", "computed_minus_theory"], comparisons)
    write_csv(output/"friction_sensitivity.csv", ["case", "darcy_factor", "speed_m_s", "first_rise_pa", "min_absolute_pa", "period_s", "fft_hz", "first_phase_crossing_s", "phase_status"], sensitivities)
    closure_rows = []
    for mode in ("bulk", "elastic"):
        for closure in (.020, .10, .25, .50):
            case = replace(shared_case(mode, estimate["darcy_factor"]), closure_time=closure)
            duration = 50*theory(case)["period_s"] if closure == .50 else 2.0
            result = simulate(case, finest, duration)
            closure_rows.append([mode, closure, duration, float(result["pressure"][:, -1].max()),
                                 result["minimum_absolute_pa"], result["first_phase_crossing_s"], result["phase_status"]])
    write_csv(output/"closure_sensitivity.csv", ["speed_mode", "closure_s", "record_duration_s", "max_valve_absolute_pa", "global_min_absolute_pa", "first_phase_crossing_s", "phase_status"], closure_rows)
    (output/"verification.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/latest"))
    parser.add_argument("--reaches", type=int, nargs="+", default=list(DEFAULT_REACHES))
    parser.add_argument("--code-version", default="uncommitted; use source_sha256")
    args = parser.parse_args()
    run(args.out, args.code_version, args.reaches)
