"""Generate Stage 1 figures and verification tables from a full benchmark run."""

import argparse
import csv
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from complete_stage1 import comparison, write_csv
from resolution import TOLERANCE_PERCENT, meets_tolerance, refinement_change
from run_stage1 import metrics
from stage1 import Case, estimate_period, simulate, spectrum, theory

MODES = ("bulk", "elastic")
COLORS = {"bulk": "#0072b2", "elastic": "#d55e00"}


def save(fig, output, filename):
    fig.savefig(output / filename, dpi=180)
    plt.close(fig)


def analytical_pressure(case, time, stations):
    """Superpose delayed finite ramps for fixed reservoir and closed valve."""
    target = theory(case)
    time, stations = time[:, None], stations[None, :]
    pressure = np.full((time.size, stations.size), case.pressure_abs)
    for m in range(int(time.max() / target["round_trip_s"]) + 1):
        delay = (2 * m + 1) * target["one_way_s"]
        incoming = np.clip((time - delay + stations / case.speed) / case.closure_time, 0, 1)
        reflected = np.clip((time - delay - stations / case.speed) / case.closure_time, 0, 1)
        pressure += (-1)**m * target["surge_pa"] * (incoming - reflected)
    return pressure


def pressure_figures(case, run, output):
    time, stations = run["time"], run["stations"]
    dense = np.linspace(0, 1, 12001)
    exact = analytical_pressure(case, dense, stations)
    crossing = run["first_phase_crossing_s"]
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(10, 6), layout="constrained")
    for i, axis in enumerate(axes.flat):
        station = stations[i]
        arrival = (case.length - station) / case.speed
        axis.plot(time, run["pressure"][:, i] / 1e6, color=f"C{i}", lw=2, label="MOC")
        axis.plot(dense, exact[:, i] / 1e6, "k--", lw=1.1, label="Analytical finite closure")
        axis.axvline(arrival, color="0.45", ls=":", lw=1, label="First arrival")
        if crossing is not None:
            axis.axvspan(crossing, 1, color="#e6b96f", alpha=.18, label="After first phase-screen crossing")
        axis.set(xlim=(0, 1), title=f"Valve ({station:g} m)" if station == case.length else f"{station:g} m from reservoir")
        axis.grid(alpha=.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, ncol=4, loc="outside lower center")
    for axis in axes[-1]:
        axis.set_xlabel("Time after closure starts [s]")
    fig.supylabel("Absolute pressure [MPa]")
    fig.suptitle("100 m benchmark | bulk wave speed, frictionless\nFormal single-phase verification")
    save(fig, output, "benchmark_pressure_overlay.png")

    fig, axis = plt.subplots(figsize=(8.5, 4), layout="constrained")
    t_return = theory(case)["round_trip_s"]
    axis.plot(time, run["pressure"][:, -1] / 1e6, color=COLORS["bulk"], lw=2, label="MOC")
    axis.plot(dense, exact[:, -1] / 1e6, "k--", lw=1, label="Analytical finite closure")
    axis.axvline(t_return, color="0.4", ls=":", label=f"First return at {t_return * 1e3:.2f} ms")
    axis.axhline(case.phase_limit / 1e6, color="#9b6000", ls="--", label="Methane surrogate phase threshold")
    if crossing is not None:
        axis.axvspan(crossing, t_return + 1.5 * case.closure_time, color="#e6b96f", alpha=.18,
                     label="After first phase-screen crossing")
    axis.set(xlim=(t_return - .4 * case.closure_time, t_return + 1.5 * case.closure_time),
             ylim=(-.15, 1.55), xlabel="Time after closure starts [s]", ylabel="Valve absolute pressure [MPa]",
             title="First reservoir reflection at the valve | formal single-phase verification")
    axis.legend(fontsize=8, loc="upper right")
    axis.grid(alpha=.25)
    save(fig, output, "valve_reflection_zoom.png")


def valve_frequency(case, item, source, output):
    mode = case.speed_mode
    saved_time, pressure = np.loadtxt(source / f"{mode}_nominal_long_valve.csv",
                                     delimiter=",", skiprows=1, usecols=(0, 1), unpack=True)
    dt = item["grids"][-1]["dt_s"]
    time = np.arange(saved_time.size) * dt
    # Reconstruct the uniform time vector from metadata after checking CSV rounding.
    if time[-1] < 17.5 or not np.allclose(time, saved_time, rtol=6e-10, atol=1e-12):
        raise ValueError("Expected the full, uniformly sampled nominal valve record.")
    long_record = item["long_record"]
    start = round(long_record["start_s"] / dt)
    official = slice(start, start + long_record["samples"])
    selected = (time >= .5) & (time < 17.5)
    frequency, amplitude, peak = spectrum(time[official], pressure[official])
    freq_check, amp_check, raw_peak = spectrum(time[selected], pressure[selected])
    i = np.argmax(amp_check[1:]) + 1
    left, centre, right = np.log(amp_check[i-1:i+2])
    # Interpolate the dominant peak using three adjacent log-amplitude bins.
    refined_peak = raw_peak + .5 * (left - right) / (left - 2 * centre + right) * freq_check[1]
    period = estimate_period(time[selected], pressure[selected] - case.pressure_abs, .5)
    target = theory(case)["fundamental_hz"]
    rows = [[mode, target, time[official][0], len(time[official]), len(time[official]) * dt,
             frequency[1], peak, time[selected][0], int(selected.sum()), selected.sum() * dt,
             freq_check[1], raw_peak, refined_peak, period, 1 / period]]
    harmonics = []
    for order in (1, 3, 5, 7):
        j = np.argmin(abs(frequency - order * target))
        harmonics.append([mode, order, order * target, frequency[j], amplitude[j]])
    if mode == "bulk":
        fig, (trace, spectral) = plt.subplots(2, 1, figsize=(9.5, 6.8), layout="constrained")
        keep = time <= 3 + dt
        trace.plot(time[keep], pressure[keep] / 1e6, color=COLORS[mode])
        trace.axhline(case.phase_limit / 1e6, color="#9b6000", ls="--", lw=1, label="Methane surrogate threshold")
        crossing = long_record["first_phase_crossing_s"]
        if crossing is not None and crossing < 3:
            trace.axvspan(crossing, 3, color="#e6b96f", alpha=.18, label=f"After first crossing ({crossing:.3f} s)")
        trace.set(xlim=(0, 3), xlabel="Time after closure starts [s]", ylabel="Valve absolute pressure [MPa]",
                  title="First 3 s of the valve record")
        trace.legend(fontsize=8, loc="upper right")
        spectral.plot(frequency, amplitude / 1e6, color=COLORS[mode])
        for _, order, harmonic_hz, _, _ in harmonics:
            spectral.axvline(harmonic_hz, color="0.45", ls=":", lw=1)
            spectral.text(harmonic_hz + (.4 if order == 1 else 0), .90 if order == 1 else .78,
                          f"{harmonic_hz:.2f} Hz", transform=spectral.get_xaxis_transform(),
                          ha="left" if order == 1 else "center", fontsize=8)
        spectral.set(xlim=(0, 20), xlabel="Frequency [Hz]", ylabel="Valve pressure amplitude [MPa]",
                     title=f"Hann-windowed spectrum | {len(time[official]) * dt:.2f} s analysed")
        for axis in (trace, spectral):
            axis.grid(alpha=.22)
        fig.suptitle("Bulk benchmark with nominal friction | formal single-phase verification")
        save(fig, output, "valve_record_and_spectrum.png")
    return rows, harmonics


def waveform_figure(source, output):
    with (source / "convergence.csv").open(newline="") as file:
        rows = list(csv.DictReader(file))
    fig, axis = plt.subplots(figsize=(8, 4.6), layout="constrained")
    for mode in MODES:
        data = [r for r in rows if r["case"] == f"{mode}_nominal" and r["waveform_difference_vs_finest_percent"]]
        axis.loglog([int(r["reaches"]) for r in data],
                    [float(r["waveform_difference_vs_finest_percent"]) for r in data],
                    "o-", color=COLORS[mode], label=mode.capitalize())
    reference = data[0]["reference_reaches"]
    axis.axhline(TOLERANCE_PERCENT, color="k", ls="--", label=f"Chosen {TOLERANCE_PERCENT:g}% criterion")
    axis.set(xlabel="Number of reaches, N [log scale]", ylabel="Maximum difference / Joukowsky rise [%]",
             title=f"Pressure histories compared with N = {reference}\nFirst second, four stations, nominal friction")
    axis.legend()
    axis.grid(which="both", alpha=.22)
    save(fig, output, "waveform_convergence.png")


def closure_figure(rows, phase_limit, output):
    fig, (peak, minimum) = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True, layout="constrained")
    for mode in MODES:
        data = [r for r in rows if r[0] == mode]
        for axis, column in ((peak, 4), (minimum, 5)):
            axis.plot([r[1] for r in data], [r[column] / 1e6 for r in data], "o-",
                      color=COLORS[mode], label=mode.capitalize())
    peak.set_title("Maximum valve pressure")
    minimum.set_title("Minimum pressure anywhere in the pipe")
    minimum.axhline(phase_limit / 1e6, color="#9b6000", ls="--", label="Methane surrogate threshold")
    for axis in (peak, minimum):
        axis.set_xlabel("Closure time [s]")
        axis.set_xticks([.02, .10, .25, .50])
        axis.grid(alpha=.22)
        axis.legend(fontsize=8)
    fig.supylabel("Absolute pressure [MPa]")
    fig.suptitle("Closure sensitivity over the first 2 s | 100 m benchmark, nominal friction")
    save(fig, output, "closure_pressure_comparison.png")


def run(source, output):
    report = json.loads((source / "verification.json").read_text())
    for name, expected in report["source_sha256"].items():
        source_bytes = Path(__file__).with_name(name).read_bytes()
        if hashlib.sha256(source_bytes).hexdigest() != expected:
            raise ValueError(f"Solver file {name} differs from the report. Regenerate the source run.")
    output.mkdir(parents=True, exist_ok=True)
    benchmark, frequency, harmonics, refinement, closure = [], [], [], [], []
    for mode in MODES:
        ideal_item = report["cases"][f"{mode}_ideal"]
        nominal_item = report["cases"][f"{mode}_nominal"]
        reaches = nominal_item["production_reaches"]
        ideal = Case(**ideal_item["parameters_SI"])
        nominal = Case(**nominal_item["parameters_SI"])
        short = simulate(ideal, reaches, 1.)
        observed = metrics(ideal, short)
        freq = ideal_item["long_record"]["peak_hz"]
        benchmark.extend([[mode] + r for r in comparison(ideal, short, observed, freq)])
        rise = short["pressure"][short["time"] < theory(ideal)["round_trip_s"], 0].max() - ideal.pressure_abs
        incident_to_reflected = 2 * short["stations"][0] / ideal.speed
        exact = theory(ideal)["surge_pa"] * min(1., incident_to_reflected / ideal.closure_time)
        benchmark.append([mode, "first_pressure_rise_x10", "Pa", exact, rise, rise - exact])
        if mode == "bulk":
            pressure_figures(ideal, short, output)
        estimates, bands = valve_frequency(nominal, nominal_item, source, output)
        frequency.extend(estimates)
        harmonics.extend(bands)
        coarse, fine = (simulate(nominal, count, 1.) for count in (reaches, 2 * reaches))
        changes = refinement_change(coarse, fine, metrics(nominal, coarse), metrics(nominal, fine),
                                    theory(nominal)["surge_pa"])
        refinement.append([mode, reaches, 2 * reaches, changes["first_valve_rise_pa_change_percent"],
                           changes["period_s_change_percent"], changes["waveform_difference_percent"],
                           TOLERANCE_PERCENT, meets_tolerance(changes)])
        for close in (.020, .10, .25, .50):
            scenario = replace(nominal, closure_time=close)
            dt = (scenario.length / reaches) / scenario.speed
            # The solver rounds up to a full step. End at the last sample at or before 2 s.
            duration = (math.floor(2 / dt) - .5) * dt
            result = simulate(scenario, reaches, duration)
            peak = result["pressure"][:, -1].max()
            closure.append([mode, close, close / theory(scenario)["round_trip_s"], result["time"][-1], peak,
                            result["minimum_absolute_pa"], result["first_phase_crossing_s"], result["phase_status"]])
    write_csv(output / "benchmark_checks.csv",
              ["speed_mode", "quantity", "unit", "theory", "computed", "computed_minus_theory"], benchmark)
    write_csv(output / "frequency_checks.csv",
              ["speed_mode", "theory_hz", "main_start_s", "main_samples", "main_dft_duration_s", "main_bin_spacing_hz",
               "main_peak_hz", "check_start_s", "check_samples", "check_dft_duration_s", "check_bin_spacing_hz",
               "check_raw_peak_hz", "check_interpolated_peak_hz", "check_period_s", "check_crossing_frequency_hz"], frequency)
    write_csv(output / "pressure_harmonics.csv", ["speed_mode", "order", "theory_hz", "fft_bin_hz", "amplitude_pa"], harmonics)
    write_csv(output / "fine_grid_check.csv",
              ["speed_mode", "production_n", "finer_n", "first_rise_difference_pct", "period_difference_pct",
               "waveform_difference_pct_of_joukowsky", "criterion_pct", "passes"], refinement)
    write_csv(output / "closure_pressure_comparison.csv",
              ["speed_mode", "closure_s", "closure_to_return_ratio", "simulation_end_s", "valve_max_absolute_pa",
               "pipe_min_absolute_pa", "first_phase_crossing_s", "phase_status"], closure)
    waveform_figure(source, output)
    closure_figure(closure, report["phase_limit_pa"], output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("results/latest"))
    parser.add_argument("--out", type=Path, default=Path("results/selected/evidence"))
    args = parser.parse_args()
    run(args.input, args.out)
