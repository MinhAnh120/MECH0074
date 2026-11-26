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
from fluid_assumptions import PHASE_LABEL
from resolution import TOLERANCE_PERCENT, meets_tolerance, refinement_change
from run_stage1 import metrics, plot_valve_spectrum
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
            axis.axvspan(crossing, 1, color="#e6b96f", alpha=.18, label="Single-phase continuation after crossing")
        axis.set(xlim=(0, 1), title=f"Valve ({station:g} m)" if station == case.length else f"{station:g} m from reservoir")
        axis.grid(alpha=.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, ncol=4, loc="outside lower center")
    for axis in axes[-1]:
        axis.set_xlabel("Time after closure starts [s]")
    fig.supylabel("Absolute pressure [MPa]")
    fig.suptitle("100 m benchmark | c = 1100 m/s, frictionless")
    save(fig, output, "benchmark_pressure_overlay.png")

    fig, axis = plt.subplots(figsize=(8.5, 4), layout="constrained")
    t_return = theory(case)["round_trip_s"]
    axis.plot(time, run["pressure"][:, -1] / 1e6, color=COLORS["bulk"], lw=2, label="MOC")
    axis.plot(dense, exact[:, -1] / 1e6, "k--", lw=1, label="Analytical finite closure")
    axis.axvline(t_return, color="0.4", ls=":", label=f"First return at {t_return * 1e3:.2f} ms")
    axis.axhline(case.phase_limit / 1e6, color="#9b6000", ls="--", label=PHASE_LABEL)
    if crossing is not None:
        axis.axvspan(crossing, t_return + 1.5 * case.closure_time, color="#e6b96f", alpha=.18,
                     label="Single-phase continuation after crossing")
    axis.set(xlim=(t_return - .4 * case.closure_time, t_return + 1.5 * case.closure_time),
             ylim=(-.15, 1.55), xlabel="Time after closure starts [s]", ylabel="Valve absolute pressure [MPa]",
             title="First reservoir reflection at the valve")
    axis.legend(fontsize=8, loc="upper right")
    axis.grid(alpha=.25)
    save(fig, output, "valve_reflection_zoom.png")


def plot_valve_trace(axis, time, pressure, case, crossing, *,
                     crossing_label="Single-phase continuation after screen crossing ({crossing:.3f} s)"):
    keep = time <= 3
    axis.plot(time[keep], pressure[keep] / 1e6, color=COLORS["bulk"])
    axis.axhline(case.phase_limit / 1e6, color="#9b6000", ls="--", lw=1, label=PHASE_LABEL)
    if crossing is not None and crossing < 3:
        axis.axvspan(crossing, 3, color="#e6b96f", alpha=.18,
                    label=crossing_label.format(crossing=crossing))
    axis.set(xlim=(0, 3), xlabel="Time after closure starts [s]", ylabel="Valve absolute pressure [MPa]",
             title="First 3 s of the valve record")
    axis.legend(fontsize=8, loc="upper right")
    axis.grid(alpha=.22)


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
        crossing = long_record["first_phase_crossing_s"]
        fig, trace = plt.subplots(figsize=(9.5, 3.8), layout="constrained")
        plot_valve_trace(trace, time, pressure, case, crossing,
                         crossing_label="Single-phase continuation from crossing pressure at 111 K ({crossing:g} s)")
        fig.suptitle("100 m benchmark | c = 1100 m/s, nominal friction")
        save(fig, output, "valve_pressure_3s.png")

        fig, (trace, spectral) = plt.subplots(2, 1, figsize=(9.5, 6.8), layout="constrained")
        plot_valve_trace(trace, time, pressure, case, crossing)
        plot_valve_spectrum(spectral, case, frequency, amplitude)
        spectral.set_title(f"Hann-windowed spectrum | {len(time[official]) * dt:.2f} s analysed")
        spectral.grid(alpha=.22)
        fig.suptitle("100 m benchmark | c = 1100 m/s, nominal friction")
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
             title="Pressure waveform sensitivity to mesh size")
    fig.supxlabel(f"Reference mesh: N = {reference}. Maximum over 0–1 s at 10, 25, 50 and 100 m.",
                   fontsize=9)
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
    minimum.axhline(phase_limit / 1e6, color="#9b6000", ls="--", label=PHASE_LABEL)
    for axis in (peak, minimum):
        axis.set_xlabel("Closure time [s]")
        axis.set_xticks([.02, .10, .25, .50])
        axis.grid(alpha=.22)
        axis.legend(fontsize=8)
    fig.supylabel("Absolute pressure [MPa]")
    fig.suptitle("Closure sensitivity over the first 2 s | 100 m benchmark, nominal friction")
    save(fig, output, "closure_pressure_comparison.png")


def theory_figure(source, output, dt):
    with (source / "theory_comparison.csv").open(newline="") as file:
        rows = {r["quantity"]: r for r in csv.DictReader(file) if r["case"] == "bulk_ideal"}

    headline = (
        ("first_valve_rise", "First valve rise [MPa]", 1e-6, 5),
        ("reservoir_reflection_return_at_valve", "First valve return [ms]", 1e3, 3),
        ("pressure_period", "Pressure period [ms]", 1e3, 3),
        ("long_record_fft_peak", "Dominant FFT peak [Hz]", 1., 3),
    )
    values = [[label, *(f"{float(rows[key][column]) * scale:.{digits}f}"
                        for column in ("theory", "computed"))]
              for key, label, scale, digits in headline]
    arrivals = (
        ("first_arrival_x100", "Valve (closure starts)"),
        ("first_arrival_x50", "50 m station"),
        ("first_arrival_x25", "25 m station"),
        ("first_arrival_x10", "10 m station"),
        ("reservoir_one_way_arrival_using_flow", "Reservoir (from flow)"),
        ("reservoir_reflection_return_at_valve", "First return to valve"),
    )
    fig, (left, timing) = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    left.axis("off")
    table = left.table(cellText=values, colLabels=["Quantity", "Theory", "MOC / FFT"],
                       colWidths=[.49, .24, .27], loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.7)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor("#dbe1e7")
        cell.set_facecolor("#e8f1f8" if row == 0 else "white")
    timing.axvline(dt * 1e3, color="#d55e00", ls="--", label=f"One time step ({dt * 1e3:.3f} ms)")
    timing.scatter([float(rows[key]["computed_minus_theory"]) * 1e3 for key, _ in arrivals],
                   range(len(arrivals)), color=COLORS["bulk"], zorder=3, label="First detected sample")
    timing.set(yticks=range(len(arrivals)), yticklabels=[label for _, label in arrivals],
               xlim=(0, dt * 1e3 * 1.45), xlabel="Detected time minus theoretical arrival [ms]")
    timing.invert_yaxis()
    timing.grid(axis="x", alpha=.22)
    timing.legend(fontsize=8, loc="lower left")
    fig.suptitle("Stage 1 benchmark: theory versus computation | 100 m, c = 1100 m/s, frictionless")
    save(fig, output, "theory_comparison.png")


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
    theory_figure(source, output, report["cases"]["bulk_ideal"]["grids"][-1]["dt_s"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("results/latest"))
    parser.add_argument("--out", type=Path, default=Path("results/selected/evidence"))
    args = parser.parse_args()
    run(args.input, args.out)
