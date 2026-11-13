"""Shared hydraulic refinement study; numerical tolerance is a group choice."""

import csv
import math
import numpy as np
from stage1 import simulate, theory

DEFAULT_REACHES = (10, 20, 40, 80, 160, 320)
TOLERANCE_PERCENT = 1.0


def validate_reaches(reaches):
    values = tuple(reaches)
    if (len(values) < 3
            or any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1
                   for n in values)
            or any(b != 2*a for a, b in zip(values, values[1:]))):
        raise ValueError("Supply at least three positive, successively doubled integer reach counts.")
    return values


def refinement_change(coarse, fine, coarse_metrics, fine_metrics, surge):
    """Reconstruct coarse histories on fine times to detect between-sample errors."""
    if not math.isfinite(surge) or surge <= 0:
        raise ValueError("Waveform normalisation needs a positive reference surge.")
    if not np.array_equal(coarse["stations"], fine["stations"]):
        raise ValueError("Compare the same physical recording stations.")
    keep = fine["time"] <= min(1.0, coarse["time"][-1])
    t = fine["time"][keep]
    difference = max(np.max(abs(np.interp(t, coarse["time"], coarse["pressure"][:, j])
                                - fine["pressure"][keep, j]))
                     for j in range(fine["stations"].size))
    changes = {"waveform_difference_percent": float(100*difference/surge)}
    for key in ("first_valve_rise_pa", "period_s"):
        a, b = coarse_metrics[key], fine_metrics[key]
        changes[key+"_change_percent"] = (
            None if a is None or b is None or b == 0 else 100*abs(b-a)/abs(b))
    return changes


def meets_tolerance(changes, tolerance=TOLERANCE_PERCENT):
    return changes is not None and all(
        value is not None and math.isfinite(value) and value < tolerance
        for value in changes.values())


def resolution_study(case, reaches, measure):
    reaches = validate_reaches(reaches)
    histories = [simulate(case, n, 1.0) for n in reaches]
    rows = [{"reaches": n, **measure(case, result)}
            for n, result in zip(reaches, histories)]
    for i, row in enumerate(rows):
        if i:
            row["change_from_previous"] = refinement_change(
                histories[i-1], histories[i], rows[i-1], row, theory(case)["surge_pa"])
        if i < len(rows)-1:
            row["change_to_next"] = refinement_change(
                histories[i], histories[i+1], row, rows[i+1], theory(case)["surge_pa"])
            row["difference_from_finest"] = refinement_change(
                histories[i], histories[-1], row, rows[-1], theory(case)["surge_pa"])
            row["meets_criterion_vs_next"] = meets_tolerance(row["change_to_next"])
            row["meets_criterion_vs_finest"] = meets_tolerance(row["difference_from_finest"])
        else:
            # The reference cannot demonstrate convergence by comparison with itself.
            row.update(change_to_next=None, difference_from_finest=None,
                       meets_criterion_vs_next=None, meets_criterion_vs_finest=None)
    passing = next((row["reaches"] for i, row in enumerate(rows[:-1])
                    if all(r["meets_criterion_vs_next"] and r["meets_criterion_vs_finest"]
                           for r in rows[i:-1])), None)
    return histories[-1], rows, passing


def save_resolution(output, cases):
    """Save the requested peak plot, error plot and machine-readable diagnostics."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    headings = ["case", "reaches", "dt_s", "first_rise_pa", "period_s", "surge_error_percent",
                "first_rise_mpa", "interpolated_stations", "reference_reaches",
                "peak_change_vs_next_percent", "period_change_vs_next_percent",
                "waveform_change_vs_next_percent", "meets_1pct_vs_next",
                "peak_difference_vs_finest_percent", "period_difference_vs_finest_percent",
                "waveform_difference_vs_finest_percent", "meets_1pct_vs_finest"]
    keys = ("first_valve_rise_pa_change_percent", "period_s_change_percent",
            "waveform_difference_percent")
    with (output/"convergence.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headings)
        for name, item in cases.items():
            for row in item["grids"]:
                following, finest = row["change_to_next"], row["difference_from_finest"]
                writer.writerow([name, row["reaches"], row["dt_s"], row["first_valve_rise_pa"],
                                 row["period_s"], row["surge_error_percent"],
                                 row["first_valve_rise_pa"]/1e6, row["interpolated_stations"],
                                 item["grids"][-1]["reaches"],
                                 *[None if following is None else following[k] for k in keys],
                                 row["meets_criterion_vs_next"],
                                 *[None if finest is None else finest[k] for k in keys],
                                 row["meets_criterion_vs_finest"]])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for mode, axis in zip(("bulk", "elastic"), axes):
        for name, item in cases.items():
            if item["parameters_SI"]["speed_mode"] != mode:
                continue
            rows = item["grids"]
            label = name.removeprefix(mode+"_").replace("_", " ")
            axis.plot([r["reaches"] for r in rows],
                      [r["first_valve_rise_pa"]/1e6 for r in rows], "o-", label=label)
        axis.set(title=f"{mode.capitalize()} wave-speed case", xlabel="Number of reaches, N",
                 ylabel="First valve pressure rise [MPa]")
        axis.set_xscale("log", base=2)
        axis.set_xticks(next(iter(cases.values()))["reaches"],
                        [str(n) for n in next(iter(cases.values()))["reaches"]])
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
        axis.grid(alpha=.3)
        axis.legend(fontsize=8, loc="lower right")
    fig.suptitle("Peak pressure rise versus reaches (horizontal scale: successive doubling)")
    fig.savefig(output/"convergence.png", dpi=160)
    plt.close(fig)

    selected = {n: c for n, c in cases.items() if n.endswith("_nominal")}
    if not selected:
        selected = {n: c for n, c in cases.items() if c["parameters_SI"]["friction"] > 0}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for name, item in selected.items():
        rows = item["grids"][:-1]
        for axis, key in zip(axes, (keys[0], keys[2])):
            values = [r["difference_from_finest"][key] for r in rows]
            axis.plot([r["reaches"] for r in rows], values, "o-", label=name)
    for axis, title in zip(axes, ("First-rise difference", "Maximum waveform difference")):
        axis.axhline(TOLERANCE_PERCENT, color="black", ls="--", label="1% criterion")
        axis.set(xlabel="Number of reaches, N", ylabel="Difference from finest grid [%]", title=title)
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        ticks = next(iter(cases.values()))["reaches"][:-1]
        axis.set_xticks(ticks, [str(n) for n in ticks])
        axis.grid(alpha=.3, which="both")
        axis.legend(fontsize=8)
    reference = next(iter(cases.values()))["reaches"][-1]
    fig.suptitle(f"Refinement diagnostics; N={reference} is the numerical reference")
    fig.savefig(output/"convergence_error.png", dpi=160)
    plt.close(fig)
