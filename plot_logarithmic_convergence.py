"""Create optional base-10 plots from the researched convergence CSV."""

import argparse
import csv
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, LogFormatterMathtext, LogLocator, NullFormatter

from resolution import PEAK_PLOT_TITLE, TOLERANCE_PERCENT, label_refinement_errors


def logarithmic_x(axis):
    axis.set_xscale("log", base=10)
    axis.xaxis.set_major_locator(LogLocator(base=10))
    axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
    axis.xaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.set_xlabel("Number of reaches, N (logarithmic scale)")
    axis.grid(alpha=.3, which="major")
    axis.grid(alpha=.12, which="minor")


def run(source, output=None):
    output = source.parent if output is None else output
    with source.open(newline="") as file:
        rows = list(csv.DictReader(file))
    cases = {}
    for row in rows:
        cases.setdefault(row["case"], []).append(row)
    for values in cases.values():
        values.sort(key=lambda row: int(row["reaches"]))
    if not {"bulk_nominal", "elastic_nominal"}.issubset(cases):
        raise ValueError("Use the researched convergence CSV containing both nominal cases.")
    output.mkdir(parents=True, exist_ok=True)
    # A CSV stores physical values, not axis settings: preserve the original data.
    shutil.copyfile(source, output/"logarithmic_convergence.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), constrained_layout=True)
    for mode, axis in zip(("bulk", "elastic"), axes):
        for name, values in cases.items():
            if not name.startswith(mode+"_"):
                continue
            axis.plot([int(r["reaches"]) for r in values],
                      [float(r["first_rise_mpa"]) for r in values], "o-",
                      label=name.removeprefix(mode+"_").replace("_", " "))
        logarithmic_x(axis)
        axis.set(title=f"{mode.capitalize()} wave-speed case",
                 ylabel="First valve pressure rise [MPa]")
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
        axis.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.19), ncol=3)
    fig.suptitle(PEAK_PLOT_TITLE)
    fig.savefig(output/"logarithmic_convergence.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    fields = ("peak_difference_vs_finest_percent", "waveform_difference_vs_finest_percent")
    values = cases["bulk_nominal"]
    for axis, field in zip(axes, fields):
        # The finest grid has no independent difference from itself.
        compared = [r for r in values if r[field] != ""]
        axis.plot([int(r["reaches"]) for r in compared],
                  [float(r[field]) for r in compared], "o-", label="Bulk nominal")
    for axis in axes:
        logarithmic_x(axis)
        axis.set_yscale("log", base=10)
        axis.axhline(TOLERANCE_PERCENT, color="black", ls="--", label=f"{TOLERANCE_PERCENT:g}% criterion")
        axis.legend(fontsize=8)
    reference = values[-1]["reference_reaches"]
    label_refinement_errors(fig, axes, reference)
    fig.savefig(output/"logarithmic_convergence_error.png", dpi=160)
    plt.close(fig)
    print(f"Saved three logarithmic-convergence files in {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("results/selected/convergence.csv"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    run(args.input, args.out)
