"""Run all Stage 2 scripts in order.

Usage (from any folder):
    python run_stage2.py

Runs, in this order:
    1. stage2_euler_bernoulli.py   analytical hand calculation + mode shapes
    2. stage2_cantilever.py        beam FEM checks (flow, mesh, Timoshenko)
    3. stage2_results.py           comparison + mesh-study tables and figures

Each script's printed output is shown on screen and also saved to
stage2_logs/<script>.log. All outputs are written next to this file.
Stops at the first script that fails.
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

REQUIRED = ["params.py", "stage2_euler_bernoulli.py", "stage2_cantilever.py", "stage2_results.py"]
SEQUENCE = ["stage2_euler_bernoulli.py", "stage2_cantilever.py", "stage2_results.py"]
EXPECTED_OUTPUTS = [
    "stage2_EB_mode_shapes.png",
    "stage2_table_comparison.csv",
    "stage2_table_mesh.csv",
    "stage2_fig_comparison.png",
    "stage2_fig_mesh.png",
]


def check_setup():
    missing = [f for f in REQUIRED if not (HERE / f).exists()]
    if missing:
        print("Missing files in", HERE)
        for f in missing:
            print("   -", f)
        print("Put all Stage 2 files in the same folder as run_stage2.py and try again.")
        sys.exit(1)
    try:
        import numpy, scipy, matplotlib  # noqa: F401
    except ImportError as err:
        print(f"Missing package: {err.name}. Run:  pip install -r requirements.txt")
        sys.exit(1)


def run(script, log_dir):
    print("\n" + "=" * 78)
    print(f" Running {script}")
    print("=" * 78)
    t0 = time.time()
    result = subprocess.run([sys.executable, script], cwd=HERE,
                            capture_output=True, text=True)
    output = result.stdout + result.stderr
    print(output)
    (log_dir / f"{Path(script).stem}.log").write_text(output, encoding="utf-8")
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"FAILED: {script} (exit code {result.returncode}) after {dt:.1f} s")
        print(f"See {log_dir / (Path(script).stem + '.log')}")
        sys.exit(result.returncode)
    print(f"Finished {script} in {dt:.1f} s")


def main():
    check_setup()
    log_dir = HERE / "stage2_logs"
    log_dir.mkdir(exist_ok=True)
    for script in SEQUENCE:
        run(script, log_dir)

    print("\n" + "=" * 78)
    print(" Stage 2 complete - outputs in", HERE)
    print("=" * 78)
    for f in EXPECTED_OUTPUTS:
        status = "OK     " if (HERE / f).exists() else "MISSING"
        print(f"  [{status}] {f}")
    print(f"  Logs: {log_dir}")


if __name__ == "__main__":
    main()
