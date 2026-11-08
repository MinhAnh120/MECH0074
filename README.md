# MECH0074: LNG pressure waves

MAN's UCL MECH0074 Engineering in Extreme Environments coursework, 2026/27.
Development branch: **stage1-pressure-waves** in **MinhAnh120/MECH0074**.

## Contents
- [Run the model](#run-the-model)
- [Cases and assumptions](#cases-and-assumptions)
- [Outputs](#outputs)
- [Verification](#verification)
- [Files and provenance](#files-and-provenance)
- [Continue in Codex](#continue-in-codex)

## Run the model
Tested with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8.
From a terminal with Python 3.12 available:

```sh
git clone --branch stage1-pressure-waves https://github.com/MinhAnh120/MECH0074.git
cd MECH0074
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_stage1.py --out results/latest
```

The activation command is for Linux/macOS. Alternatively run the Python commands in an environment containing the pinned dependencies. No SciPy, Jupyter or LAMMPS dependency is needed for this implementation.

```sh
python run_stage1.py --reaches 100 200 400 --periods 50 --friction 0.020 --out results/run-custom
```

`--phase-limit-pa` optionally accepts a justified LNG absolute phase pressure. No value is supplied by default. Record its composition/temperature/source in PROJECT_CONTEXT.md before interpreting it physically.

## Cases and assumptions
All runs use the issued 100 m benchmark, Di=0.79394 m, Do=0.813 m, rho=455 kg/m3, u0=1.5 m/s and prescribed 0.020 s linear velocity/flow closure.

| Case | Wave speed | Darcy friction | Purpose |
|---|---|---|---|
| elastic_ideal | 997.899867 m/s from lecture 6 p.39 | 0 | Elastic-wall analytical verification |
| bulk_ideal | 1100 m/s supplied bulk speed | 0 | Explicit fixed-area alternative |
| elastic_friction_sensitivity | Lecture elastic speed | Assumed 0.020 by default | Quasi-steady friction sensitivity |

The brief does not unambiguously settle the bulk/effective-speed interpretation. These cases retain the distinction; the lecturer's intended reference remains to be confirmed. E=213.233 GPa and wall thickness=0.00953 m enter the elastic speed calculation. Poisson/axial-motion refinements are outside the simplified lecture formula.

The pipe is horizontal with z=0; the reservoir boundary is assigned 0.7 MPa absolute. H is pressure head above the fixed external pressure 101325 Pa: p_abs=p_external+rho*g*H. The frictional initial H profile includes steady Darcy loss, so initial valve pressure is slightly below the reservoir pressure. The brief does not locate the initial pressure in a frictional profile. Friction 0.020 is taken only as a teaching assumption, not as LNG property data.

The solver adapts notebook cells6-7: fixed head reservoir, prescribed-flow valve, old-to-new characteristics and signed explicit friction. All output stations are measured FROM the reservoir. Benchmark stations lie exactly on all default grids. Other stations use explicitly reported linear interpolation; the solver never silently changes wave speed.

**Physical limits:** constant-property, single-phase, quasi-steady-friction hydraulics; no cavitation, bend-reflection or moving-pipe feedback model. The bulk case reaches negative absolute pressure. Elastic cases also require an LNG phase-pressure limit before physical validity can be established. Curves beyond a phase limit are formal numerical continuations, not approved loads for Beau.

## Outputs
`results/latest/` is generated locally and ignored by Git. For each case:
- `*_histories.csv`: time, absolute pressure and flow at 10,25,50,100 m; the unaltered CFL grid extends to at least 1 s.
- `*_spectrum.csv`: valve single-sided Hann-windowed pressure amplitude (Pa) against frequency (Hz).
- `*.png`: 0-1 s absolute and incremental pressures plus longer-record FFT.
- `verification.json`: all parameters, source-code hashes, versions, analytical values, per-grid metrics, refinement differences, timing/phase diagnostics and FFT window details.
- `convergence.csv`: first valve rise and period versus reach count.

`results/selected/` contains a reproducible snapshot of selected numerical verification from the initial implementation. Rerunning generates all histories/spectra, not just these selected files. A benchmark station CSV is not an actual loop bend CSV.

The longer FFT uses 50 ideal periods by default, starting one ideal period plus the closure ramp after t=0. For elastic cases the DFT window is 20.042091 s; bulk is 18.181818 s. M=4*N*periods samples; duration=M*dt, not (M-1)*dt. The fundamental lies at bin 50 in the default ideal cases. Frequency recovery must also agree with the independently measured time-domain period. Sampling/detection of the leading ramp first becomes nonzero one dt after its theoretical onset; reports retain this difference instead of relabelling the sample.

## Verification
11 standard-library unittest checks cover:
- supplied/derived analytical values and invalid inputs;
- frictional steady state and pressure/flow boundary conditions;
- the independently calculated first valve and interior updates;
- full finite-ramp frictionless histories at all four stations using travelling/reflected waves;
- period and arrival times, including 10 m overlap for a=1100;
- explicit off-grid station sampling;
- independent sinusoid FFT frequency/amplitude;
- unknown, entered and crossed pressure-limit diagnostics;
- frictional 100/200/400 reach comparisons.

The independent wave comparison applies before and after reflections mathematically; it does not make an invalid single-phase trough physically real. On a matching characteristic grid the frictionless solution can agree to roundoff; this is expected, not evidence that all physical approximations are exact. In the default friction sensitivity, 200-to-400 reach peak change is about 0.0004695%, and maximum difference sampled at common times is about 0.0002853% of the ideal surge. These checks and the proposed <1% refinement target are not assessment tolerances or experimental validation.

| Selected result at N=400 | First valve rise (MPa) | Measured period (s) | FFT peak (Hz) |
|---|---:|---:|---:|
| elastic_ideal | 0.681067 | 0.400841821 | 2.494750 |
| bulk_ideal | 0.750750 | 0.363636364 | 2.750000 |
| elastic_friction_sensitivity | 0.682310 | 0.400841829 | 2.494750 |

The slightly higher frictional first maximum includes loss recovery; exact Joukowsky equality is checked using the separate frictionless cases. Code does not clip negative pressure or substitute the notebook's water vapour pressure.

## Files and provenance
| File | Role |
|---|---|
| stage1.py | Case data, speed, streaming MOC, pressure diagnostics, period and FFT functions |
| run_stage1.py | Three reference cases, grid studies, plots and data export |
| tests/test_stage1.py | Independent numerical/physical checks |
| requirements.txt | Tested dependencies |
| AGENTS.md | Persistent scope, physics and code-delivery rules |
| docs/PROJECT_CONTEXT.md | Verified coursework context, dependencies, decisions and unresolved questions |
| docs/SOURCES.md | Exact source pages/cells, original hashes, attribution and peer-reviewed extensions |

The supplied notebook is an attributed baseline, not an assertion that all defaults are coursework data. Detailed differences and known lecture discrepancies are retained in the context files. Original course PDFs, ZIP and full notebook are not redistributed; provide them separately for source-level review. No overall licence is assigned to the course-derived material by this repository.

## Continue in Codex
A separate application is not needed to develop this model here. For another Codex session, open this repository and checkout `stage1-pressure-waves`. Read AGENTS.md, PROJECT_CONTEXT.md, SOURCES.md and selected verification before editing. Repository guidance travels with the files; complete past chats and attachments do not automatically transfer.

Suggested handoff prompt:

> Work only in MinhAnh120/MECH0074 on stage1-pressure-waves. Read AGENTS.md, docs/PROJECT_CONTEXT.md and docs/SOURCES.md. Run the tests and inspect results/selected/verification.json. Explain the current physical assumptions and unresolved phase/friction/wave-speed choices, then implement my requested change. Return complete changed code in copyable fields with exact paths and commit.

The actual loop route, bend stations and structural response remain later work. Confirm dimensions/boundaries before replacing benchmark histories with bend histories. Discussion of applicable codes and design improvements is still required for the coursework and is not delivered by passing these numerical checks.
