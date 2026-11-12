# MECH0074: LNG pressure waves

MAN / Minh-Anh’s UCL MECH0074 coursework, 2026/27. Use only MinhAnh120/MECH0074, branch stage1-pressure-waves. No merge has been requested.

## Contents
- [Run](#run)
- [Cases](#cases)
- [Outputs](#outputs)
- [Interpretation](#interpretation)
- [Beau handoff](#beau-handoff)
- [Verification and sources](#verification-and-sources)

## Run
From a terminal with Python and the pinned requirements:

```sh
git clone --branch stage1-pressure-waves https://github.com/MinhAnh120/MECH0074.git
cd MECH0074
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python complete_stage1.py --out results/latest --code-version YOUR_CODE_COMMIT
```

Activation shown for Linux/macOS; Codespaces or another Python environment also works. Replace YOUR_CODE_COMMIT with the revision of the code being executed. Source SHA256 values are recorded independently, including when the working tree differs. These files are scripts: viewing them on GitHub does not execute them.

The original three-case teaching runner remains available:

```sh
python run_stage1.py --out results/run-legacy --phase-limit-pa 95873.898949157
```

Without --phase-limit-pa that legacy runner intentionally retains its old unknown-phase behaviour. Use complete_stage1.py for the researched workflow. It reads the supplied/shared params.py. The legacy numerical core remains independently testable and does not silently acquire new physical defaults.

## Cases
The required baseline remains100m, bore0.79394m, OD0.813m, rho455kg/m3, U0=1.5m/s, reservoir0.7MPa absolute, external0.101325MPa absolute, linear velocity closure0.020s. Initial frictional pressure includes the steady Darcy gradient. Assigning0.7MPa to the reservoir is a documented interpretation of an unspecified pressure station.

Both speeds are explicit: supplied bulk1100m/s and Lecture6 p39 elastic997.8998673m/s. The supplied USE_KORTEWEG=False flag is still a TBC group selection; no final group decision is invented. Both are run consistently with their own analytical targets.

For each speed, run ideal f=0, nominal Colebrook f=0.01010783190, assumed-envelope limits0.00881783056/0.01163976282, and separate teaching f=0.020. The provisional inputs are methane-like viscosity120microPa.s and equivalent roughness15micrometres. They are not measured properties of the actual mixture/pipe. Read [the research and qualifications](docs/FRICTION_AND_PHASE.md).

## Outputs
The complete runner generates:

| File | Meaning |
|---|---|
| *_histories.csv | Common time coordinates, four benchmark absolute pressures and signed volume flows,0–1s |
| *_spectrum.csv | Mean-removed, Hann-windowed valve pressure amplitude versus frequency |
| *_long_valve.csv | Full long valve histories for ideal/nominal cases |
| *.png | Required short histories plus longer-record spectrum, phase threshold and labels |
| convergence.csv / convergence.png | First rise and period for N=100,200,400 |
| theory_comparison.csv | Theory and measured first crest, arrivals, reservoir travel/return, period and FFT |
| friction_sensitivity.csv | Friction, surge, minima, phase flags, period and frequency |
| closure_sensitivity.csv | Additional0.020/0.10/0.25/0.50s closure scenarios; last case screened over50 periods |
| verification.json | Configuration, source hashes, code version, versions, input corners, per-grid metrics, phase flags and FFT details |

Generated results/latest is ignored by Git. results/selected stores current compact reports, ideal/nominal benchmark CSVs and selected plots. Initial results remain in Git history. Rerun to produce all larger histories and spectra.

Long FFT records use50 ideal periods, starting one period after the closure ramp. Their analysed windows are18.181818s bulk and20.042091s elastic, with bin spacing0.055/0.049895Hz. M=80000 at N400; DFT duration=M*dt, not(M-1)*dt. Time-domain period and an independent sinusoid check supplement this coherently chosen FFT record. Peak recovery is numerical verification, not experimental evidence.

## Interpretation
The auxiliary methane saturation threshold at111K is95873.89895Pa absolute. Actual LNG composition is unknown. All prescribed rapid-baseline cases cross this screen: full pressure traces, periods and spectra are **formal single-phase verification**. No pressure clipping, cavity model or physical approval is implied. The first positive crest occurs before the detected threshold crossing.

The model is a horizontal uniform-area, constant-property MOC with quasi-steady signed Darcy resistance, fixed reservoir and prescribed valve flow. It omits local bend losses/reflections, elevation, cavity dynamics and moving-pipe feedback. The viscosity/roughness envelope is an engineering scenario, not an actual-fluid confidence interval.

## Beau handoff
Read [BEAU_HANDOFF.md](docs/BEAU_HANDOFF.md). Benchmark stations10,25,50,100m are **not actual bends**. Supply the full reservoir-to-valve hydraulic length, four ordered bend stations, geometry revision and selected speed in a copy of route_template.json, then:

```sh
python export_bends.py route.json --out results/latest/bends --code-version YOUR_CODE_COMMIT
```

Null geometry is deliberately rejected. The exporter preserves a common time vector and each station’s wave delay and writes complete metadata/phase flags. Physical interpretation still needs actual composition and appropriate hydraulic modelling. Geometry-specific convergence must be checked after route changes.

## Verification and sources
Sixteen unittest cases cover the original eleven independent wave/steady/boundary/FFT/refinement checks plus correlation fixed points/residuals, shared input mapping, measured boundary timings, and the export contract using synthetic test geometry only. See [executed validation](docs/VALIDATION.md).

At Courant1 the ideal wave solution can match independent reflected finite ramps to roundoff. Grid convergence does not establish phase validity or code compliance. A proposed <1% refinement criterion is an engineering choice, not a UCL mark threshold. A slight first-crest increase with friction includes initial pressure-loss recovery.

Read AGENTS.md, [PROJECT_CONTEXT.md](docs/PROJECT_CONTEXT.md), [SOURCES.md](docs/SOURCES.md), [FRICTION_AND_PHASE.md](docs/FRICTION_AND_PHASE.md), and selected verification before continuing in another coding session. Earlier chats/attachments do not automatically transfer. The issued brief pp3–4 and supplied notebook control scope and teaching conventions; scientific extensions are attributed. Original PDFs, ZIP and notebook are not redistributed here. No licence over supplied course material is asserted.

Discussion of codes, actual structural response, damping sensitivity, and combined feasible designs remains group work. The longer0.50s pressure run is not the Stage4 preferred-design structural assessment.
