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
| convergence.csv / convergence.png | First valve pressure rise in Pa and MPa for N=5,10,20,40,80,160,320; two speed panels |
| convergence_error.png | Nominal first-rise difference as a percentage of the reference rise and maximum waveform difference as a percentage of the Joukowsky rise. Both compare against N = 320. |
| theory_comparison.csv | Theory and measured first crest, arrivals, reservoir travel/return, period and FFT |
| friction_sensitivity.csv | Friction, surge, minima, phase flags, period and frequency |
| closure_sensitivity.csv | Additional0.020/0.10/0.25/0.50s closure scenarios; last case screened over50 periods |
| verification.json | Configuration, source hashes, code version, versions, input corners, per-grid metrics, phase flags and FFT details |

Generated results/latest is ignored by Git. results/selected stores current compact reports, ideal/nominal benchmark CSVs and selected plots. Initial results remain in Git history. Rerun to produce all larger histories and spectra.

An optional base-10 presentation is generated separately, without replacing any original convergence files:

```sh
python plot_logarithmic_convergence.py --input results/selected/convergence.csv
```

This adds `logarithmic_convergence.csv` (an identical data copy), `logarithmic_convergence.png` (base-10 logarithmic x, linear pressure y), and `logarithmic_convergence_error.png` (base-10 logarithmic x/y, retaining the1% criterion). The original base-2 plots remain available. Powers-of-ten ticks change the presentation, not the numerical data or relative point spacing. For a newly simulated run, use `--input results/latest/convergence.csv`; outputs default to the input directory. Re-run this plotting command whenever its input CSV changes. It does not execute the hydraulic solver.

The default study successively doubles N=5,10,20,40,80,160,320. Both runners accept `--reaches 5 10 20 40 80 160 320`; at least three successively doubled positive integers are required. All final histories, long records and closure scenarios use the finest grid, N=320 by default. At N=5, the 10, 25 and 50 m stations need spatial interpolation; at N=10, only 25 m does. All four stations lie on nodes for N=20–320. The CSV records this explicitly.

Long FFT records use50 ideal periods, starting one period after the closure ramp. Their analysed windows are18.181818s bulk and20.042091s elastic, with bin spacing0.055/0.049895Hz. M=64000 at N320; DFT duration=M*dt, not(M-1)*dt. Time-domain period and an independent sinusoid check supplement this coherently chosen FFT record. Peak recovery is numerical verification, not experimental evidence.

The chosen strict <1% criterion checks first-rise and period changes relative to the finer result, plus maximum absolute waveform difference divided by that case's Joukowsky rise. Coarse histories are linearly reconstructed on the finer time grid over the common first second at all four stations. This includes sampling/interpolation effects. Each candidate is compared with the next grid and with N=320, and must keep passing on subsequent tested grids. N=320 is a numerical reference, not an exact solution or a self-certified pass. The resulting smallest tested passing N is160 for bulk and20 for elastic in all five friction settings. Retain320 for a common production resolution. The peak alone is nearly unchanged even on coarse grids and cannot establish full-waveform convergence.

The peak-rise figure reports the first valve rise above its initial local pressure before the reservoir reflection returns. The error figure shows two of the three convergence metrics. The period check also appears in the CSV. The waveform metric compares each coarse trace with the finer trace at matching times and takes the largest absolute difference across all four stations.

## Interpretation
The figure label **Methane saturation pressure (111 K, proxy)** refers to 95,873.89895 Pa absolute. The actual LNG composition is unknown. All prescribed rapid-baseline cases cross this screen: full pressure traces, periods and spectra are **formal single-phase verification**. No pressure clipping, cavity model or physical approval is implied. The first positive crest occurs before the detected threshold crossing.

The model is a horizontal uniform-area, constant-property MOC with quasi-steady signed Darcy resistance, fixed reservoir and prescribed valve flow. It omits local bend losses/reflections, elevation, cavity dynamics and moving-pipe feedback. The viscosity/roughness envelope is an engineering scenario, not an actual-fluid confidence interval.

## Beau handoff
Read [BEAU_HANDOFF.md](docs/BEAU_HANDOFF.md). Benchmark stations10,25,50,100m are **not actual bends**. Supply the full reservoir-to-valve hydraulic length, four ordered bend stations, geometry revision and selected speed in a copy of route_template.json, then:

```sh
python export_bends.py route.json --out results/latest/bends --code-version YOUR_CODE_COMMIT
```

Null geometry is deliberately rejected. The exporter preserves a common time vector and each station’s wave delay and writes complete metadata/phase flags. Physical interpretation still needs actual composition and appropriate hydraulic modelling. Geometry-specific convergence must be checked after route changes.

## Verification and sources

The [Stage 1 evidence folder](results/selected/evidence/) contains plots and tables for the presentation and model appendix. Regenerate them after the full benchmark run:

```sh
python stage1_evidence.py --input results/latest --out results/selected/evidence
```

| Slide | Evidence to use |
|---|---|
| 13 | `benchmark_pressure_overlay.png` shows the four 1 s pressure histories with the analytical finite-closure solution. |
| 14 | [Theory comparison](results/selected/evidence/theory_comparison.png) compares the rise, return, period and FFT and isolates first-sample timing offsets. The [reflection zoom](results/selected/evidence/valve_reflection_zoom.png) shows when the reflected front starts and when the methane proxy threshold is crossed. |
| 15 | Use the [peak-rise plot](results/selected/logarithmic_convergence.png) with the separate [mesh-error plot](results/selected/logarithmic_convergence_error.png). The [standalone 3 s trace](results/selected/evidence/valve_pressure_3s.png) and [combined trace and spectrum](results/selected/evidence/valve_record_and_spectrum.png) are available for showing the longer valve response. |
| 32 | `closure_pressure_comparison.png` compares the recorded maxima and minima over the first 2 s. The cases are benchmark scenarios for discussion. |

The full time histories and spectra are generated in `results/latest` for the reproducible submission ZIP. See [executed validation](docs/VALIDATION.md) for the frequency estimator, the check at N = 640 and the limitations of the phase screen.

Twenty unittest cases cover the original eleven independent wave/steady/boundary/FFT/refinement checks, five researched-model/export tests, and four checks of between-sample waveform error, grid/station sampling, strict tolerance and invalid grid inputs. Synthetic geometry is used only in export tests. See [executed validation](docs/VALIDATION.md).

At Courant1 the ideal wave solution can match independent reflected finite ramps to roundoff. Grid convergence does not establish phase validity or code compliance. A proposed <1% refinement criterion is an engineering choice, not a UCL mark threshold. A slight first-crest increase with friction includes initial pressure-loss recovery.

Read AGENTS.md, [PROJECT_CONTEXT.md](docs/PROJECT_CONTEXT.md), [SOURCES.md](docs/SOURCES.md), [FRICTION_AND_PHASE.md](docs/FRICTION_AND_PHASE.md), and selected verification before continuing in another coding session. Earlier chats/attachments do not automatically transfer. The issued brief pp3–4 and supplied notebook control scope and teaching conventions; scientific extensions are attributed. Original PDFs, ZIP and notebook are not redistributed here. No licence over supplied course material is asserted.

Discussion of codes, actual structural response, damping sensitivity, and combined feasible designs remains group work. The longer0.50s pressure run is not the Stage4 preferred-design structural assessment.
