# Executed validation — 16 November 2026

Code revision: dd410277c465aeb68836565b45040b8024bff6ca. Reports store this revision plus seven independent source SHA256 values, including resolution.py. The subsequent results commit contains data/documentation only.

Executed successfully:

```sh
python -m unittest discover -s tests -v
python complete_stage1.py --out results/latest --code-version dd410277c465aeb68836565b45040b8024bff6ca
python run_stage1.py --out results/run-legacy-seven-grid --phase-limit-pa 95873.898949157
```

20 tests passed: the original11 numerical tests, five researched-model/export tests and four resolution tests. Python3.12.14, NumPy2.3.5, Matplotlib3.10.8. The new between-sample test independently demonstrates that comparing only at coarse times can miss waveform differences. Export tests use synthetic temporary geometry, not the actual loop. Null route input remains rejected.

The complete runner finished10 cases: two speeds times five friction settings. Each case used N=5,10,20,40,80,160,320 for the first second. All production histories, long records and closure scenarios use320. Each FFT window still spans50 ideal pressure periods, now64000 samples. The windows last18.181818s bulk and20.042091s elastic. All seven source hashes were checked against the executed files. The legacy runner also completed the seven-grid workflow.

| N=320 result | Bulk ideal | Elastic ideal | Bulk nominal | Elastic nominal |
|---|---:|---:|---:|---:|
| First valve rise Pa |750750.000000|681066.659443|751375.738385|681694.616157|
| Time-domain period s |0.363636364|0.400841821|0.363636365|0.400841823|
| FFT peak Hz |2.750000000|2.494749668|2.750000000|2.494749668|
| Minimum absolute Pa, first second |−50750.000000|18933.340557|−50098.301764|19584.932106|
| First detected global surrogate phase crossing s |0.200000000|0.219523529|0.200000000|0.219523529|

## Chosen numerical criterion and corrected waveform comparison
The strict <1% acceptance rule is a project engineering choice, not a coursework marking threshold. For each candidate grid, check against both the next grid and the finest numerical reference:

- Absolute first-rise difference divided by the finer-grid first rise.
- Absolute time-domain period difference divided by the finer-grid period.
- Maximum absolute pressure-history difference across all four physical recording stations over the common first second, divided by the same-speed analytical Joukowsky rise.

For the waveform check, linearly reconstruct the coarse history on the finer time coordinates. This measures the reproducibility of the reported histories, including temporal and station interpolation; it is not an isolated spatial truncation error. The previous comparison interpolated the finer result only at coarse times, which could miss between-sample differences. Its old waveform numbers must not be reused as the current convergence evidence.

A candidate is labelled the smallest tested passing grid only if it and every subsequently tested non-reference grid pass both comparisons for all three quantities. The finest320 grid is labelled a reference with blank self-comparison/acceptance fields, rather than passing by comparison with itself. This finite study does not establish an exact solution, asymptotic order or a universal minimum grid.

| Nominal case | Smallest tested passing N | Peak difference160→320, % | Maximum waveform difference160→320, % |
|---|---:|---:|---:|
| Bulk |160|0.0002690844|0.5678293005|
| Elastic |20|0.0002963966|0.2104142421|

The same smallest tested passing N values apply to all five friction settings. The elastic20 result is specific to these properties, closure law and output checks; repeat the study after changing them. N=320 is retained as a common production resolution for both speeds. The ideal peak is flat to roundoff even on coarse grids; frictional peak changes are very small. A flat peak alone does not establish full-waveform convergence.

At N=5 the10,25,50m stations need spatial interpolation; at N=10 only25m does. All stations fall on nodes for20–320. The CSV records this flag, first rise in both Pa/MPa, periods, analytical surge errors, next-grid and finest-reference comparisons. It contains70 rows, ten cases times seven grids. MPa values were checked against Pa/1e6.

The peak-rise PNG and separate error PNG were visually reviewed. The peak plot has separate bulk/elastic panels, MPa labels and explicit successive-doubling horizontal scale. Legends lie below the data. The diagnostic plot retains the1% line and excludes self-comparison at320. Nominal history/spectrum plots retain time labels, longer-record titles, phase thresholds and unclipped curves. Benchmark fronts are independently measured; fixed-reservoir arrival uses flow. Reflection return is measured in frictionless cases only, not fabricated for frictional cases.

## Physical and handoff limits
Physical validity has NOT passed: all0.020s baseline cases cross the sourced pure-methane surrogate threshold. First positive crests precede the crossing; complete later histories remain formal single-phase continuation. Actual LNG composition, mixture bubble point, pipe finish, cavity dynamics, actual bend-route convergence, experimental validation and code compliance remain unestablished. Closure scenarios are pressure studies only, not approved operating instructions or structural designs.

The selected benchmark history CSVs have been regenerated at320 and checked to be unchanged from the previous six-grid run; adding N=5 does not alter their time coordinates or values. The verification report records the new executed code revision. Beau must replace older400-grid files as a complete time/history set, check verification.json's code revision, and preserve the pressure convention/phase flags. Actual bend histories still require agreed route geometry and speed. Only stage1-pressure-waves was updated; no other branch was merged or modified.

## Separate base-10 convergence presentation — 16 November 2026
Plotting code revision:5355290aaae173beb262e6a8447eb93decd0cba8. Executed `python plot_logarithmic_convergence.py --input results/selected/convergence.csv` and `python -m unittest discover -s tests -v`; all20 existing tests passed. This is postprocessing of the stored70-row convergence table, not a new hydraulic run. The hydraulic code, selected verification.json and its recorded solver revision remain unchanged.

Created logarithmic_convergence.csv, logarithmic_convergence.png and logarithmic_convergence_error.png in results/selected. The new CSV is an identical copy of convergence.csv. The peak figure uses base-10 logarithmic x and linear pressure y; the diagnostic figure uses base-10 logarithmic x/y and preserves the1% criterion. Both figures were visually checked for powers-of-ten ticks, readable labels/legends, correct series and reference exclusion. All original convergence files and benchmark pressure/flow histories are retained unchanged.

Re-run the separate plotting command after changing the input convergence CSV. Source identification (CSV canonicalised to UTF-8 with LF newlines):

- Input CSV SHA256: `dde326281e1d88312c104a5590c7efa3a0f0e6bd5049b28e146eb264070013f1`
- Plotting script SHA256: `d1f6c2df11da40dab3c4267c81f726c6c6deb841b96bebd601b732d8f99097fb`

## Presentation evidence — 19 November 2026

The optional `stage1_evidence.py` runner at revision `5ffa8df50b7c1a5ea026a7789e00c85af59f6bce` reads the full researched benchmark outputs. The supplied benchmark inputs, MOC solver, production resolution and oscillation records are the same as in the seven-grid run above.

```sh
python -m unittest discover -s tests -v
python complete_stage1.py --out results/latest --code-version dd410277c465aeb68836565b45040b8024bff6ca
python stage1_evidence.py --input results/latest --out results/selected/evidence
```

The full rerun and evidence generator completed with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8. All 20 tests passed. The rerun reproduced the tracked verification report, convergence, friction, closure and theory tables, and both nominal pressure/flow CSVs. The evidence generator checks the reported source hashes before combining stored and new calculations.

| Check | Bulk | Elastic |
|---|---:|---:|
| Ideal 10 m first pressure rise, MPa | 0.682500 | 0.681067 |
| Nominal fundamental from the approximately 17 s window, interpolated spectrum, Hz | 2.7491 | 2.4954 |
| Nominal fundamental from repeated pressure crossings in that window, Hz | 2.7500 | 2.4947 |
| Nominal 320-to-640 maximum waveform difference, % of Joukowsky rise | 0.56776 | 0.21033 |

The frictionless MOC and analytical finite ramps agree at the recorded N = 320 samples to within 1e-8 Pa. Dotted arrival lines indicate the mathematical start of each ramp. The table reports the first detectable numerical sample, one time step later in these aligned cases. At 10 m the incident and reflected fronts are 18.18 ms apart for bulk speed, shorter than the 20 ms closure. At elastic speed the separation is 20.04 ms. The analytical overlay and its reflected-wave zoom show the frictionless verification case. The diagnostic at N = 640 is separate from the original seven-grid series. Its first-rise and period differences also meet the strict 1% benchmark criterion.

The plotted spectrum uses the original 50-period Hann-windowed, mean-removed record. The 0–3 s panel is a visual crop. The second FFT checks the same nominal simulations over 0.5–17.5 s with approximately 0.059 Hz native bin spacing. Its peak is refined with a parabola through three adjacent log-amplitude bins. The crossing estimate uses the median time between equivalent upward crossings of the reservoir pressure. Harmonic amplitudes are single-sided valve-pressure amplitudes at the bins nearest the first four odd harmonic targets. The FFT precision represents estimator agreement for this formal simulation. See Harris (1978) for finite records and spectral windows, as cited in SOURCES.md.

The closure figure compares each model over a common 2 s window, ending at the last full time step. The original longer 0.50 s scenario remains in the existing sensitivity table. The 0.02 s and 0.10 s nominal scenarios cross the pure-methane surrogate at 111 K. The 0.25 s and 0.50 s scenarios stay above that surrogate over the comparison window. Route-specific loading and an LNG-composition phase assessment are still needed to establish an operating closure time. Use the valve spectrum as a benchmark when interpreting the eventual local bend-force spectra.

## Figure updates — 25 November 2026

The peak-rise plot reports the first valve pressure rise before the reflected wave returns. The mesh-error panels define both denominators. First-rise difference is relative to the reference mesh's first rise. The maximum waveform difference is relative to the Joukowsky rise and covers matching times in the first second at all four stations. For the nominal bulk case, N = 160 differs from N = 320 by 0.00027% in peak rise and 0.568% in the waveform metric. These are numerical reference comparisons. The frictionless theory comparison is separate and shows how first-detection times differ from the theoretical front times by one time step.

Bulk figures use c for the supplied sound speed and elastic figures use a for the effective wave speed. The spectral markers show the first four predicted odd harmonics. The 0–1 MPa display includes every peak in the ten benchmark cases. The methane proxy label retains the temperature and phase basis. The figure notes and threshold-crossing shading identify the single-phase continuation.

All 20 tests passed. The full researched and legacy runners and both presentation plotting commands completed using Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8. The tracked CSV contents, benchmark case results and phase flags reproduced exactly. The updated verification report records the current source hashes and code revision.
