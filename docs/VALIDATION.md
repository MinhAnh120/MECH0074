# Executed validation — 13 November 2026

Code revision: c9a4b1b341e4350c43381d5d4f822787c7160327. Reports store this revision plus seven independent source SHA256 values, including resolution.py. The subsequent results commit contains data/documentation only.

Executed successfully:

```sh
python -m unittest discover -s tests -v
python complete_stage1.py --out results/latest --code-version c9a4b1b341e4350c43381d5d4f822787c7160327
python run_stage1.py --out results/run-legacy-six-grid --phase-limit-pa 95873.898949157
```

20 tests passed: the original11 numerical tests, five researched-model/export tests and four resolution tests. Python3.12.14, NumPy2.3.5, Matplotlib3.10.8. The new between-sample test independently demonstrates that comparing only at coarse times can miss waveform differences. Export tests use synthetic temporary geometry, not the actual loop. Null route input remains rejected.

The complete runner finished10 cases: two speeds times five friction settings. Each case used N=10,20,40,80,160,320 for the first second. All production histories, long records and closure scenarios use320. Each FFT window still spans50 ideal pressure periods, now64000 samples. The windows last18.181818s bulk and20.042091s elastic. All seven source hashes were checked against the executed files. The legacy runner also completed the six-grid workflow.

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

Only the25m station on N=10 needs spatial interpolation; all stations fall on nodes for20–320. The CSV records this flag, first rise in both Pa/MPa, periods, analytical surge errors, next-grid and finest-reference comparisons. It contains60 rows, ten cases times six grids. MPa values were checked against Pa/1e6.

The peak-rise PNG and separate error PNG were visually reviewed. The peak plot has separate bulk/elastic panels, MPa labels and explicit successive-doubling horizontal scale. Legends lie below the data. The diagnostic plot retains the1% line and excludes self-comparison at320. Nominal history/spectrum plots retain time labels, longer-record titles, phase thresholds and unclipped curves. Benchmark fronts are independently measured; fixed-reservoir arrival uses flow. Reflection return is measured in frictionless cases only, not fabricated for frictional cases.

## Physical and handoff limits
Physical validity has NOT passed: all0.020s baseline cases cross the sourced pure-methane surrogate threshold. First positive crests precede the crossing; complete later histories remain formal single-phase continuation. Actual LNG composition, mixture bubble point, pipe finish, cavity dynamics, actual bend-route convergence, experimental validation and code compliance remain unestablished. Closure scenarios are pressure studies only, not approved operating instructions or structural designs.

The selected benchmark history CSVs have been regenerated at320. Beau must replace older400-grid files as a complete time/history set, check verification.json's code revision, and preserve the pressure convention/phase flags. Actual bend histories still require agreed route geometry and speed. Only stage1-pressure-waves was updated; no other branch was merged or modified.
