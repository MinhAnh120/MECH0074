# Executed validation — 12 November 2026

Code revision:166e671f3b2f9f66e47fb68b33469401683e0b5f. Reports store this code revision plus six independent source SHA256 values. The following results commit adds data/documentation only.

Executed successfully:

```sh
python -m unittest discover -s tests -v
python complete_stage1.py --out results/latest --code-version 166e671f3b2f9f66e47fb68b33469401683e0b5f
python run_stage1.py --out results/run-legacy --phase-limit-pa 95873.898949157
```

16 tests passed. Python3.12.14, NumPy2.3.5, Matplotlib3.10.8. Both the original11 numerical tests and five new researched-model/export tests pass. Exact sampled pressures/flows are checked against the exporter using synthetic geometry confined to temporary test files; this is not actual loop verification. Null route input is rejected.

The full new runner completed10 cases: two speeds times five friction settings, with N100/200/400,50-period FFT windows and four closure-time scenarios per speed. The0.50s closure case additionally ran for50 pressure periods. All6 report source hashes were checked against the executed files. The legacy three-case runner also completed with the methane surrogate phase threshold.

| N400 result | Bulk ideal | Elastic ideal | Bulk nominal | Elastic nominal |
|---|---:|---:|---:|---:|
| First rise Pa |750750|681066.659443|751376.143968|681695.021538|
| Time-domain period s |0.363636364|0.400841821|0.363636365|0.400841823|
| FFT peak Hz |2.75|2.494749668|2.75|2.494749668|
| Minimum absolute Pa |−50750|18933.340557|−50098.506491|19584.717854|
| First detected global surrogate phase crossing s |0.200000|0.219460897|0.200000|0.219460897|

The nominal200→400 reach first-rise changes are0.000215596% bulk and0.000237492% elastic. Maximum sampled waveform differences are0.000135375%/0.000144402% of Joukowsky surge. The proposed1% numerical criterion is satisfied. Frictionless agreement is to floating-point roundoff. All five independently detected station/boundary arrival fronts agree within one dt; reservoir arrival uses flow because its pressure is fixed. Return-front timing is independently measured in the frictionless cases; nominal friction cases do not fabricate a separate return detector.

Both nominal plots and convergence plot were visually reviewed: readable legends outside, top Time[s] label, longer-record spectrum title, phase-threshold line, no clipping of pressure curves. The convergence curves are nearly flat, as confirmed by the numerical differences.

Physical validity has NOT passed: all0.020s baseline cases cross the sourced methane surrogate threshold. The first positive crest precedes the crossing; the complete later history is only a formal single-phase continuation. Actual LNG composition/roughness, mixture bubble-point calculation, physical cavity dynamics, actual bend-route convergence, experimental validation and code compliance are not established. The0.25/0.50s modifications remain above the surrogate threshold over their reported durations only; they are not approved operating instructions or structural designs.
