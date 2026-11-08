# Executed verification - 8 November 2026

Environment: Python 3.12.14; NumPy 2.3.5; Matplotlib 3.10.8.

Commands executed successfully:

```sh
python3 -m unittest discover -s tests -v
python3 run_stage1.py --out results/latest
```

11 unittest cases passed. All three reference cases completed N=100/200/400 hydraulic studies and 50-period FFT windows. Source-code hashes in selected verification match the files used to generate the report.

Additional inspection of the generated report confirmed: ideal first surges agree with theory; period errors are below 1%; detected arrivals are within one time step of theoretical onset; successive peak, period and common-time waveform differences are below the proposed 1% criterion; FFT peaks are within one frequency bin of the analytical fundamental. Figures were visually reviewed for all three cases; the final legend/title layout was checked after the formatting revision.

The independently superposed finite-ramp wave test checks the complete frictionless pressure history at all four requested stations, including reflection signs and delay. Exact characteristic-grid agreement to roundoff is expected for that ideal case. The physical model still excludes phase change, transient wall-shear refinements and pipe-motion feedback.

All results are formal single-phase teaching/verification calculations. Bulk-speed minimum is -50750 Pa absolute; the elastic ideal minimum is 18933.34 Pa absolute. The LNG phase limit remains unknown, so no physical safety, code compliance, final design or actual-loop load validation is claimed.

The selected numerical report is in results/selected/. The command regenerates the visually reviewed figures in results/latest/. Reproduction commands and pending coursework tasks are in README.md and PROJECT_CONTEXT.md.
