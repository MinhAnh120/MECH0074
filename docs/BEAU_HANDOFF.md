# Stage 1 handoff to Beau

## Ready now
Use `results/selected/theory_comparison.csv`, `friction_sensitivity.csv`, `closure_sensitivity.csv` and `verification.json`. The selected bulk/elastic nominal benchmark history files have time_s, four absolute pressures and four signed volume flows. Their stations are **10,25,50,100 m benchmark locations, NOT the four actual bends**. Each row provides the common physical time for that case. Do not combine rows by index between speed cases: their time steps differ.

`complete_stage1.py` generates all ten cases, short histories, longer valve histories for ideal/nominal cases, spectra and plots. Every result has source hashes and a code-version field. The report identifies the code commit; the results commit is a later snapshot. Original brief/notebook material remains outside this public repository.

## Shared conventions
- SI: time s, hydraulic station m from reservoir along centreline, p_abs Pa, Q m3/s.
- External pressure = **101325 Pa absolute**. Structural pressure above the exterior is **p_abs - 101325**, not delta-p from the initial condition. Initial nominal reservoir gauge pressure = 598675 Pa.
- Positive flow runs reservoir to valve; reverse flow remains signed. Velocity=u=Q/A, A=0.495068446630 m2. rho=455 kg/m3; initial u=1.5 m/s, Q0=0.742602669945 m3/s.
- Reservoir fixed at 700000 Pa absolute; initial pressure decreases downstream by the steady Darcy loss. The brief does not specify the pressure station, so this is a documented interpretation.
- Valve flow: Q_valve=A U0 max(1-t/0.020,0). This prescribes velocity/flow closure, not valve-area closure.
- Explicit alternatives: bulk a=1100 m/s; simplified lecture elastic a=997.8998673 m/s. Bulk is selected by the supplied params flag, but group choice remains TBC; do not mix speeds in comparisons.
- Nominal constant Darcy f=0.01010783190; assumed-input sensitivity 0.00881783056–0.01163976282; f=0 and 0.020 are separate verification/teaching cases. See FRICTION_AND_PHASE.md for the qualified justification.
- Auxiliary phase screen: pure methane at 111 K, p_sat=95873.89895 Pa absolute. Actual LNG composition is unknown. All rapid-baseline cases cross this threshold and are formal single-phase verification.

## Geometry needed before actual bend histories can exist
Provide the total reservoir-to-valve hydraulic centreline length, four ordered bend stations, a geometry revision/source, and the selected speed mode. Four right-angle bends of radius 1.4 m contribute 4*(pi/2)*1.4=8.79646 m to total arc length, in addition to straight runs. Agree whether a representative station is the bend midpoint; a compact single-pressure bend is an approximation. A resolved control-volume calculation could require inlet/outlet pressures, flows and momentum storage instead.

Fill `route_template.json`; nulls are deliberate and cannot be run as real geometry. Specify duration/reaches appropriate to the *actual* route, and refine its hydraulic grid. The 100 m benchmark convergence does not establish convergence for every new route. No arbitrary four bend stations are invented.

```sh
python export_bends.py route.json --out results/latest/bends --code-version CODE_COMMIT_SHA
```

The exporter writes `bend_histories.csv`: one time_s column, p_abs_bend1..4_pa, Q_bend1..4_m3_s, and before_global_phase_crossing. Metadata contains route/source, complete case, exact dt/dx, interpolation flag, phase basis/crossing time, conventions and source hashes. An externally justified phase threshold can be entered as phase_limit_pa plus phase_basis in route.json; this does not implement cavitation.

The exporter preserves every station's wave delay through the shared MOC solution; do not apply a second manual delay. It assumes a uniform horizontal pipe with compact sampling stations, no local bend losses/reflections and no moving-pipe feedback. Changes to network boundaries, elevation or area require a revised hydraulic model, not merely changed station positions.

`before_global_phase_crossing=1` means only “before the first detected global threshold crossing,” not approved_physical_load=true. Invalid continuation after that crossing must be labelled. A small positive pressure does not by itself establish liquid-phase validity.

## Coupling boundary
Use each bend's own p(t), Q(t), geometry and directions; check the complete uniform-pressure force AND moment balance including relevant terminal tractions. Zero global resultant does not imply zero individual bend loads or zero local stress/deformation. First-crest agreement and numerical convergence do not certify the full transient as an LNG design load. Actual structural damping, modes, stress allowables, code assessment and feasible combined modifications belong to the group Stage4 model.
