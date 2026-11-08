# Source provenance

## Supplied sources
1. `MECH0074_CW1_2026_27_v5.pdf`, pp.3-4: mandatory technical scope, inputs and justification of missing data; pp.1-2,6-7: submission, AI policy and rubric. SHA256 `caed2203a6651208bc1fc657e28417129cc8fab1e3f25ad149ea2ae5e8729eb5`.
2. `MECH0074_Slides_Lecture6_v7.pdf`, p.39: effective speed; pp.38-40: elastic surge; pp.41-44: propagation/reflections; pp.45-47: closure; p.48: protection examples. The full loop support/FSI model is not specified by these slides.
3. `MECH0074_Slides_Lecture5_v10.pdf`, pp.40-44: wave equation, mixed-boundary frequencies and pressure-flow equations; p.51: impedance.
4. `water_hammer_valve_closure_and_bend_forces2.ipynb`, 21 cells, zero-based indexing. SHA256 `46795cb3ba37fd26fff04aa31fa213861291908e2d904b71bae115fc6de05ddf`. Cells3,6-7: hand targets, head-flow MOC, initial gradient, boundary equations; cell9: steady/Joukowsky checks; cells10-11: plots; cells13,20: phase limits; cells15-18: compact bend-force postprocessing.
5. `IMECE2012-87585.pdf`, Flora et al. (2012), "Evaluation of Numerical Methods Applied in the Analysis of the Transient Flow of Pipeline Systems", pp.2-3 equations/initial conditions/MOC; pp.5-8 method comparisons. Supplied ASME conference paper.
6. `MECH0074_Source_Register_v2.md` and `MECH0074_Project_Context_Audit.md`: prior source inventory and discrepancies, not independent authority over the issued brief. Original files remain outside this repository.

## Attribution of implementation
Adapted from notebook cells6-7: characteristic coefficients, signed explicit quasi-steady friction, steady head gradient, fixed reservoir and prescribed-flow valve. The arithmetic is recognisably the supplied teaching formulation, not claimed as a newly invented solver.

New for this repository: parameter validation; lecture-derived speed and labelled alternative; streaming storage at requested stations; global pressure diagnostics; benchmark CSVs, FFT and period estimates; hydraulic refinement comparisons; source/decision records; tests against an independently superposed finite-ramp wave; portable handoff documentation. Plotting retains the teaching notebook's Matplotlib approach. Tests independently verify adaptation rather than trusting the notebook's baseline defaults.

The source notebook's full bytes and course PDFs/ZIP are not redistributed here. Place copies in an ignored local `references/` directory when using another environment.

## Peer-reviewed extensions
- Ghidaoui, Zhao, McInnis and Axworthy (2005), "A Review of Water Hammer Theory and Practice", Applied Mechanics Reviews 58,49-76. DOI 10.1115/1.1828050. Sections2-4: model restrictions, compliance and wall shear.
- Bergant, Simpson and Tijsseling (2006), "Water hammer with column separation: A historical review", Journal of Fluids and Structures 22,135-171. DOI 10.1016/j.jfluidstructs.2005.08.008. Phase change and collapse can invalidate single-phase continuation.
- Harris (1978), "On the use of windows for harmonic analysis with the discrete Fourier transform", Proceedings of the IEEE 66,51-83. DOI 10.1109/PROC.1978.10837. Uniform sampling, DFT bins and windowing. Hann window use and the chosen 50-period record are explicit analysis choices, not requirements stated by the coursework.
