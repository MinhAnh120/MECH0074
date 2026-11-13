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

## Research added 12 November 2026
- Sotiriadou, Antoniadis, Assael, Martinek & Huber (2026), International Journal of Thermophysics47:18, DOI10.1007/s10765-025-03690-7; online2025. Table8 p25 and liquid-uncertainty discussion pp1,27. https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=960665
- Colebrook (1939), J.Inst.CivilEngineers11(4),133–156, DOI10.1680/ijoti.1939.13150. Publisher metadata/preview accessible; full publisher text paywalled. Equation also cross-checked in published transient review/conference literature.
- Botros & Shouman (2016), J.NaturalGasScience&Engineering35,197–202, DOI10.1016/j.jngse.2016.08.044. Publisher indexed abstract/introduction accessible; full article direct retrieval403. Claims here are limited to the accessible stainless-pipe findings and correlation; no measured roughness interval is attributed to it.
- Setzmann & Wagner (1991), JPCRD20,1061–1155, DOI10.1063/1.555898. Eq3.2 p1082, Table26 p1086 (PDF pp23,27 including cover). https://srd.nist.gov/jpcrdreprint/1.555898.pdf
- Friend (1990), Cryogenics30,105–112, DOI10.1016/0011-2275(90)90254-A; accessible publisher abstract for composition dependence only.
- Thol, Richter, May, Lemmon & Span (2019), JPCRD48,033102, DOI10.1063/1.5093800. NIST publication page abstract accessed. Its journal metadata incorrectly says Journal of Chemical and Engineering Data; the DOI/paper identify JPCRD. It supports mixture-specific thermodynamics, not a numerical bubble-point result for an unknown mixture.
- Supplied params.py, 4233bytes, attachment9November: retained shared supplied/derived constants, corrected provenance wording only. New assumptions live in fluid_assumptions.py, not disguised as supplied parameters.
- MECH0074 Group Presentation Layout.pdf,36pages: pp11–15 Stage1 allocation; p15 requests peak pressure rise in MPa against reaches and the longer-record spectrum; p32 Minh-Anh recommendations; p35 model appendix. A planning document, not authority for prefilled numerical answers.

The13November resolution.py study, chosen N=10,20,40,80,160,320 sequence and strict1% criterion are project analysis choices. The maximum waveform error uses linear reconstruction of coarse output on finer times over the common first second, normalised by the same-speed Joukowsky rise. These choices are documented in VALIDATION.md; no coursework or literature source is claimed to prescribe this particular acceptance rule.

No source supplied an actual LNG assay or pipe roughness. The nominal roughness and uncertainty envelopes are explicitly engineering choices. See FRICTION_AND_PHASE.md for access limitations, source-specific claims and conditional validity.
