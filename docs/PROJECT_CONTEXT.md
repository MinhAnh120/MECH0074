# Project context and Codex handoff

## Identity, scope and authority
MAN's separate UCL MECH0074 2026/27 group coursework concerns LNG water hammer in cryogenic expansion loops. This repository is `MinhAnh120/MECH0074`. The user owns Part B Stage 1: "What pressure rise and forcing frequencies follow valve closure?" Python and the supplied teaching notebook's MOC formulation are the agreed starting point. Full historical chats are not available here; this document transfers verified relevant context, not an assertion of complete chat access.

Authority order: issued CW1 v5 brief; latest supplied lectures; supplied technical references/notebook checked independently; clearly identified engineering choices and external peer-reviewed literature. See SOURCES.md for precise references.

## Assessed work
- Part A: cryogenic expansion-loop review, material/size/flow/contraction, comparison with steam, anchors/guides, design guidance.
- Stage 1: independent 100 m benchmark; pressure histories at 10, 25, 50 m from the reservoir and the valve; 1 s transient plot; pressure rise, travel time and 4L/a period against theory; longer frequency record; hydraulic resolution; industrial context and design improvements.
- Stage 2: fluid-filled cantilevers 4, 10, 20 m; first six bending frequencies against Euler-Bernoulli theory; mesh sensitivity.
- Stage 3: first six in-plane loop modes/frequencies; justify fluid mass and restraints; vary dimensions and guides.
- Stage 4: derive 90-degree bend force direction; map each bend's own hydraulic history with delay; uniform-pressure force AND moment balance including terminal loads; compare forcing and modes. Cryogenic baseline plus at least two modifications; compare thermal movement, surge displacement, support force and nominal stress; structural/hydraulic resolution, damping sensitivity and longer preferred run.
- Group deliverables: 2-page formative; UCL-template slides PDF with references/model appendix; >=20-minute face-visible MP4 with speaker identification and captions/transcript; reproducible code/configuration/output ZIP; individual IPAC. Dates remain on Moodle. Rubric: organisation 25%, technical content 50%, critical analysis 25%. Numerical AI category 3; presentation category 2: consult the actual brief for policy detail.

## Parameters from brief p. 4
Do=0.813 m; Di=0.79394 m; LNG rho=455 kg/m3; bulk sound speed=1100 m/s; initial velocity=1.5 m/s; initial absolute pressure=0.7 MPa; external absolute pressure=0.101325 MPa; prescribed linear velocity closure=0.020 s. Steel 304: rho=8000 kg/m3, E=213.233 GPa, thermal expansivity=1.4e-5 /K. Uniform cooldown 293 to 111 K. Four tangent bends, radius=1.4 m. Straight runs, full hydraulic length/stations and support details are unspecified.

Derived: wall=0.00953 m; A=0.49506844663 m2; Q0=0.74260266995 m3/s. Lecture 6 p.39 formula 1/a^2=1/c_bulk^2+rho*D/(wall*E) gives a=997.8998673 m/s. Ideal surge=0.681066659 MPa; L/a=0.100210455 s; 2L/a=0.200420911 s; period=0.400841821 s; f1=2.494749668 Hz. Fixed-area a=1100 reference gives 0.75075 MPa, 0.090909091 s, 0.181818182 s, 0.363636364 s, 2.75 Hz.

## Decisions versus unresolved questions
Implemented choices: uniform horizontal reservoir-pipe-valve benchmark; negligible convection; constant effective speed; quasi-steady Darcy friction; external pressure fixed; 0.7 MPa assigned to the reservoir boundary; old-to-new MOC updates. Cases distinguish elastic/fixed-area speed and frictionless/assumed f=0.020. These choices do not resolve the following:
- Does the lecturer intend 1100 directly in the benchmark or the lecture's pipe-compliance correction? Show both consistently; do not silently conflate speeds.
- The pressure location becomes important with friction: the brief does not explicitly locate the supplied 0.7 MPa in a frictional steady state.
- LNG friction, viscosity, roughness, composition and phase-pressure limit are not supplied. f=0.020 comes from the water notebook only.
- Ideal minima are 0.018933341 MPa absolute (elastic) and -0.05075 MPa (fixed-area). The first needs an LNG phase check; the second crosses zero. No result is labelled physically validated when the phase limit is unknown. No cavitation or structural feedback is implemented.
- Actual loop dimensions/stations, losses/reflections, boundary network, supports, damping and detailed cryogenic property/strength data remain to be agreed. Numerical tolerances in code are proposed engineering checks, not UCL thresholds.

## Member dependencies
User: Stage 1; group discussion names Wei for Stage 2, Kabir for Stage 3, Beau for Stage 4. The subsequently supplied 36-slide layout identifies Rali for Part A; Minh-Anh, Wei Yit, Kabir and Beau present Stages1–4 respectively. Stages 1 and 2 can begin independently; they share geometry/material/fluid data. Loop geometry/support definition can start in parallel. Final Stage 3 relies on structural verification. Stage 4 load-transfer derivation can start early, but final response needs local Stage 1 histories and a verified cooled structural model/modes. Do not describe Stage 2 as needing Stage 1 pressure results.

## Known source discrepancies
- Notebook hand exercise uses 0.25 s; solver baseline defaults to 0.01 s; assessed closure is 0.020 s.
- Notebook uses water, 600 m, 0.30 m, rho=998, a=1000 and water vapour pressure 2340 Pa. These are not LNG coursework inputs.
- Lecture 6 solution p.4 prints rho*U*deltaU and omits density from its compact wave-speed expression. Use brief/lecture 6 pp.39-40 and independently consistent dimensions, recording disagreement.
- Lecture 5 worked beam notes have mass/sign/root problems. Stage 2 needs positive steel-plus-LNG line mass and exact cantilever roots, not an odd-half-pi approximation for the first mode.
- The notebook bend calculation samples a straight-pipe solution; moving the bend sample does not generate hydraulic reflections.
- API674 supplied is historical reciprocating-liquid-pump guidance; API618 explanation concerns gas compressors. Neither establishes general LNG pipe compliance or a universal separation margin.

## LAMMPS compartment
Prior archive review identified C60, CNT, notched iron, hydrogen grain-boundary and illustrative Fe/Fe-C atomistic examples. These are not LNG hydraulics or 304 stainless cryogenic properties. The saved fifth example was incomplete; no LAMMPS rerun was verified. The archive and detailed review must be provided separately if future work needs them. They are not dependencies for Stage 1.

## Continue in Codex
Open/clone this repository, select `stage1-pressure-waves`, and read AGENTS.md, this file, SOURCES.md, README.md and results/selected/verification.json. Run the documented tests and benchmark before changing the model. Provide the original brief, relevant lectures and notebook separately if source-level reassessment is needed: repository summaries do not replace inaccessible attachments. Never assume that a new Codex session can see this chat.

Suggested first prompt: "Work only in MinhAnh120/MECH0074 on stage1-pressure-waves. Read AGENTS.md and docs/PROJECT_CONTEXT.md, inspect SOURCES.md and selected verification, run the tests, then explain the current model and unresolved physical assumptions before making the requested change. Return complete changed code in copyable fields."

## 12 November 2026 researched extension
Read FRICTION_AND_PHASE.md and BEAU_HANDOFF.md for current decisions. The original unresolved items above describe the initial implementation. A provisional nominal Darcy factor 0.01010783190 and assumed-input range 0.00881783056–0.01163976282 are now calculated using Colebrook. Methane surrogate phase pressure 95873.89895 Pa at111 K is now supplied for screening; actual LNG composition and pipe finish remain unknown. All prescribed rapid-baseline cases cross that screen, so neither elastic nor bulk full histories are physically approved. complete_stage1.py reads the supplied params.py and runs both speed interpretations. No actual bend histories exist until geometry is supplied; export_bends.py enforces this. See current selected results and executed validation.

## 13 November 2026 resolution update
Both runners now share resolution.py and the requested N=10,20,40,80,160,320 sequence. Production histories and sensitivities use320; the route template uses320 only as an initial value pending route-specific checks. Slide15 of the supplied presentation layout calls for peak pressure rise in MPa against reaches. The two speed interpretations remain separate panels, with a supplementary waveform-error plot and explicit1% diagnostics in CSV/JSON. The tolerance and comparison method are engineering choices, not prescribed UCL assessment thresholds.

The old waveform comparison sampled the fine result only at coarse times and could miss differences between samples. Current verification reconstructs the coarse output on finer times and reports next-grid and finest-reference errors across all four benchmark stations over the first second, normalised by the same-speed Joukowsky rise. First rise and period are checked separately. The smallest tested passing N is160 for bulk and20 for elastic; this does not establish a universal minimum or phase validity. Earlier N=400 results remain in Git history and must not be mixed with the new320 time coordinates. Actual geometry and mixture composition remain outstanding.
