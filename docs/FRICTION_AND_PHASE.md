# Friction and phase-pressure decision record — 12 November 2026

## Status and authority
The coursework brief `MECH0074_CW1_2026_27_v5.pdf`, pp.3–4, supplies density, bulk sound speed, geometry, pressure, temperature and closure law. It supplies neither viscosity, internal surface roughness, nor LNG composition. The original notebook cells6–7 use quasi-steady Darcy friction; its water factor 0.020 and water vapour pressure 2340 Pa are not LNG measurements. None of the new assumptions below changes the issued density 455 kg/m3 or bulk sound speed 1100 m/s.

These are reproducible **provisional engineering assumptions**, not newly discovered coursework inputs or a characterised industrial LNG. All numerical comparisons keep bulk and lecture elastic wave speeds distinct. The supplied params.py is now included; complete_stage1.py reads its shared values. Its USE_KORTEWEG=False is still a TBC group choice, so both modes are run explicitly.

## Viscosity: what the evidence supports
Sotiriadou et al., International Journal of Thermophysics (2026) 47:18, published online 8 December 2025, DOI [10.1007/s10765-025-03690-7](https://doi.org/10.1007/s10765-025-03690-7), Table8 p.25 and §3.1: saturated-liquid methane viscosity is 121.15 microPa s at 110 K and 98.180 at 120 K. Linear interpolation to 111 K gives 118.853 microPa s. Adopt a rounded **120 microPa s = 0.000120 Pa s** as a methane-like screening proxy. This interpolation is not an exact evaluation at 111 K, 0.7 MPa. The paper's estimated 3% liquid correlation uncertainty does not quantify unknown LNG composition, pressure correction or interpolation error.

An assumed **±20% model-input envelope**, 96–144 microPa s, is used for sensitivity. It is broader than the correlation uncertainty and deliberately allows proxy error. It is an engineering scenario, not a measured composition envelope or a 95% confidence interval. Without composition data, no finite range can be certified to encompass the actual LNG. Friend (1990), Cryogenics 30,105–112, DOI 10.1016/0011-2275(90)90254-A, independently demonstrates composition-dependent methane–ethane mixture viscosity; it does not identify this coursework mixture.

## Roughness: avoid false material certainty
Botros and Shouman (2016), Journal of Natural Gas Science and Engineering 35,197–202, DOI [10.1016/j.jngse.2016.08.044](https://doi.org/10.1016/j.jngse.2016.08.044), tested seven commercial stainless pipes. Their equivalent sand-grain roughness depends on surface topography, and differs from a simple arithmetic roughness Ra. Their stainless relation is ks=2.2907 Rq+0.1029 Rq² in micrometres. Therefore “304 stainless” alone cannot determine Colebrook roughness.

Use **epsilon=15 micrometres** as an explicit low-roughness nominal scenario and **0–50 micrometres** as an exploratory envelope from hydraulically smooth to a rougher metal surface. These numbers are **not quoted measured limits from Botros**, and no pipe-specific roughness evidence has been supplied. The nominal value corresponds to Rq≈5.3 micrometres if their empirical conversion is used; transferring that gas-flow relation to this diameter/finish/cryogenic service is itself approximate. Obtain an internal finish specification or pressure-loss calibration before treating it as established. No supplier chart or non-peer-reviewed web calculator was silently used as evidence.

## Calculation and use
Re=rho U0 D/mu=4,515,533.75 for the adopted nominal inputs. Solve Colebrook:

`1/sqrt(f_D) = -2 log10[epsilon/(3.7 D) + 2.51/(Re sqrt(f_D))]`.

Colebrook (1939), Journal of the Institution of Civil Engineers 11(4),133–156, DOI 10.1680/ijoti.1939.13150, is the primary correlation source. The implementation uses bisection, checks the residual and the independently tabulated smooth-pipe result at Re=100000. Use the **Darcy** convention, f_D=4 f_Fanning.

Nominal f_D=**0.01010783190**. Evaluate all four viscosity/roughness corners: the resulting envelope is **0.00881783056–0.01163976282**. This propagates the assumed input envelope; it does not establish the actual pipe's uncertainty. Keep f_D=0 for exact wave verification and 0.020 as a separate teaching comparison. At the nominal Re, f_D=0.020 implies equivalent roughness approximately 0.844 mm, far above the adopted 0–0.050 mm range; thus the notebook value is not selected as the nominal smooth-pipe estimate.

The nominal initial 100 m Darcy pressure loss is about 652 Pa. This explains why friction causes a small change to the first surge on the MPa scale. Friction can slightly increase the first *increment from the initial local pressure* through recovery of initial pressure loss; do not assert that every first crest must decrease. Use a consistent steady initial pressure gradient, signed Q|Q| resistance and a constant f_D selected at initial Re. This last approximation is not a validated treatment of flow reversal or near-zero-Re flow.

Ghidaoui et al. (2005), Applied Mechanics Reviews 58,49–76, DOI 10.1115/1.1828050, §4.1 pp.52–53, distinguish quasi-steady wall shear from transient contributions. The supplied quasi-steady model is retained for tractable teaching verification; its long-time damping is not experimentally validated. Structural damping is a different Stage4 parameter.

## Phase pressure: a conditional answer, not a universal LNG constant
Adopt **100 mol% methane only as an auxiliary screening surrogate**. Setzmann and Wagner (1991), Journal of Physical and Chemical Reference Data 20,1061–1155, DOI [10.1063/1.555898](https://doi.org/10.1063/1.555898), Eq.(3.2) p.1082 and Table26 p.1086, provide:

`theta=1-T/190.564; ln(p_sat/4.5992e6)=(190.564/T)*[-6.036219 theta+1.409353 theta^1.5-0.4945199 theta^2-1.443048 theta^4.5]`.

At **111 K**, this gives **95,873.89895 Pa absolute = 0.0958739 MPa absolute**. Tests independently check the published triple point and normal boiling point. The formula is not an Antoine fit outside its stated range.

This benchmark is intentionally hybrid: its supplied rho=455 and sound speed=1100 remain unchanged, while methane is used only for auxiliary viscosity/phase screening. It is not a thermodynamically self-consistent model of pure methane. Saturated methane density near 111 K is about 423 kg/m3, so the brief's density must not be cited as proof of pure composition.

An LNG mixture requires the **bubble-point pressure at the specified temperature and liquid composition**, using a suitable mixture equilibrium model; methane vapour pressure is not a universal conservative bound. Thol et al. (2019), *Journal of Physical and Chemical Reference Data* 48,033102, DOI 10.1063/1.5093800, develop EOS-LNG for multicomponent cryogenic thermodynamics. Exact composition-specific bubble pressure remains unresolved until composition is supplied and the selected model is assessed for that state.

## Results and permitted interpretation
The N=400 ideal bulk case reaches -50,750 Pa absolute and the ideal elastic case 18,933.34 Pa. Both are below the methane surrogate threshold. The nominal-friction cases also cross it. Their first global detected crossings are about 0.2000 s (bulk) and 0.21946 s (elastic). Thus use **formal single-phase verification**, including the recovered mathematical frequencies; do not transfer the complete post-crossing traces as validated physical LNG loads.

The first positive valve crest occurs before the first phase crossing in these calculations. Analytical agreement still verifies the numerical wave model, but not the later physical waveform. Bergant, Simpson and Tijsseling (2006), Journal of Fluids and Structures 22,135–171, DOI 10.1016/j.jfluidstructs.2005.08.008, pp.135–137: vapour cavities and their collapse can alter subsequent pressures. No pressure clipping or cavity dynamics has been introduced.

Resolving the old “unknown phase limit” label means replacing it with an explicitly sourced *surrogate screen*, not declaring physical validity. A cavitation model, or justified operating changes that avoid the relevant phase boundary, is needed for a physical continuation. Changing initial pressure/velocity/closure is permitted only as a separately labelled modification; the assessed 0.7 MPa/1.5 m/s/0.020 s baseline is retained.

## Remaining evidence needed
Actual LNG mole fractions (including nitrogen and heavier components), operating temperature specification, and pipe internal finish/roughness or measured pressure-loss data. The group must also settle the wave-speed interpretation. Formal benchmark verification can proceed without these; industrial physical validation cannot.
