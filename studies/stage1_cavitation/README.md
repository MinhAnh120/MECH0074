# MECH0074 · Separate Stage 1 cavity sensitivity

This exercise covers the 100 m Stage 1 benchmark only. Use the existing single-phase benchmark and the original exported histories for the group work. This package neither exports bend loads nor updates the shared code or presentation.

## Reproduce

Run the commands from the repository root using Python 3.12. The run used NumPy 2.3.5 and Matplotlib 3.10.8.

```sh
python -m pip install -r requirements.txt
python -m studies.stage1_cavitation.cavitation_sensitivity --out studies/stage1_cavitation/results
```

The new script calls the original single-phase solver for paired checks. The repository-root `stage1.py` and `fluid_assumptions.py` match the original Stage 1 branch at commit `be7f68ab1348566e677c8b58c22f372f017e7f39`. SHA-256 source fingerprints and input parameters are in `results/run_metadata.json`. Stage 1 course materials are referenced in the original repository and are excluded from this archive.

## How to read the outputs

- `cavitation_summary.csv` gives six meshes at the nominal methane proxy and selected meshes at a rounded pressure threshold and alternate cavity volume integration weight. “Later peak” is the maximum sampled valve pressure between first valve collapse and the end of the 1 s record. It is not necessarily the pulse emitted at first collapse.
- `valve_history_N640.csv` gives the plotted matched-grid time histories and the discrete cavity volume at the valve. A cavity volume is associated with that node. It is sensitive to the spatial representation and is not a measured vapour inventory.
- Discrete nodal cavities use a fixed vapour pressure and separate incoming and outgoing discharges. When integrated volume reaches zero, continuity is restored on the time grid. Collapse and reopening within a time step are unresolved. The optional counter in the CSV counts updates reset to a new cavity because the end-of-step liquid solution still falls below the proxy.
- The 0.100 MPa comparison is a rounded threshold, not a composition uncertainty bound. The main threshold is pure methane at 111 K and the benchmark still uses the prescribed LNG liquid density and bulk wave speed. The model does not include vapour transport, composition-dependent phase equilibrium, thermodynamic feedback or structural response. Values after cavity formation are exploratory.
- Validation checks cover the same original pressure histories before onset, a cavity-free slower-closure control, the analytical frictionless valve crossing and the cavity-volume balance during persistent cavity steps. There is no experimental validation of collapse pressure.

## References

- Setzmann and Wagner (1991), doi:10.1063/1.555898, Eq. (3.2) and Table 26.
- Simpson and Bergant (1994), doi:10.1061/(ASCE)0733-9429(1994)120:3(361).
- Bergant, Simpson and Tijsseling (2006), doi:10.1016/j.jfluidstructs.2005.08.008.
- Thol et al. (2019), doi:10.1063/1.5093800.
