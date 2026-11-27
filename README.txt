MECH0074 coursework - Stage 2 (cantilever verification)

Natural bending frequencies of the LNG-filled pipe as a cantilever (L = 4, 10, 20 m),
from Euler-Bernoulli theory, a beam FEM, and ANSYS, with mesh studies and a comparison.
Everything is in one script: stage2_all.py (ANSYS results are typed into it, no
spreadsheet needed).


1. One-off setup (Windows PowerShell, in this folder):

   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   python -m pip install -r requirements.txt

   (If activation is blocked: Set-ExecutionPolicy -Scope CurrentUser RemoteSigned)


2. Run everything:

   python stage2_all.py              saves all tables/figures, then opens the figure windows
   python stage2_all.py --no-show    saves all tables/figures only, no windows

   Without activating the venv:  .\.venv\Scripts\python.exe stage2_all.py --no-show

   Or open stage2_all.py in VS Code (with the .venv interpreter selected) and press Run.
   Outputs always go into this folder, whatever the working directory.
   At the end the script lists every expected output file as [OK] or [MISSING].
   If a figure fails to save with "Invalid argument", close it in any viewer/preview
   tab (the file is locked) and run again.


3. Output files (22, all written to this folder)

   Euler-Bernoulli theory and beam FEM
     stage2_EB_mode_shapes.png             mode shapes 1-6 vs x/L
     stage2_table_frequencies.csv          report table: theory + FEM, modes 1-6, each length
     stage2_table_comparison.csv           theory vs FEM (64 elements) with % difference
     stage2_table_mesh.csv                 full mesh study: frequency, % error, % change per doubling
     stage2_fig_mesh_EB.png                FEM error vs number of elements (all points)
     stage2_fig_refinement_EB.png          FEM % change per mesh doubling (all points)
     stage2_fig_mesh_EB_zeroed.png         as above, round-off points set to 0 (hollow markers)
     stage2_fig_refinement_EB_zeroed.png   as above, round-off points set to 0 (hollow markers)

   ANSYS mesh convergence (one set per length: 4m, 10m, 20m)
     ansys_<L>_frequencies.png             bending modes 1-6 vs number of elements
     ansys_<L>_change_Hz.png               change per mesh doubling, Hz
     ansys_<L>_change_pct.png              change per mesh doubling, %, with 0.1 % line
     ansys_bending_modes.csv               bending frequencies (duplicates removed), all meshes
     ansys_mesh_changes.csv                change per doubling (Hz, %), measured / upper bound

   ANSYS vs Euler-Bernoulli theory vs FEM
     stage2_table_ansys_vs_EB_FEM.csv      all three methods + % differences, each length/mode
     stage2_fig_ansys_vs_EB_FEM_freq.png   frequency vs mode, one panel per length
     stage2_fig_ansys_vs_EB_FEM_pct.png    % difference from Euler-Bernoulli theory vs mode


Notes for figure captions
   - Mesh-study errors are identical for 4, 10 and 20 m, so those figures plot one length.
   - Zeroed figures: hollow markers = values below floating-point round-off (~1e-6 %), set to 0.
   - Modes 5-6 have no result at 2 elements (too few degrees of freedom).
