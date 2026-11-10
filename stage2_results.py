"""
Stage 2 deliverable: "Compare 4, 10 and 20 m cantilevers with analytical
frequencies; check mesh sensitivity."

Uses the beam FEM in stage2_cantilever.py (keep in same folder):
  - EB FEM  : Euler-Bernoulli Hermite elements, consistent mass
  - TB FEM  : Timoshenko elements (Friedman & Kosmatka 1993), shear + rotary inertia
Outputs:
  stage2_table_comparison.csv   FEM vs analytical, 3 lengths x 6 modes
  stage2_table_mesh.csv         mesh study, both element types, all lengths
  stage2_fig_mesh.png           mesh convergence (log-log)
  stage2_fig_comparison.png     FEM / analytical ratio per mode and length
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stage2_cantilever import f_EB, f_FEM, m, EI, I, As, Do

LENGTHS = [4.0, 10.0, 20.0]
MESHES = [2, 4, 8, 16, 32, 64, 128]
NREF = 256                     # reference mesh for Timoshenko convergence
NCONV = 64                     # mesh reported as "converged"
r = np.sqrt(I / As)

# ---------------- mesh study ----------------
mesh_rows = []
mesh_err_EB, mesh_err_TB = {}, {}
for L in LENGTHS:
    fa = f_EB(L)
    f_ref_TB = f_FEM(L, NREF, timoshenko=True)
    mesh_err_EB[L], mesh_err_TB[L] = [], []
    for ne in MESHES:
        fe = f_FEM(L, ne)
        ft = f_FEM(L, ne, timoshenko=True)
        eEB = (fe - fa) / fa * 100            # vs exact analytical
        eTB = (ft - f_ref_TB) / f_ref_TB * 100  # vs finest Timoshenko mesh
        mesh_err_EB[L].append(eEB); mesh_err_TB[L].append(eTB)
        for n in range(6):
            mesh_rows.append([L, ne, n + 1, fe[n], eEB[n], ft[n], eTB[n]])

with open("stage2_table_mesh.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["L_m", "elements", "mode", "f_EB_FEM_Hz", "err_vs_analytical_pct",
                "f_TB_FEM_Hz", "change_vs_ne256_pct"])
    for row in mesh_rows:
        w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

# ---------------- comparison table ----------------
comp_rows = []
print("=" * 86)
print("COMPARISON: analytical Euler-Bernoulli vs FEM (ne = %d)" % NCONV)
print("=" * 86)
for L in LENGTHS:
    fa = f_EB(L)
    fe = f_FEM(L, NCONV)
    ft = f_FEM(L, NCONV, timoshenko=True)
    print(f"\nL = {L:.0f} m   (L/D = {L/Do:.1f}, slenderness L/r = {L/r:.1f}, "
          f"total mass = {m*L:,.0f} kg)")
    print(f"{'mode':>4} {'EB analytical':>14} {'EB FEM':>10} {'diff %':>8} "
          f"{'Timoshenko FEM':>15} {'diff %':>8}")
    for n in range(6):
        d1 = (fe[n] - fa[n]) / fa[n] * 100
        d2 = (ft[n] - fa[n]) / fa[n] * 100
        print(f"{n+1:>4} {fa[n]:>14.3f} {fe[n]:>10.3f} {d1:>8.3f} {ft[n]:>15.3f} {d2:>8.2f}")
        comp_rows.append([L, n + 1, fa[n], fe[n], d1, ft[n], d2])

with open("stage2_table_comparison.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["L_m", "mode", "f_EB_analytical_Hz", "f_EB_FEM_Hz", "EB_FEM_diff_pct",
                "f_Timoshenko_FEM_Hz", "Timoshenko_vs_EB_pct"])
    for row in comp_rows:
        w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

# ---------------- mesh table print (compact) ----------------
print("\n" + "=" * 86)
print("MESH SENSITIVITY (EB FEM, % error vs analytical; identical for all L)")
print("=" * 86)
print(f"{'ne':>5}" + "".join(f"{'mode '+str(n):>11}" for n in range(1, 7)))
for i, ne in enumerate(MESHES):
    print(f"{ne:>5}" + "".join(
        ("        n/a" if np.isnan(e) else f"{e:>11.4f}") for e in mesh_err_EB[10.0][i]))
print("\nMESH SENSITIVITY (Timoshenko FEM, %% change vs ne = %d)" % NREF)
for L in LENGTHS:
    print(f" L = {L:.0f} m")
    for i, ne in enumerate(MESHES):
        print(f"{ne:>5}" + "".join(
            ("        n/a" if np.isnan(e) else f"{e:>11.4f}") for e in mesh_err_TB[L][i]))

# convergence rate (EB, mode 6, between 16 and 32 elements)
e16, e32 = mesh_err_EB[10.0][3][5], mesh_err_EB[10.0][4][5]
print(f"\nObserved convergence order (mode 6, ne 16->32): {np.log2(e16/e32):.2f}  (theory: 4)")

# ---------------- figures ----------------
COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MUTED = "#8a8984"
def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
ax = axs[0]
for n in range(6):
    e = np.array([mesh_err_EB[10.0][i][n] for i in range(len(MESHES))])
    ok = ~np.isnan(e) & (e > 1e-5)   # below ~1e-5 % is floating-point round-off
    ax.loglog(np.array(MESHES)[ok], e[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"mode {n+1}")
ax.axhline(0.1, color=MUTED, ls="--", lw=0.8); ax.text(2.1, 0.12, "0.1 %", color=MUTED, fontsize=8)
style(ax)
ax.set_xlabel("number of elements"); ax.set_ylabel("error vs analytical (%)")
ax.set_title("(a) Euler–Bernoulli FEM vs exact solution\n(same for 4, 10, 20 m)", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8, ncol=2)

ax = axs[1]
for n in range(6):
    e = np.array([abs(mesh_err_TB[4.0][i][n]) for i in range(len(MESHES))])
    ok = ~np.isnan(e) & (e > 1e-9)
    ax.loglog(np.array(MESHES)[ok], e[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"mode {n+1}")
ax.axhline(0.1, color=MUTED, ls="--", lw=0.8)
style(ax)
ax.set_xlabel("number of elements"); ax.set_ylabel(f"change vs {NREF}-element mesh (%)")
ax.set_title("(b) Timoshenko FEM, L = 4 m\n(converges to its own limit)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("stage2_fig_mesh.png", dpi=160)

fig, ax = plt.subplots(figsize=(8, 4.3))
modes = np.arange(1, 7)
for L, c, mk in zip(LENGTHS, ["#2a78d6", "#eb6834", "#1baf7a"], ["o", "s", "^"]):
    fa = f_EB(L); ft = f_FEM(L, NCONV, timoshenko=True); fe = f_FEM(L, NCONV)
    ax.plot(modes, ft / fa, mk + "-", color=c, lw=2, ms=6, label=f"Timoshenko FEM, L = {L:.0f} m (L/r = {L/r:.0f})")
ax.plot(modes, f_FEM(10.0, NCONV) / f_EB(10.0), "k--", lw=1, label="Euler–Bernoulli FEM, all L")
ax.set_ylim(0, 1.1); style(ax)
ax.set_xlabel("mode number"); ax.set_ylabel("f_FEM / f_analytical (EB)")
ax.set_title("FEM vs Euler–Bernoulli analytical frequencies, liquid-filled cantilevers", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.tight_layout(); fig.savefig("stage2_fig_comparison.png", dpi=160)
print("\nSaved CSV tables and figures.")
