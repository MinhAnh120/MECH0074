"""
Stage 2 deliverable: "Compare 4, 10 and 20 m cantilevers with analytical
frequencies; check mesh sensitivity."

Uses the beam FEM in stage2_cantilever.py (keep in same folder):
  - EB FEM  : Euler-Bernoulli Hermite elements, consistent mass
  - TB FEM  : Timoshenko elements (Friedman & Kosmatka 1993), shear + rotary inertia
  - exact Timoshenko solution (transfer matrix) to verify the TB element and as the
    reference for its mesh study (TB converges as h^2, EB as h^4 - see order check)
Outputs:
  stage2_table_frequencies.csv  report table: first six frequencies at each length
  stage2_table_comparison.csv   FEM vs analytical, 3 lengths x 6 modes (+ TB verification)
  stage2_table_mesh.csv         mesh study, both element types, all lengths
  stage2_fig_mesh.png           mesh convergence (log-log)
  stage2_fig_comparison.png     FEM / analytical ratio per mode and length
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from stage2_cantilever import f_EB, f_FEM, f_TB_exact, m, EI, I, As, Do

LENGTHS = [4.0, 10.0, 20.0]
MESHES = [2, 4, 8, 16, 32, 64, 128]
NCONV = 64                     # mesh reported as "converged" (EB < 0.001 %, TB < 0.25 % vs exact)
r = np.sqrt(I / As)

# ---------------- mesh study ----------------
mesh_rows = []
mesh_err_EB, mesh_err_TB = {}, {}
for L in LENGTHS:
    fa = f_EB(L)
    fx = f_TB_exact(L)
    mesh_err_EB[L], mesh_err_TB[L] = [], []
    for ne in MESHES:
        fe = f_FEM(L, ne)
        ft = f_FEM(L, ne, timoshenko=True)
        eEB = (fe - fa) / fa * 100            # vs exact Euler-Bernoulli solution
        eTB = (ft - fx) / fx * 100            # vs exact Timoshenko solution
        mesh_err_EB[L].append(eEB); mesh_err_TB[L].append(eTB)
        for n in range(6):
            mesh_rows.append([L, ne, n + 1, fe[n], eEB[n], ft[n], eTB[n]])

with open("stage2_table_mesh.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["L_m", "elements", "mode", "f_EB_FEM_Hz", "err_vs_analytical_pct",
                "f_TB_FEM_Hz", "err_vs_exact_timoshenko_pct"])
    for row in mesh_rows:
        w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

# ---------------- comparison table ----------------
comp_rows = []
res = {}                       # L -> (EB theory, EB FEM, TB FEM, TB exact)
print("=" * 98)
print("COMPARISON: Euler-Bernoulli theory vs FEM (ne = %d)" % NCONV)
print("=" * 98)
for L in LENGTHS:
    fa = f_EB(L)
    fe = f_FEM(L, NCONV)
    ft = f_FEM(L, NCONV, timoshenko=True)
    fx = f_TB_exact(L)
    res[L] = (fa, fe, ft, fx)
    print(f"\nL = {L:.0f} m   (L/D = {L/Do:.1f}, slenderness L/r = {L/r:.1f}, "
          f"total mass = {m*L:,.0f} kg)")
    print(f"{'mode':>4} {'EB theory':>10} {'EB FEM':>10} {'diff %':>8} "
          f"{'Timo. FEM':>10} {'vs EB %':>8} {'Timo. exact':>12} {'FEM err %':>10}")
    for n in range(6):
        d1 = (fe[n] - fa[n]) / fa[n] * 100
        d2 = (ft[n] - fa[n]) / fa[n] * 100
        d3 = (ft[n] - fx[n]) / fx[n] * 100
        # round(...) + 0.0 turns round-off "-0.000" into "0.000" (display only; CSV keeps raw value)
        print(f"{n+1:>4} {fa[n]:>10.3f} {fe[n]:>10.3f} {round(d1, 3) + 0.0:>8.3f} "
              f"{ft[n]:>10.3f} {d2:>8.2f} {fx[n]:>12.3f} {round(d3, 3) + 0.0:>10.3f}")
        comp_rows.append([L, n + 1, fa[n], fe[n], d1, ft[n], d2, fx[n], d3])

with open("stage2_table_comparison.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["L_m", "mode", "f_EB_analytical_Hz", "f_EB_FEM_Hz", "EB_FEM_diff_pct",
                "f_Timoshenko_FEM_Hz", "Timoshenko_vs_EB_pct",
                "f_Timoshenko_exact_Hz", "Timoshenko_FEM_vs_exact_pct"])
    for row in comp_rows:
        w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

# ---------------- report table: first six bending frequencies at each length ----------------
with open("stage2_table_frequencies.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["mode"] + [h.format(L=f"{L:.0f}") for L in LENGTHS for h in
               ("EB_theory_{L}m_Hz", "FEM_EB_{L}m_Hz", "FEM_Timoshenko_{L}m_Hz", "Timoshenko_vs_EB_{L}m_pct")])
    for n in range(6):
        row = [n + 1]
        for L in LENGTHS:
            fa, fe, ft, _ = res[L]
            row += [f"{fa[n]:.3f}", f"{fe[n]:.3f}", f"{ft[n]:.3f}", f"{(ft[n]-fa[n])/fa[n]*100:.2f}"]
        w.writerow(row)

print("\n" + "=" * 98)
print("REPORT TABLE - first six bending frequencies (Hz), FEM ne = %d" % NCONV)
print("=" * 98)
print(f"{'':>5}" + "".join(f"{f'L = {L:.0f} m':^30}" for L in LENGTHS))
print(f"{'mode':>5}" + f"{'EB':>9}{'FEM-EB':>9}{'FEM-TB':>9}{'':>3}" * len(LENGTHS))
for n in range(6):
    print(f"{n+1:>5}" + "".join(f"{res[L][0][n]:>9.2f}{res[L][1][n]:>9.2f}{res[L][2][n]:>9.2f}{'':>3}"
                                for L in LENGTHS))
worst = max(abs(row[8]) for row in comp_rows)
print(f"\nTimoshenko FEM vs exact Timoshenko solution: max |error| = {worst:.3f} % (element verified)")

# ---------------- mesh table print (compact) ----------------
print("\n" + "=" * 86)
print("MESH SENSITIVITY (EB FEM, % error vs analytical; identical for all L)")
print("=" * 86)
print(f"{'ne':>5}" + "".join(f"{'mode '+str(n):>11}" for n in range(1, 7)))
for i, ne in enumerate(MESHES):
    print(f"{ne:>5}" + "".join(
        ("        n/a" if np.isnan(e) else f"{e:>11.4f}") for e in mesh_err_EB[10.0][i]))
print("Note: errors below ~1e-4 % (64+ elements) are floating-point round-off, not discretisation error.")
print("\nMESH SENSITIVITY (Timoshenko FEM, % error vs exact Timoshenko solution)")
for L in LENGTHS:
    print(f" L = {L:.0f} m")
    for i, ne in enumerate(MESHES):
        print(f"{ne:>5}" + "".join(
            ("        n/a" if np.isnan(e) else f"{e:>11.4f}") for e in mesh_err_TB[L][i]))

# Observed convergence order p = log2(e_ne / e_2ne), mode 6, ne 16 -> 32 (asymptotic range,
# above round-off). Expected rates (eigenvalue error ~ square of energy-norm error):
#   EB Hermite cubic : curvature w'' linear per element -> energy error O(h^2) -> f error O(h^4)
#   TB interdependent: shear strain (w' - psi) constant per element -> energy error O(h)
#                      -> f error O(h^2); recovers O(h^4) only as shear stiffness -> infinity
i16, i32 = MESHES.index(16), MESHES.index(32)
print("\nObserved convergence order, mode 6, ne 16 -> 32:")
print(f"  EB FEM (all L)  : {np.log2(mesh_err_EB[10.0][i16][5] / mesh_err_EB[10.0][i32][5]):.2f}  (theory: 4)")
for L in LENGTHS:
    p = np.log2(mesh_err_TB[L][i16][5] / mesh_err_TB[L][i32][5])
    print(f"  TB FEM, L = {L:>2.0f} m : {p:.2f}  (theory: 2, constant shear strain per element)")
print("Note: points with fewer than ~2 elements per half-wavelength (e.g. ne = 2-4 for modes 4-6)\n"
      "are pre-asymptotic, so their apparent orders are not meaningful.")

# ---------------- figures ----------------
COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MUTED = "#8a8984"
def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)

# plain decimal tick labels on log axes (0.001, 10, 100 instead of 10^-3, 10^1, 10^2)
PLAIN = FuncFormatter(lambda v, _: np.format_float_positional(v, trim="-"))
def plain_log_ticks(ax):
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(PLAIN)
    ax.set_xticks(MESHES); ax.xaxis.set_minor_locator(plt.NullLocator())   # tick at the meshes run

fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
ax = axs[0]
for n in range(6):
    e = np.array([mesh_err_EB[10.0][i][n] for i in range(len(MESHES))])
    ok = ~np.isnan(e) & (e > 1e-5)   # below ~1e-5 % is floating-point round-off
    # stop the line once the error rises again (round-off floor reached)
    idx = np.flatnonzero(ok)
    rising = np.flatnonzero(np.diff(e[idx]) > 0)
    if rising.size:
        ok[idx[rising[0] + 1]:] = False
    ax.loglog(np.array(MESHES)[ok], e[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"Mode {n+1}")
ax.axhline(0.1, color=MUTED, ls="--", lw=0.8); ax.text(2.1, 0.12, "0.1 %", color=MUTED, fontsize=8)
style(ax); plain_log_ticks(ax)
ax.set_xlabel("Number of Elements"); ax.set_ylabel("Error vs Analytical (%)")
ax.set_title("(a) Euler–Bernoulli FEM vs Exact Solution\n(Same for 4, 10, 20 m; Error ÷16 Each Time Elements Double)", loc="left", fontsize=10)

ax = axs[1]
for n in range(6):
    e = np.array([abs(mesh_err_TB[4.0][i][n]) for i in range(len(MESHES))])
    ok = ~np.isnan(e) & (e > 1e-9)
    ax.loglog(np.array(MESHES)[ok], e[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"Mode {n+1}")
ax.axhline(0.1, color=MUTED, ls="--", lw=0.8)
style(ax); plain_log_ticks(ax)
ax.set_xlabel("Number of Elements"); ax.set_ylabel("Error vs Exact Timoshenko (%)")
ax.set_title("(b) Timoshenko FEM, $L$ = 4 m (Slowest-Converging Length)\n(Error ÷4 Each Time Elements Double, Slowed by Shear)",
             loc="left", fontsize=10)
fig.tight_layout(rect=(0, 0.08, 1, 1))
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=9, frameon=True,
           fancybox=False, edgecolor=MUTED, framealpha=1)
fig.savefig("stage2_fig_mesh.png", dpi=160, bbox_inches="tight")

fig, ax = plt.subplots(figsize=(8, 4.3))
modes = np.arange(1, 7)
for L, c, mk in zip(LENGTHS, ["#2a78d6", "#eb6834", "#1baf7a"], ["o", "s", "^"]):
    fa = f_EB(L); ft = f_FEM(L, NCONV, timoshenko=True)
    ax.plot(modes, ft / fa, mk + "-", color=c, lw=2, ms=6, label=f"Timoshenko, $L$ = {L:.0f} m ($L/r$ = {L/r:.0f})")
ax.plot(modes, f_FEM(10.0, NCONV) / f_EB(10.0), "k--", lw=1, label="Euler–Bernoulli, All Lengths")
ax.set_ylim(0, 1.1); style(ax)
ax.set_xlabel("Mode Number"); ax.set_ylabel(r"$f_\mathrm{FEM}\,/\,f_\mathrm{analytical}$ (EB)")
ax.set_title("FEM vs Euler–Bernoulli Analytical Frequencies, Liquid-Filled Cantilevers", loc="left", fontsize=10)
ax.legend(title="FEM Model", title_fontsize=8,
          loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=8, frameon=True,
          fancybox=False, edgecolor=MUTED, framealpha=1)
fig.tight_layout(); fig.savefig("stage2_fig_comparison.png", dpi=160, bbox_inches="tight")
print("\nSaved CSV tables and figures.")
