"""MECH0074 LNG expansion loop - all Stage 2 work in one self-contained file.
Usage:
    python stage2_all.py             # run everything, then open all figure windows
    python stage2_all.py --no-show   # run everything, save outputs only (no windows)

All outputs (PNG figures, CSV tables) are written next to this file.
Requires: numpy, scipy, matplotlib

Goal (Stage 2): natural bending frequencies of the LNG pipe modelled as a cantilever
(fixed at one end, free at the other) for L = 4, 10 and 20 m, found two ways -
exact Euler-Bernoulli theory and a beam finite-element model (FEM) - then compared,
with a mesh study showing how many elements the FEM needs.

File layout (sections run in this order by main()):
  0. Package check     - stop early with a clear message if numpy/scipy/matplotlib are missing
  1. Shared parameters - all inputs from the brief + derived section properties
  2. EB hand calc      - analytical frequencies + mode shapes, step by step as in the hand calc
  3. Beam FEM          - element matrices, assembly, eigenvalue solve; quick mesh check
  4. Results           - comparison tables (CSV), mesh / refinement study, own-FEM figures
  5. ANSYS             - ANSYS mesh convergence (data typed in), per-length figures + CSVs
  6. Comparison        - ANSYS vs EB theory vs own FEM, table + two figures
  MAIN                 - runs 1-6, then lists which output files exist
"""
import argparse
import csv
import importlib
import math
import sys
from pathlib import Path


# ##########################################################################################
# 0. PACKAGE CHECK
# ##########################################################################################
# Runs before numpy/scipy/matplotlib are imported, so a missing package gives one clear
# message (with the pip command to fix it) instead of an ImportError traceback.
REQUIRED_PACKAGES = ["numpy", "scipy", "matplotlib"]


def check_packages():
    """Try to import every required package; report all missing ones and exit if any fail."""
    print(f"Python {sys.version.split()[0]} ({sys.executable})")
    missing = []
    for name in REQUIRED_PACKAGES:
        try:
            mod = importlib.import_module(name)
            print(f"  [OK     ] {name} {getattr(mod, '__version__', '')}")
        except ImportError:
            print(f"  [MISSING] {name}")
            missing.append(name)
    if missing:
        print(f"\nMissing package(s): {', '.join(missing)}")
        print("Activate the project venv or run:")
        print(f'  "{sys.executable}" -m pip install {" ".join(missing)}')
        sys.exit(1)


check_packages()

# Safe to import now that the check has passed.
# brentq = root finder (beta_n L roots), eigh = generalised symmetric eigen-solver (FEM K, M).
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent      # outputs go next to this file, whatever the cwd


# ##########################################################################################
# 1. SHARED PARAMETERS  (was params.py)
# ##########################################################################################
# Supplied values are taken from the brief. Anything marked TBC is still to be
# agreed by the group; log any change and its reason in the assumptions log on Teams.

# ---------------------------------------------------------------- pipe geometry (brief)
D_OUT = 0.813            # outer diameter, m
D_IN = 0.79394           # inner diameter, m
T_WALL = (D_OUT - D_IN) / 2          # wall thickness, m (~9.53 mm)
BEND_RADIUS = 1.4        # radius of curvature of each of the four bends, m
N_BENDS = 4

# ---------------------------------------------------------------- material: 304 stainless (brief)
RHO_STEEL = 8000.0       # kg/m^3
E_STEEL = 213.233e9      # Young's modulus (NIST fit), Pa
NU_STEEL = 0.29          # Poisson's ratio, 304 SS (TBC: typical datasheet value, confirm source)
ALPHA = 1.4e-5           # thermal expansivity, 1/K (confirmed on Moodle)
T_INSTALL = 293.0        # K
T_OPERATE = 111.0        # K
DELTA_T = T_INSTALL - T_OPERATE      # 182 K, uniform slow cooldown

# ---------------------------------------------------------------- fluid: LNG (brief)
RHO_LNG = 455.0          # kg/m^3
C_BULK = 1100.0          # bulk sound speed, m/s
U0 = 1.5                 # initial flow velocity, m/s

# ---------------------------------------------------------------- pressure and valve (brief)
P0_ABS = 0.7e6           # initial absolute pressure, Pa
P_EXT_ABS = 0.101325e6   # external absolute pressure, Pa
T_CLOSE = 0.02           # linear valve closure time, s
L_BENCHMARK = 100.0      # Stage 1 benchmark length, m

# ---------------------------------------------------------------- wave speed choice (TBC, group decision)
USE_KORTEWEG = False
K_LNG = RHO_LNG * C_BULK**2                      # bulk modulus, Pa
C_KORTEWEG = C_BULK / math.sqrt(1 + K_LNG * D_IN / (E_STEEL * T_WALL))   # ~998 m/s
C_WAVE = C_KORTEWEG if USE_KORTEWEG else C_BULK

# ---------------------------------------------------------------- section properties (derived)
A_STEEL = math.pi / 4 * (D_OUT**2 - D_IN**2)      # m^2 (~0.02406)
A_BORE = math.pi / 4 * D_IN**2                    # m^2 (~0.4951)
I_SECTION = math.pi / 64 * (D_OUT**4 - D_IN**4)  # m^4 (~1.941e-3)
M_STEEL = RHO_STEEL * A_STEEL                     # kg/m (~192.4)
M_LNG = RHO_LNG * A_BORE                          # kg/m (~225.3)
M_TOTAL = M_STEEL + M_LNG                         # kg/m (~417.7), fluid as added mass

# ---------------------------------------------------------------- hand-check values (derived)
DP_JOUKOWSKY = RHO_LNG * C_WAVE * U0              # surge, Pa (~0.75 MPa with bulk c)
T_RETURN = 2 * L_BENCHMARK / C_WAVE               # wave return time 2L/c, s
T_PERIOD = 4 * L_BENCHMARK / C_WAVE               # pressure period 4L/c, s
F_HYDRAULIC = 1 / T_PERIOD                        # fundamental, Hz
DL_PER_100M = ALPHA * 100.0 * DELTA_T             # free contraction per 100 m, m (~0.255)
SIGMA_RESTRAINED = E_STEEL * ALPHA * DELTA_T      # fully restrained stress, Pa (~543 MPa)

# ---------------------------------------------------------------- loop geometry (TBC, Kabir)
LOOP_HEIGHT = None       # leg length of the U-loop, m
LOOP_WIDTH = None        # distance between legs, m
GUIDE_POSITIONS = []     # distances of guides from the loop, m

# ---------------------------------------------------------------- damping (TBC, Beau)
DAMPING_RATIO = None     # modal damping ratio for Stage 4 baseline


def print_params():
    """Step 1 output: print the key derived values so they can be checked against the hand calc."""
    print(f"Wall thickness     {T_WALL*1e3:.2f} mm")
    print(f"Mass per length     {M_TOTAL:.1f} kg/m (steel {M_STEEL:.1f}, LNG {M_LNG:.1f})")
    print(f"Second moment       {I_SECTION:.4e} m^4")
    print(f"Wave speed used     {C_WAVE:.1f} m/s (Korteweg would give {C_KORTEWEG:.1f})")
    print(f"Joukowsky surge     {DP_JOUKOWSKY/1e6:.3f} MPa")
    print(f"2L/c, 4L/c          {T_RETURN:.3f} s, {T_PERIOD:.3f} s ({F_HYDRAULIC:.2f} Hz)")
    print(f"Contraction/100 m   {DL_PER_100M*1e3:.1f} mm")
    print(f"Restrained stress   {SIGMA_RESTRAINED/1e6:.1f} MPa")


# ##########################################################################################
# 2. EULER-BERNOULLI HAND CALCULATION  (was stage2_euler_bernoulli.py)
# ##########################################################################################
# Follows the hand-calc steps one-to-one:
#   Step 4  roots of 1 + cos(bL) cosh(bL) = 0      (bisection, + asymptotic check)
#   Step 6  section properties, mass per length, sqrt(EI/m)
#   Step 7  worked example (L = 10 m, mode 1)
#   Step 8  frequency table  f_n = C_L * (beta_n L)^2,  C_L = sqrt(EI/m) / (2 pi L^2)
#   Step 9  mode shapes W_n(x) and node count
#   Step 10 checks (1/L^2 scaling, mode ratios, fluid-mass effect)
# References: Rao, Mechanical Vibrations (beam vibration chapter);
#             Blevins, Formulas for Natural Frequency and Mode Shape (cantilever table).

def run_euler_bernoulli():
    """Analytical cantilever frequencies and mode shapes, printed step by step.
    Self-contained (own local I, m, EI) so it can be checked against the hand calc on its own.
    Output: stage2_EB_mode_shapes.png + printed tables."""
    import matplotlib.pyplot as plt   # imported here so main() can pick the backend first

    Do, Di, E, rho_s, rho_f = D_OUT, D_IN, E_STEEL, RHO_STEEL, RHO_LNG
    LENGTHS = [4.0, 10.0, 20.0]   # m
    N_MODES = 6

    # ---------------- STEP 4 - roots of g(x) = 1 + cos(x) cosh(x) = 0 ----------------
    # Clamped-free boundary conditions reduce to this characteristic equation. Its roots
    # beta_n L are the only "mode-dependent" numbers; they are the same for every length.
    def g(x):
        return 1.0 + math.cos(x) * math.cosh(x)

    def bisection(f, a, b, tol=1e-12, max_iter=200):
        """Plain bisection: needs f(a) and f(b) of opposite sign."""
        fa = f(a)
        if fa * f(b) > 0:
            raise ValueError("No sign change in bracket")
        for _ in range(max_iter):
            c = 0.5 * (a + b)
            fc = f(c)
            if abs(b - a) < tol:
                break
            if fa * fc < 0:
                b = c
            else:
                a, fa = c, fc
        return 0.5 * (a + b)

    def cantilever_roots(n):
        # Each root lies near (2k-1)*pi/2, so bracket +/- 0.6 around it.
        return [bisection(g, (2*k - 1) * math.pi / 2 - 0.6, (2*k - 1) * math.pi / 2 + 0.6)
                for k in range(1, n + 1)]

    BL = cantilever_roots(N_MODES)

    print("STEP 4 - Roots of 1 + cos(bL) cosh(bL) = 0")
    print(f"{'n':>3} {'beta_n L':>10} {'(2n-1)pi/2':>12} {'(beta_n L)^2':>14}")
    for n, b in enumerate(BL, 1):
        print(f"{n:>3} {b:>10.4f} {(2*n-1)*math.pi/2:>12.4f} {b*b:>14.4f}")
    print("Trial values for root 1:  g(1.8) = %.4f, g(1.875) = %.5f, g(1.8751) = %.6f"
          % (g(1.8), g(1.875), g(1.8751)))

    # ---------------- STEP 6 - section properties ----------------
    # Hollow circular section. Only the steel gives stiffness (EI); the LNG inside moves
    # with the pipe, so it adds mass per length but no stiffness.
    I   =math.pi * (Do**4 - Di**4) / 64          # m^4
    As  = math.pi * (Do**2 - Di**2) / 4           # m^2, steel
    Af  = math.pi * Di**2 / 4                     # m^2, fluid
    m_s = rho_s * As                              # kg/m
    m_f = rho_f * Af                              # kg/m
    m   = m_s + m_f                               # kg/m  (fluid = added mass, no stiffness)
    EI  = E * I
    root_EI_m = math.sqrt(EI / m)

    print("\nSTEP 6 - Section properties")
    print(f"  Do^4 = {Do**4:.6f} m4,  Di^4 = {Di**4:.6f} m4,  difference = {Do**4 - Di**4:.6f}")
    print(f"  I         = {I:.4e} m4")
    print(f"  EI        = {EI:.4e} N m2")
    print(f"  m_steel   = {m_s:.2f} kg/m   (A_s = {As:.6f} m2)")
    print(f"  m_LNG     = {m_f:.2f} kg/m   (A_f = {Af:.5f} m2)")
    print(f"  m total   = {m:.2f} kg/m")
    print(f"  EI/m      = {EI/m:.4e} m4/s2")
    print(f"  sqrt(EI/m)= {root_EI_m:.2f} m2/s")

    # ---------------- STEP 5/7/8 - frequencies ----------------
    # f_n = (beta_n L)^2 / (2 pi L^2) * sqrt(EI/m). C_L holds everything except (beta_n L)^2,
    # so each length needs one number; freqs() can swap the mass (used for the steel-only check).
    def C_L(L):
        return root_EI_m / (2 * math.pi * L**2)

    def freqs(L, mass=m):
        return [b**2 / (2 * math.pi * L**2) * math.sqrt(EI / mass) for b in BL]

    print("\nSTEP 7 - Worked example, L = 10 m, mode 1")
    print(f"  f1 = {BL[0]**2:.4f} / (2 pi x 10^2) x {root_EI_m:.2f}"
          f" = {BL[0]**2/(2*math.pi*100):.7f} x {root_EI_m:.2f} = {freqs(10)[0]:.3f} Hz")

    print("\nSTEP 8 - Frequencies f_n = C_L (beta_n L)^2   [Hz]")
    for L in LENGTHS:
        print(f"  C_L({L:>4.0f} m) = {C_L(L):.5f} Hz")
    header = f"{'mode':>5} {'(bL)^2':>9}" + "".join(f"{f'L={L:.0f} m':>12}" for L in LENGTHS)
    print(header)
    table = {L: freqs(L) for L in LENGTHS}
    for n in range(N_MODES):
        row = f"{n+1:>5} {BL[n]**2:>9.4f}" + "".join(f"{table[L][n]:>12.3f}" for L in LENGTHS)
        print(row)

    # ---------------- STEP 9 - mode shapes ----------------
    # W_n(x) = cosh - cos - sigma_n (sinh - sin), plotted against x/L so one figure covers all
    # lengths. Node count check: mode n should cross zero n-1 times (excluding the clamp).
    def sigma(bL):
        return (math.cosh(bL) + math.cos(bL)) / (math.sinh(bL) + math.sin(bL))

    def mode_shape(n, xi):
        """Normalised mode shape W_n at xi = x/L (tip value scaled to +1)."""
        b = BL[n - 1]
        s = sigma(b)
        z = b * xi
        W = np.cosh(z) - np.cos(z) - s * (np.sinh(z) - np.sin(z))
        return W / W[-1]

    xi = np.linspace(0, 1, 2001)          # x/L, fine enough for smooth curves and node counting
    print("\nSTEP 9 - Mode shape constants and nodes")
    for n in range(1, N_MODES + 1):
        W = mode_shape(n, xi)
        nodes = np.sum(np.diff(np.sign(W[5:])) != 0)     # skip the clamped end itself
        print(f"  mode {n}: sigma = {sigma(BL[n-1]):.4f}, interior nodes = {nodes} (expected {n-1})")

    # Figure: first six mode shapes (stage2_EB_mode_shapes.png)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors =["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
    for n in range(1, N_MODES + 1):
        ax.plot(xi, mode_shape(n, xi), lw=2, color=colors[n-1], label=f"Mode {n}")
    ax.axhline(0, color="#8a8984", lw=0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(r"$x/L$  (0 = Fixed End, 1 = Free End)")
    ax.set_ylabel(r"Normalised Deflection, $W_n(x)\,/\,W_n(L)$")
    ax.set_title(r"Euler–Bernoulli Cantilever Mode Shapes (Independent of $L$)", loc="left")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True,
              fancybox=False, edgecolor="#8a8984", framealpha=1, fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "stage2_EB_mode_shapes.png", dpi=150, bbox_inches="tight")

    # ---------------- STEP 10 - checks ----------------
    # Sanity checks on the results:
    #   1/L^2 scaling  - doubling L should quarter f; 10 m -> 4 m should multiply f by 6.25
    #   mode ratios    - (beta_n/beta_1)^2, standard cantilever values 1 : 6.27 : 17.55 ...
    #   fluid effect   - empty pipe is stiffer per unit mass: f rises by sqrt(m_total/m_steel)
    print("\nSTEP 10 - Checks")
    f10, f20, f4 = table[10][0], table[20][0], table[4][0]
    print(f"  1/L^2 scaling: f1(10)/4 = {f10/4:.3f} vs f1(20) = {f20:.3f};"
          f"  f1(10)*6.25 = {f10*6.25:.2f} vs f1(4) = {f4:.2f}")
    print("  Mode ratios f_n/f_1 :", " : ".join(f"{(b/BL[0])**2:.3f}" for b in BL))
    f_steel = freqs(10, mass=m_s)[0]
    print(f"  Steel only, L = 10 m: f1 = {f_steel:.3f} Hz;  ratio to filled = {f_steel/f10:.3f}"
          f"  (sqrt(m/m_s) = {math.sqrt(m/m_s):.3f})")
    print("\nSaved stage2_EB_mode_shapes.png")


# ##########################################################################################
# 3. BEAM FEM  (was stage2_cantilever.py)
# ##########################################################################################
# Analytical : Euler-Bernoulli cantilever, f_n = (beta_n L)^2 / (2 pi L^2) * sqrt(EI/m)
#              beta_n L roots of 1 + cos(bL) cosh(bL) = 0  (Blevins; Rao)
# FEM        : 2-node Euler-Bernoulli (Hermite cubic) beam elements, consistent mass.
# Fluid      : internal LNG treated as mass moving with the pipe (added to m, not to EI).
# These are module-level because section 4 (results) uses f_EB, f_FEM, m and Do.

Do = D_OUT                               # outer diameter, m (used for L/D in printouts)
m = M_TOTAL                              # mass per unit length (pipe + contents)
EI = E_STEEL * I_SECTION                 # bending stiffness, N m^2 (steel only)


# ---- analytical reference (same theory as section 2, compact form used by the FEM comparison) ----
def bl_roots(n=6):
    """First n roots beta_n L of 1 + cos cosh = 0, each bracketed around (k - 1/2) pi."""
    g = lambda x: 1 + np.cos(x) * np.cosh(x)
    return np.array([brentq(g, (k - 0.5) * np.pi - 0.6, (k - 0.5) * np.pi + 0.6) for k in range(1, n + 1)])


BL = bl_roots()                          # beta_n L, modes 1-6 (1.875, 4.694, 7.855, ...)


def f_EB(L):
    """Exact Euler-Bernoulli cantilever frequencies (Hz), modes 1-6, for length L."""
    return BL**2 / (2 * np.pi * L**2) * np.sqrt(EI / m)


def element(Le):
    """Euler-Bernoulli Hermite cubic element: stiffness and consistent mass.
    DOFs per element: [w1, theta1, w2, theta2] (deflection and rotation at each end node).
    k  = standard beam stiffness, EI/L^3 * [...]
    mt = consistent mass, mL/420 * [...] (same cubic shape functions as k, not lumped)."""
    L = Le
    k = EI / L**3 * np.array([
        [12,    6*L,    -12,   6*L],
        [6*L,   4*L**2, -6*L,  2*L**2],
        [-12,  -6*L,     12,  -6*L],
        [6*L,   2*L**2, -6*L,  4*L**2]])
    mt = m * L / 420 * np.array([
        [156,    22*L,    54,    -13*L],
        [22*L,   4*L**2,  13*L,  -3*L**2],
        [54,     13*L,    156,   -22*L],
        [-13*L, -3*L**2, -22*L,   4*L**2]])
    return k, mt


def f_FEM(L, ne, nmodes=6):
    """FEM cantilever frequencies (Hz) for length L split into ne equal elements."""
    # 1. Global matrices: ne + 1 nodes x 2 DOFs (w, theta)
    nd = 2 * (ne + 1)
    K, M = np.zeros((nd, nd)), np.zeros((nd, nd))
    ke, me = element(L / ne)              # all elements identical (uniform pipe)
    # 2. Assembly: element e joins nodes e and e+1 -> global DOFs 2e .. 2e+3 (overlap = shared node)
    for e in range(ne):
        s = slice(2*e, 2*e + 4)
        K[s, s] += ke; M[s, s] += me
    # 3. Boundary conditions: delete the 2 DOFs of node 0 (clamped); tip node is left free
    K, M = K[2:, 2:], M[2:, 2:]           # clamp node 0 (w = 0, theta = 0)
    # 4. Solve K phi = omega^2 M phi, keep the lowest modes
    w2 = eigh(K, M, eigvals_only=True)[:nmodes]
    # 5. omega^2 -> f in Hz. A coarse mesh has only 2*ne DOFs, so it cannot give all nmodes:
    #    missing modes are left as NaN (e.g. ne = 2 -> modes 5 and 6 = NaN).
    f = np.full(nmodes, np.nan); f[:len(w2)] = np.sqrt(w2) / (2 * np.pi)
    return f


def run_cantilever_checks():
    """Quick FEM checks printed to screen (no files):
    - section properties used by the FEM (should match section 2)
    - flow-velocity check: dimensionless velocity u* << 1 means flow does not change the
      frequencies, which justifies treating the LNG as added mass only
    - short mesh study per length, % error vs the exact frequencies"""
    print("---- Section properties ----")
    print(f"A_steel = {A_STEEL:.5f} m2   A_fluid = {A_BORE:.5f} m2")
    print(f"I       = {I_SECTION:.5e} m4   EI = {EI:.4e} N m2")
    print(f"m_steel = {M_STEEL:.2f} kg/m  m_fluid = {M_LNG:.2f} kg/m  m_total = {m:.2f} kg/m")
    print("beta_n L =", np.round(BL, 4))

    print("\n---- Flow-velocity check (Paidoussis dimensionless velocity) ----")
    for L in (4, 10, 20):
        u = U0 * L * np.sqrt(M_LNG / EI)
        print(f"L = {L:2d} m : u* = U L sqrt(M_f/EI) = {u:.2e}")

    for L in (4, 10, 20):
        fa = f_EB(L)
        print(f"\n==== L = {L} m  (L/D = {L/Do:.1f}) ====")
        print("Analytical EB (Hz):", np.round(fa, 3))
        print("Mesh study, EB FEM: % error vs analytical")
        for ne in (2, 4, 8, 16, 32, 64):
            ff = f_FEM(L, ne)
            err = (ff - fa) / fa * 100
            print(f"  ne={ne:3d}  " + "  ".join(("      n/a" if np.isnan(e) else f"{e:+9.4f}") for e in err))


# ##########################################################################################
# 4. RESULTS: COMPARISON + MESH STUDY  (was stage2_results.py)
# ##########################################################################################
# Deliverable: "Compare 4, 10 and 20 m cantilevers with analytical frequencies;
# check mesh sensitivity."
# Outputs:
#   stage2_table_frequencies.csv  report table: first six frequencies at each length
#   stage2_table_comparison.csv   FEM vs analytical, 3 lengths x 6 modes
#   stage2_table_mesh.csv         mesh study, all lengths
#   stage2_fig_mesh_EB.png        mesh convergence, Euler-Bernoulli FEM (log-log)
#   stage2_fig_refinement_EB.png  % change per mesh doubling, Euler-Bernoulli FEM
#   stage2_fig_mesh_EB_zeroed.png        as mesh_EB, round-off points set to 0 (hollow markers)
#   stage2_fig_refinement_EB_zeroed.png  as refinement_EB, round-off points set to 0 (hollow markers)

def run_results():
    """Deliverable tables and figures.
    Order: mesh study -> comparison table -> report table -> printed summaries -> figures.
    Note: % error and % change are identical for every L (both exact and FEM frequencies
    scale as 1/L^2, so the ratio cancels L). Figures therefore plot L = 10 m only and are
    titled "Same for 4, 10, 20 m"."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    LENGTHS = [4.0, 10.0, 20.0]
    MESHES = [2, 4, 8, 16, 32, 64, 128, 256]   # element counts, doubling each time
    NCONV = 64                     # mesh reported as "converged" (EB < 0.001 % vs analytical)

    # ---------------- mesh study ----------------
    # For every length and mesh, store two measures per mode:
    #   eEB = signed % error vs exact theory (only possible because an exact answer exists)
    #   cEB = % change from the previous, half-as-fine mesh (the practical convergence check)
    # Results -> stage2_table_mesh.csv
    mesh_rows = []
    mesh_err_EB = {}
    mesh_chg_EB = {}               # % change from the previous (half-as-fine) mesh
    for L in LENGTHS:
        fa = f_EB(L)
        mesh_err_EB[L] = []
        mesh_chg_EB[L] = []
        fe_prev = None
        for ne in MESHES:
            fe = f_FEM(L, ne)
            eEB = (fe - fa) / fa * 100            # vs exact Euler-Bernoulli solution
            # before vs after refinement: needs no exact solution
            cEB = np.full(6, np.nan) if fe_prev is None else np.abs(fe - fe_prev) / fe * 100
            mesh_err_EB[L].append(eEB)
            mesh_chg_EB[L].append(cEB)
            fe_prev = fe
            for n in range(6):
                mesh_rows.append([L, ne, n + 1, fe[n], eEB[n], cEB[n]])

    # write the mesh study: one row per (length, mesh, mode); 6 significant figures
    with open(HERE / "stage2_table_mesh.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["L_m", "elements", "mode", "f_EB_FEM_Hz", "err_vs_analytical_pct",
                    "EB_change_vs_previous_mesh_pct"])
        for row in mesh_rows:
            w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

    # ---------------- comparison table ----------------
    # Exact theory vs FEM at the converged mesh (NCONV elements), each length and mode.
    # Printed and saved -> stage2_table_comparison.csv. `res` is reused for the report table.
    comp_rows = []
    res = {}                       # L -> (EB theory, EB FEM)
    print("=" * 60)
    print("COMPARISON: Euler-Bernoulli theory vs FEM (ne = %d)" % NCONV)
    print("=" * 60)
    for L in LENGTHS:
        fa = f_EB(L)
        fe = f_FEM(L, NCONV)
        res[L] = (fa, fe)
        print(f"\nL = {L:.0f} m   (L/D = {L/Do:.1f}, total mass = {m*L:,.0f} kg)")
        print(f"{'mode':>4} {'EB theory':>10} {'EB FEM':>10} {'diff %':>8}")
        for n in range(6):
            d1 = (fe[n] - fa[n]) / fa[n] * 100
            # round(...) + 0.0 turns round-off "-0.000" into "0.000" (display only; CSV keeps raw value)
            print(f"{n+1:>4} {fa[n]:>10.3f} {fe[n]:>10.3f} {round(d1, 3) + 0.0:>8.3f}")
            comp_rows.append([L, n + 1, fa[n], fe[n], d1])

    with open(HERE / "stage2_table_comparison.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["L_m", "mode", "f_EB_analytical_Hz", "f_EB_FEM_Hz", "EB_FEM_diff_pct"])
        for row in comp_rows:
            w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

    # ---------------- report table: first six bending frequencies at each length ----------------
    # Wide format (one row per mode, theory + FEM columns per length) ready to paste into the
    # report -> stage2_table_frequencies.csv, also printed below.
    with open(HERE / "stage2_table_frequencies.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mode"] + [h.format(L=f"{L:.0f}") for L in LENGTHS for h in
                   ("EB_theory_{L}m_Hz", "FEM_EB_{L}m_Hz")])
        for n in range(6):
            row = [n + 1]
            for L in LENGTHS:
                fa, fe = res[L]
                row += [f"{fa[n]:.3f}", f"{fe[n]:.3f}"]
            w.writerow(row)

    print("\n" + "=" * 60)
    print("REPORT TABLE - first six bending frequencies (Hz), FEM ne = %d" % NCONV)
    print("=" * 60)
    print(f"{'':>5}" + "".join(f"{f'L = {L:.0f} m':^21}" for L in LENGTHS))
    print(f"{'mode':>5}" + f"{'EB':>9}{'FEM-EB':>9}{'':>3}" * len(LENGTHS))
    for n in range(6):
        print(f"{n+1:>5}" + "".join(f"{res[L][0][n]:>9.2f}{res[L][1][n]:>9.2f}{'':>3}"
                                    for L in LENGTHS))

    # ---------------- mesh table print (compact) ----------------
    # Screen summary of the mesh study, L = 10 m only (identical for all L, see docstring).
    print("\n" + "=" * 86)
    print("MESH SENSITIVITY (EB FEM, % error vs analytical; identical for all L)")
    print("=" * 86)
    print(f"{'ne':>5}" + "".join(f"{'mode '+str(n):>11}" for n in range(1, 7)))
    for i, ne in enumerate(MESHES):
        print(f"{ne:>5}" + "".join(
            ("        n/a" if np.isnan(e) else f"{e:>11.4f}") for e in mesh_err_EB[10.0][i]))
    print("Note: errors below ~1e-4 % (64+ elements) are floating-point round-off, not discretisation error.")

    # Observed convergence order p = log2(e_ne / e_2ne), mode 6, ne 16 -> 32 (asymptotic range,
    # above round-off). Expected rate (eigenvalue error ~ square of energy-norm error):
    #   EB Hermite cubic : curvature w'' linear per element -> energy error O(h^2) -> f error O(h^4)
    i16, i32 = MESHES.index(16), MESHES.index(32)
    print("\nObserved convergence order, mode 6, ne 16 -> 32:")
    print(f"  EB FEM (all L)  : {np.log2(mesh_err_EB[10.0][i16][5] / mesh_err_EB[10.0][i32][5]):.2f}  (theory: 4)")
    print("Note: points with fewer than ~2 elements per half-wavelength (e.g. ne = 2-4 for modes 4-6)\n"
          "are pre-asymptotic, so their apparent orders are not meaningful.")

    # ---------------- refinement study: before vs after each mesh doubling ----------------
    # Change = |f(ne) - f(ne/2)| / f(ne). Needs no exact solution, so it is the check used
    # in practice. For a method of order p the change overstates the remaining error of the
    # finer mesh by (2^p - 1): x15 for EB (p = 4).
    def print_change(title, rows):
        """Print the % change table, one line per doubling step (row 0 has no previous mesh)."""
        print("\n" + title)
        print(f"{'ne':>11}" + "".join(f"{'mode '+str(n):>11}" for n in range(1, 7)))
        for i in range(1, len(MESHES)):
            print(f"{MESHES[i-1]:>4} -> {MESHES[i]:<4}" + "".join(
                ("        n/a" if np.isnan(c) else f"{c:>11.4f}") for c in rows[i]))

    print("\n" + "=" * 86)
    print("REFINEMENT STUDY: % change in frequency from previous mesh (before vs after doubling)")
    print("=" * 86)
    print_change("EB FEM (identical for all L)", mesh_chg_EB[10.0])

    # ---------------- figures ----------------
    # Four figures, all L = 10 m (same for every L):
    #   stage2_fig_mesh_EB.png               |error| vs mesh, every point plotted
    #   stage2_fig_refinement_EB.png         % change per doubling, every point plotted
    #   stage2_fig_mesh_EB_zeroed.png        as above, round-off points set to 0
    #   stage2_fig_refinement_EB_zeroed.png  as above, round-off points set to 0
    # Shared styling helpers first, then each figure.
    COLS =["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
    MUTED = "#8a8984"

    def style(ax):
        """Common look: no top/right spines, light grid behind the data."""
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)

    # plain decimal tick labels on log axes (0.001, 10, 100 instead of 10^-3, 10^1, 10^2)
    PLAIN = FuncFormatter(lambda v, _: np.format_float_positional(v, trim="-"))

    def plain_log_ticks(ax):
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(PLAIN)
        ax.set_xticks(MESHES); ax.xaxis.set_minor_locator(plt.NullLocator())   # tick at the meshes run

    # NOTE: trim_roundoff is no longer called (the figures now plot every point, or zero the
    # round-off points via roundoff_mask below). Kept in case the old "hide round-off" look is wanted.
    def trim_roundoff(e, floor=1e-5, rise_after=0):
        """Mask of points to plot: stop each line for good at the round-off floor, i.e. at the
        first value below `floor` or the first rise after `rise_after` elements (round-off
        grows as the matrices get worse conditioned, so later points are not discretisation error)."""
        x = np.array(MESHES[-len(e):])
        ok = ~np.isnan(e)
        below = np.flatnonzero(ok & (e <= floor))
        if below.size:
            ok[below[0]:] = False
        idx = np.flatnonzero(ok)
        rising = idx[1:][np.diff(e[idx]) > 0]
        rising = rising[x[rising] > rise_after]
        if rising.size:
            ok[rising[0]:] = False
        return ok

    def mode_legend(ax):
        """One-row legend for modes 1-6, placed below the axes so it never covers data."""
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=6, fontsize=8, frameon=True,
                  fancybox=False, edgecolor=MUTED, framealpha=1)

    # ---- Figure: mesh convergence (stage2_fig_mesh_EB.png) ----
    # |% error| vs number of elements, log-log. Straight lines of slope -4 = 4th-order convergence
    # (error /16 per doubling). Points at 0 error or NaN (ne = 2, modes 5-6) cannot be drawn.
    # At fine meshes the lines turn up: floating-point round-off, not real error.
    # Dashed line = 0.1 % accuracy target.
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for n in range(6):
        e = np.array([mesh_err_EB[10.0][i][n] for i in range(len(MESHES))])
        e = np.abs(e)                    # plot all points; |error| so round-off negatives show on log axis
        ok = np.isfinite(e) & (e > 0)
        ax.loglog(np.array(MESHES)[ok], e[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"Mode {n+1}")
    ax.axhline(0.1, color=MUTED, ls="--", lw=0.8); ax.text(2.1, 0.12, "0.1 %", color=MUTED, fontsize=8)
    style(ax); plain_log_ticks(ax)
    ax.set_xlabel("Number of Elements"); ax.set_ylabel("Error vs Analytical (%)")
    ax.set_title("Euler–Bernoulli FEM vs Exact Solution\n(Same for 4, 10, 20 m)", loc="left", fontsize=10)
    mode_legend(ax)
    fig.tight_layout(); fig.savefig(HERE / "stage2_fig_mesh_EB.png", dpi=160, bbox_inches="tight")

    # ---- Figure: refinement study (stage2_fig_refinement_EB.png) ----
    # % change in frequency from the previous mesh, plotted at each doubling step ("a→b").
    # Same behaviour as the error plot, but needs no exact solution. First step for modes 5-6
    # is missing because ne = 2 has no result for them.
    def plot_change(rows, title, fname):
        fig, ax = plt.subplots(figsize=(7, 4.3))
        x = np.array(MESHES[1:])
        for n in range(6):
            c = np.array([rows[i][n] for i in range(1, len(MESHES))])
            ok = np.isfinite(c) & (c > 0)    # plot all points, including round-off at fine meshes
            ax.loglog(x[ok], c[ok], "o-", color=COLS[n], lw=2, ms=4, label=f"Mode {n+1}")
        ax.axhline(0.1, color=MUTED, ls="--", lw=0.8)
        style(ax); plain_log_ticks(ax)
        ax.set_xticks(x); ax.set_xlim(x[0] / 1.25, x[-1] * 1.25)
        ax.text(1.01, 0.1, "0.1 %", color=MUTED, fontsize=8, ha="left", va="center",
                transform=ax.get_yaxis_transform())   # just outside the right edge: never overlaps data
        ax.set_xticklabels([f"{a}→{b}" for a, b in zip(MESHES[:-1], MESHES[1:])], fontsize=8)
        ax.set_xlabel("Mesh Refinement (Number of Elements, Before → After)")
        ax.set_ylabel("Change in Frequency (%)")
        ax.set_title(title, loc="left", fontsize=10)
        mode_legend(ax)
        fig.tight_layout(); fig.savefig(HERE / fname, dpi=160, bbox_inches="tight")

    plot_change(mesh_chg_EB[10.0],
                "Euler–Bernoulli FEM: Change in Frequency per Mesh Doubling\n(Same for 4, 10, 20 m)",
                "stage2_fig_refinement_EB.png")

    # ---- round-off zeroed versions ----
    # Discretisation error falls ~16x per doubling (order 4). Once past `after` elements, the first
    # point that rises or falls by less than 4x is round-off; it and all later points are set to 0.
    # NOTE (round-off zeroing): the figures carry NO in-plot label or title text for this.
    #   - Hollow markers on the "0" row = values zeroed as floating-point round-off (not real error).
    #   - Dotted line = link from the last genuine point to the first zeroed point.
    #   - Explain this in the report figure caption, e.g. "Hollow markers: values below
    #     floating-point round-off (~1e-6 %), set to 0."
    #   - Zeroed with current data: error  - mode 1 at 64/128/256, mode 2 at 128/256, mode 3 at 256
    #                               change - mode 1 at 64→128/128→256, mode 2 at 128→256
    #   - Original unzeroed values are still in stage2_table_mesh.csv and the non-"_zeroed" figures.
    ZERO_Y = 1e-7                          # where "0" is drawn on the log axis

    def roundoff_mask(v, x, after=8, min_ratio=4.0):
        """True for points treated as round-off. Scans each line left to right; past `after`
        elements, the first step that falls by less than `min_ratio` (or rises) marks that point
        and all later ones. `after` = 8 skips coarse pre-asymptotic meshes, where uneven drops are real."""
        bad =np.zeros(len(v), bool)
        idx = np.flatnonzero(np.isfinite(v))
        for a, b in zip(idx[:-1], idx[1:]):
            if x[b] > after and v[b] > 0 and v[a] / v[b] < min_ratio:
                bad[b:] = True
                break
        return bad

    def plot_zeroed(x, series, xlabels, xlabel, ylabel, title, fname):
        """Log-log plot of six series (one per mode). Genuine points: solid line + filled markers.
        Round-off points: hollow markers on a "0" row at ZERO_Y (log axes cannot show 0), joined
        by a dotted line. Y ticks are set by hand so the bottom tick reads "0"."""
        fig, ax = plt.subplots(figsize=(7, 4.6))
        jit = np.linspace(-0.06, 0.06, 6)                      # spread overlapping zero markers
        for n, v in enumerate(series):
            bad = roundoff_mask(v, x)
            good = np.isfinite(v) & ~bad
            ax.loglog(x[good], v[good], "o-", color=COLS[n], lw=2, ms=4, label=f"Mode {n+1}")
            if bad.any():
                xb = x[bad] * 2 ** jit[n]
                last = np.flatnonzero(good)[-1]
                ax.plot([x[last], xb[0]], [v[last], ZERO_Y], ":", color=COLS[n], lw=1.2)
                ax.plot(xb, np.full(len(xb), ZERO_Y), "o", mfc="white", mec=COLS[n], mew=1.5, ms=5)
        ax.axhline(0.1, color=MUTED, ls="--", lw=0.8)
        style(ax)
        ax.set_xticks(x); ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xticklabels(xlabels, fontsize=8 if "→" in xlabels[0] else 10)
        ax.set_xlim(x[0] / 1.25, x[-1] * 1.25)
        yt = [ZERO_Y] + [10.0 ** k for k in range(-6, 3)]
        ax.set_yticks(yt); ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.set_yticklabels(["0"] + [np.format_float_positional(t, trim="-") for t in yt[1:]])
        ax.set_ylim(ZERO_Y / 2.2, 300)
        ax.text(1.01, 0.1, "0.1 %", color=MUTED, fontsize=8, ha="left", va="center",
                transform=ax.get_yaxis_transform())
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontsize=10)
        mode_legend(ax)
        fig.tight_layout(); fig.savefig(HERE / fname, dpi=160, bbox_inches="tight")

    # ---- Figures: zeroed versions (stage2_fig_mesh_EB_zeroed.png, stage2_fig_refinement_EB_zeroed.png) ----
    xm = np.array(MESHES)
    plot_zeroed(xm, [np.abs([mesh_err_EB[10.0][i][n] for i in range(len(MESHES))]) for n in range(6)],
                [str(m_) for m_ in MESHES], "Number of Elements", "Error vs Analytical (%)",
                "Euler–Bernoulli FEM vs Exact Solution\n(Same for 4, 10, 20 m)",
                "stage2_fig_mesh_EB_zeroed.png")
    plot_zeroed(xm[1:], [np.array([mesh_chg_EB[10.0][i][n] for i in range(1, len(MESHES))]) for n in range(6)],
                [f"{a}→{b}" for a, b in zip(MESHES[:-1], MESHES[1:])],
                "Mesh Refinement (Number of Elements, Before → After)", "Change in Frequency (%)",
                "Euler–Bernoulli FEM: Change in Frequency per Mesh Doubling\n(Same for 4, 10, 20 m)",
                "stage2_fig_refinement_EB_zeroed.png")
    print("\nSaved CSV tables and figures.")


# ##########################################################################################
# 5. ANSYS MESH CONVERGENCE  (was plot_ansys_mesh_convergence.py)
# ##########################################################################################
# ANSYS modal results for the same three LNG-filled cantilevers, typed in from
# ANSYS_mesh_convergence.xlsx so no spreadsheet is needed.
#
# Duplicate modes: a round pipe bends equally in two perpendicular planes, so ANSYS reports
# each bending frequency twice; only one of each pair is kept. Modes that appear once are
# torsional or axial; they are kept separately and not counted as bending modes.
#
# Outputs, for each length (10 m, 4 m, 20 m):
#   ansys_<L>_frequencies.png   bending modes 1-6 vs element count
#   ansys_<L>_change_Hz.png     bending modes 1-6: |f(N) - f(N/2)| in Hz
#   ansys_<L>_change_pct.png    the same change as a percentage of f(N)
# plus ansys_bending_modes.csv (de-duplicated frequencies) and ansys_mesh_changes.csv.
# Steps with no change visible at the exported precision (5 significant figures) are plotted
# as hollow markers at that precision limit (an upper bound).

ANSYS_ELEMENTS = [2, 4, 8, 16, 32, 64, 128, 256]

# Rows = ANSYS modes in solver order (duplicates included, as exported); columns = ANSYS_ELEMENTS
ANSYS_FREQ_HZ = {
    "10 m": [
        [5.4425, 5.4948, 5.5052, 5.5076, 5.5082, 5.5084, 5.5084, 5.5084],
        [5.4425, 5.4948, 5.5052, 5.5076, 5.5082, 5.5084, 5.5084, 5.5084],
        [34.891, 33.628, 32.86, 32.645, 32.59, 32.576, 32.573, 32.572],
        [34.891, 33.628, 32.86, 32.645, 32.59, 32.576, 32.573, 32.572],
        [80.375, 80.356, 80.355, 80.355, 80.355, 80.355, 80.355, 80.355],
        [89.361, 88.032, 86.637, 84.837, 84.369, 84.251, 84.222, 84.214],
        [120.94, 92.336, 86.637, 84.837, 84.369, 84.251, 84.222, 84.214],
        [120.94, 92.336, 87.7, 87.616, 87.596, 87.59, 87.589, 87.589],
        [235.59, 172.17, 158.73, 152.21, 150.49, 150.05, 149.94, 149.91],
        [235.59, 172.17, 158.73, 152.21, 150.49, 150.05, 149.94, 149.91],
        [245.05, 241.37, 241.08, 229.97, 225.65, 224.55, 224.28, 224.21],
        [303.78, 274.67, 245.79, 229.97, 225.65, 224.55, 224.28, 224.21],
        [449.07, 329.36, 245.79, 241.07, 241.06, 241.06, 241.06, 241.06],
        [575.84, 329.36, 265.76, 263.51, 262.95, 262.81, 262.78, 262.77],
        [654.72, 405.24, 344.68, 314.9, 306.17, 303.95, 303.4, 303.26],
        [721.09, 413.59, 344.68, 314.9, 306.17, 303.95, 303.4, 303.26],
        [1274.1, 413.59, 402.02, 401.79, 389.93, 386.01, 385.02, 384.78],
        [1274.1, 490.01, 449.67, 405.3, 389.93, 386.01, 385.02, 384.78],
        [1280.5, 500.99, 449.67, 405.3, 401.78, 401.77, 401.77, 401.77],
        [1280.5, 500.99, 451.76, 441.41, 438.81, 438.16, 438, 437.96],
    ],
    "4 m": [
        [32.394, 32.644, 32.691, 32.702, 32.705, 32.705, 32.706, 32.706],
        [32.394, 32.644, 32.691, 32.702, 32.705, 32.705, 32.706, 32.706],
        [170.9, 162.68, 158.99, 158, 157.74, 157.68, 157.66, 157.66],
        [170.9, 162.68, 158.99, 158, 157.74, 157.68, 157.66, 157.66],
        [200.94, 200.89, 200.89, 200.89, 200.89, 200.89, 200.89, 200.89],
        [223.4, 220.08, 219.25, 219.04, 218.99, 218.98, 218.98, 218.97],
        [487.09, 384.37, 360.28, 353.17, 351.36, 350.9, 350.79, 350.76],
        [487.09, 384.37, 360.28, 353.17, 351.36, 350.9, 350.79, 350.76],
        [612.61, 603.41, 583.7, 560.9, 554.99, 553.51, 553.14, 553.04],
        [659.2, 634.61, 583.7, 560.9, 554.99, 553.51, 553.14, 553.04],
        [659.2, 634.61, 602.71, 602.66, 602.66, 602.66, 602.66, 602.66],
        [759.46, 686.67, 664.41, 658.79, 657.39, 657.04, 656.95, 656.93],
        [1122.7, 1013.1, 826.97, 776.88, 763.36, 759.95, 759.1, 758.88],
        [1344.6, 1075.9, 826.97, 776.88, 763.36, 759.95, 759.1, 758.88],
        [1344.6, 1075.9, 1005.1, 991.79, 967.05, 960.75, 959.17, 958.77],
        [1439.6, 1212.4, 1074.3, 991.79, 967.05, 960.75, 959.17, 958.77],
        [1624.4, 1212.4, 1074.3, 1004.5, 1004.4, 1004.4, 1004.4, 1004.4],
        [1624.4, 1225, 1129.4, 1103.5, 1097, 1095.4, 1095, 1094.9],
        [1636.8, 1329.9, 1295.1, 1203.1, 1166.1, 1156.3, 1153.8, 1153.2],
        [1802.7, 1329.9, 1295.1, 1203.1, 1166.1, 1156.3, 1153.8, 1153.2],
    ],
    "20 m": [
        [1.3707, 1.3843, 1.387, 1.3877, 1.3878, 1.3879, 1.3879, 1.3879],
        [1.3707, 1.3843, 1.387, 1.3877, 1.3878, 1.3879, 1.3879, 1.3879],
        [9.1528, 8.8452, 8.6419, 8.5844, 8.5697, 8.5659, 8.565, 8.5648],
        [9.1528, 8.8452, 8.6419, 8.5844, 8.5697, 8.5659, 8.565, 8.5648],
        [34.213, 25.733, 24.122, 23.6, 23.463, 23.429, 23.42, 23.418],
        [34.213, 25.733, 24.122, 23.6, 23.463, 23.429, 23.42, 23.418],
        [40.188, 40.178, 40.177, 40.177, 40.177, 40.177, 40.177, 40.177],
        [44.68, 44.016, 43.85, 43.808, 43.798, 43.795, 43.795, 43.794],
        [92.095, 51.165, 47.179, 45.142, 44.598, 44.46, 44.425, 44.416],
        [92.095, 51.165, 47.179, 45.142, 44.598, 44.46, 44.425, 44.416],
        [122.52, 105.85, 77.83, 72.561, 71.097, 70.725, 70.632, 70.608],
        [151.89, 105.85, 77.83, 72.561, 71.097, 70.725, 70.632, 70.608],
        [224.54, 120.68, 115.65, 105.14, 102, 101.2, 101, 100.95],
        [287.92, 137.33, 115.65, 105.14, 102, 101.2, 101, 100.95],
        [327.36, 149.83, 120.54, 120.53, 120.53, 120.53, 120.53, 120.53],
        [360.55, 149.83, 132.88, 131.76, 131.48, 131.41, 131.39, 131.38],
        [830.05, 202.62, 158.98, 142.25, 136.47, 134.98, 134.6, 134.51],
        [830.05, 211.53, 158.98, 142.25, 136.47, 134.98, 134.6, 134.51],
        [1166.8, 211.53, 201.01, 183.42, 173.79, 171.29, 170.66, 170.51],
        [1166.8, 245, 201.8, 183.42, 173.79, 171.29, 170.66, 170.51],
    ],
}

ANSYS_N_BENDING = 6
ANSYS_RTOL = 1e-6                 # two frequencies this close count as the same (duplicate) mode
ANSYS_SIG_FIGS = 5                # ANSYS values were exported to 5 significant figures
ANSYS_COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
ANSYS_MUTED, ANSYS_GRID = "#8a8984", "#e6e5e0"


def ansys_split_modes(column):
    """Remove duplicates from one mesh's frequency list.
    Returns (bending, other): bending = frequencies that occur as pairs (one kept),
    other = frequencies that occur once (torsional/axial). Both ascending."""
    vals = np.sort(np.asarray(column, float))
    bending, other, i = [], [], 0
    while i < len(vals):
        if i + 1 < len(vals) and np.isclose(vals[i], vals[i + 1], rtol=ANSYS_RTOL):
            bending.append(vals[i]); i += 2
        else:
            other.append(vals[i]); i += 1
    return np.array(bending), np.array(other)


def ansys_pad(rows, n):
    """Stack per-mesh lists into an (n modes x meshes) array, NaN where a mesh has fewer modes."""
    out = np.full((n, len(rows)), np.nan)
    for j, r in enumerate(rows):
        k = min(n, len(r)); out[:k, j] = r[:k]
    return out


def ansys_style(ax):
    """Common look: no top/right spines, light grid behind the data."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color=ANSYS_GRID, lw=0.6); ax.set_axisbelow(True)


def ansys_resolution(f):
    """Smallest change visible in a value exported to ANSYS_SIG_FIGS significant figures."""
    return 10.0 ** (np.floor(np.log10(np.abs(f))) - (ANSYS_SIG_FIGS - 1))


def ansys_plain(v, _pos=None):
    """Log-axis tick label as a plain number (100, 1, 0.001) instead of 10^n."""
    if v <= 0:
        return ""
    d = max(0, -int(np.floor(np.log10(v) + 1e-9)))
    return f"{v:.{d}f}"


def ansys_save(fig, ax, name, ncol=6):
    """Style, add legend below the axes, save at 400 dpi, then close (9 figures; keeps windows down)."""
    import matplotlib.pyplot as plt
    ansys_style(ax)
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=ncol, fontsize=9,
                    frameon=True, fancybox=False, framealpha=1.0, facecolor="white",
                    edgecolor="black")
    leg.get_frame().set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(HERE / name, dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}")


def ansys_analyse(label, n_el, raw):
    """One length: de-duplicate modes, print tables, save the three figures.
    Returns (bending freqs, change Hz, change %, mask of 'below precision' steps)."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    # ---- split each mesh's modes into bending (pairs) and torsion/axial (singles) ----
    raw = np.asarray(raw, float)
    split = [ansys_split_modes(raw[:, j]) for j in range(raw.shape[1])]
    B = ansys_pad([s[0] for s in split], ANSYS_N_BENDING)          # bending 1-6, rows = mode

    print(f"\n==================== ANSYS {label} ====================")
    print(f"Distinct modes at finest mesh ({int(n_el[-1])} elements):")
    print("  bending      :", np.round(split[-1][0], 3))
    print("  torsion/axial:", np.round(split[-1][1], 3))
    print("\nBending modes (Hz), duplicates removed")
    print(f"{'N':>6}" + "".join(f"{'B'+str(k+1):>11}" for k in range(ANSYS_N_BENDING)))
    for j, ne in enumerate(n_el):
        print(f"{int(ne):>6}" + "".join(f"{B[k, j]:>11.4f}" for k in range(ANSYS_N_BENDING)))

    # ---- change per mesh doubling ----
    dF = np.abs(np.diff(B, axis=1))
    dP = dF / B[:, 1:] * 100
    # where no change is visible at the exported precision, use the precision
    # limit as an upper bound so every step can still be plotted and reported
    zero = dF == 0
    resF = ansys_resolution(B[:, 1:])
    dF_plot = np.where(zero, resF, dF)
    dP_plot = np.where(zero, resF / B[:, 1:] * 100, dP)
    n_fine = n_el[1:]
    print("\nChange from previous mesh |f(N) - f(N/2)|: Hz (%)")
    print(f"{'N':>6}" + "".join(f"{'B'+str(k+1):>20}" for k in range(ANSYS_N_BENDING)))
    for j, ne in enumerate(n_fine):
        print(f"{int(ne):>6}" + "".join(
            (f"  <{dF_plot[k, j]:>8.4f} (<{dP_plot[k, j]:5.3f}%)" if zero[k, j]
             else f"{dF[k, j]:>11.4f} ({dP[k, j]:6.3f}%)") for k in range(ANSYS_N_BENDING)))
    print("  '<' = no change visible at the exported precision; value shown is the upper bound")

    tag = label.replace(" ", "")
    Lm = label.replace(" m", "")
    xpos = np.arange(len(n_fine))
    xlab = [f"{int(a)}→{int(b)}" for a, b in zip(n_el[:-1], n_el[1:])]

    # ---- figure: natural frequencies vs mesh ----
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for k in range(ANSYS_N_BENDING):
        ax.plot(n_el, B[k], "-o", color=ANSYS_COLS[k], lw=2, ms=4, label=f"Mode {k+1}")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(n_el); ax.set_xticklabels([str(int(v)) for v in n_el])
    ax.yaxis.set_major_formatter(FuncFormatter(ansys_plain))
    ax.set_xlabel("Number of Elements"); ax.set_ylabel("Natural Frequency (Hz)")
    ax.set_title(f"ANSYS, $L$ = {Lm} m: Natural Frequency vs Mesh\n"
                 f"(Bending Modes)", loc="left", fontsize=11)
    ansys_save(fig, ax, f"ansys_{tag}_frequencies.png", ncol=6)

    # ---- figures: change per mesh doubling, Hz and % ----
    for data, unit, ylab, fname in (
            (dF_plot, "Hz", "Change in Frequency (Hz)", f"ansys_{tag}_change_Hz.png"),
            (dP_plot, "%", "Change in Frequency (%)", f"ansys_{tag}_change_pct.png")):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        for k in range(ANSYS_N_BENDING):
            ax.semilogy(xpos, data[k], "-o", color=ANSYS_COLS[k], lw=2, ms=4, label=f"Mode {k+1}")
            ax.semilogy(xpos[zero[k]], data[k][zero[k]], "o", mfc="white", mec=ANSYS_COLS[k],
                        mew=1.5, ms=6, zorder=5)
        if zero.any():
            ax.plot([], [], "o", mfc="white", mec=ANSYS_MUTED, mew=1.5, ms=6,
                    label="Below precision (upper bound)")
        if unit == "%":
            ax.axhline(0.1, color="black", ls=":", lw=1.4, zorder=1,
                       label="Convergence criterion (0.1 %)")
            ax.annotate("0.1 %", xy=(1.0, 0.1), xycoords=("axes fraction", "data"),
                        xytext=(6, 0), textcoords="offset points",
                        va="center", ha="left", color="black", fontsize=9)
        ax.set_xticks(xpos); ax.set_xticklabels(xlab, fontsize=8)
        ax.yaxis.set_major_formatter(FuncFormatter(ansys_plain))
        ax.set_xlabel("Mesh Refinement (Number of Elements, Before → After)")
        ax.set_ylabel(ylab)
        ax.set_title(f"ANSYS, $L$ = {Lm} m: Change in Frequency per Mesh Doubling\n"
                     f"(Bending Modes, {unit})", loc="left", fontsize=11)
        ansys_save(fig, ax, fname, ncol=4)
    return B, dF_plot, dP_plot, zero


def run_ansys():
    """Section 5 driver: run ansys_analyse for each length (prints + 3 figures each), then
    collect every length's results into two long-format CSV tables:
      ansys_bending_modes.csv - one row per (length, mesh): bending modes 1-6 in Hz
      ansys_mesh_changes.csv  - one row per (length, mesh step, mode): change in Hz and %,
                                flagged 'measured' or 'upper bound' (below export precision)"""
    n_el = np.array(ANSYS_ELEMENTS, float)
    freq_rows, change_rows = [], []
    for label, raw in ANSYS_FREQ_HZ.items():          # lengths in data order: 10 m, 4 m, 20 m
        B, dF, dP, zero = ansys_analyse(label, n_el, raw)
        # rows for ansys_bending_modes.csv
        for j, ne in enumerate(n_el):
            freq_rows.append([label, int(ne)] + [f"{B[k, j]:.6g}" for k in range(ANSYS_N_BENDING)])
        # rows for ansys_mesh_changes.csv (one per doubling step and mode)
        for j in range(dF.shape[1]):
            step = f"{int(n_el[j])}->{int(n_el[j + 1])}"
            for k in range(ANSYS_N_BENDING):
                change_rows.append([label, step, k + 1, f"{dF[k, j]:.6g}", f"{dP[k, j]:.6g}",
                                    "upper bound (below export precision)" if zero[k, j] else "measured"])
    # ---- write both tables ----
    with open(HERE / "ansys_bending_modes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["length", "elements"] + [f"bending_{k+1}_Hz" for k in range(ANSYS_N_BENDING)])
        w.writerows(freq_rows)
    with open(HERE / "ansys_mesh_changes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["length", "mesh_step", "mode", "change_Hz", "change_pct", "note"])
        w.writerows(change_rows)
    print("\nSaved ansys_bending_modes.csv and ansys_mesh_changes.csv")


# ##########################################################################################
# 6. COMPARISON: ANSYS vs EULER-BERNOULLI THEORY vs OWN FEM
# ##########################################################################################
# For each length, bending modes 1-6 from:
#   EB theory - exact Euler-Bernoulli cantilever (f_EB, section 3)
#   FEM       - own Euler-Bernoulli beam FEM at the converged mesh (f_FEM, CMP_NE elements)
#   ANSYS     - finest ANSYS mesh (256 elements), duplicates removed
# % differences are all relative to EB theory (ANSYS vs FEM also given in the CSV).
# Reading the result: FEM vs EB is ~0 (same beam theory, converged mesh). ANSYS sits below EB,
# more so for higher modes and shorter pipes (smaller L/D). That pattern is what shear
# deformation and rotary inertia (Timoshenko beam / pipe element) produce, which EB theory
# ignores - check the ANSYS element formulation before quoting this as the cause.
# Outputs:
#   stage2_table_ansys_vs_EB_FEM.csv     all frequencies + % differences
#   stage2_fig_ansys_vs_EB_FEM_freq.png  frequencies per mode, one panel per length
#   stage2_fig_ansys_vs_EB_FEM_pct.png   % difference from EB theory per mode

CMP_NE = 64                       # own-FEM mesh used for comparison (same as NCONV in section 4)


def run_comparison():
    """Section 6 driver: ANSYS vs Euler-Bernoulli theory vs own FEM.
    Order: collect frequencies + % differences (printed) -> CSV -> frequency figure -> % figure.
    Uses only the finest ANSYS mesh, so it does not depend on section 5 having run."""
    import matplotlib.pyplot as plt   # imported here so main() can pick the backend first

    lengths = [4.0, 10.0, 20.0]
    res, rows = {}, []                # res: L -> (f_EB, f_FEM, f_ANSYS, FEM-EB %, ANSYS-EB %)

    # ---- collect the three sets of frequencies and the % differences ----
    # dFE = FEM vs EB theory      (discretisation check, should be ~0)
    # dAN = ANSYS vs EB theory    (main result: effect of the ANSYS model's extra physics)
    # dAF = ANSYS vs own FEM      (CSV only)
    print("=" * 86)
    print(f"ANSYS (256 el.) vs EB theory vs own FEM ({CMP_NE} el.) - bending modes, Hz")
    print("=" * 86)
    for L in lengths:
        fa = f_EB(L)
        fe = f_FEM(L, CMP_NE)
        raw = np.asarray(ANSYS_FREQ_HZ[f"{L:.0f} m"], float)
        fn = ansys_split_modes(raw[:, -1])[0][:6]               # finest mesh, bending only
        dFE = (fe - fa) / fa * 100
        dAN = (fn - fa) / fa * 100
        dAF = (fn - fe) / fe * 100
        res[L] = (fa, fe, fn, dFE, dAN)
        print(f"\nL = {L:.0f} m   (L/D = {L/Do:.1f})")
        print(f"{'mode':>4} {'EB':>10} {'FEM':>10} {'ANSYS':>10} {'FEM-EB %':>10} {'ANSYS-EB %':>11}")
        for n in range(6):
            # round(...) + 0.0 turns round-off "-0.000" into "0.000" (display only; CSV keeps raw value)
            print(f"{n+1:>4} {fa[n]:>10.3f} {fe[n]:>10.3f} {fn[n]:>10.3f}"
                  f" {round(dFE[n], 3) + 0.0:>10.3f} {dAN[n]:>11.2f}")
            rows.append([L, n + 1, fa[n], fe[n], fn[n], dFE[n], dAN[n], dAF[n]])

    # ---- table: stage2_table_ansys_vs_EB_FEM.csv (one row per length and mode) ----
    with open(HERE / "stage2_table_ansys_vs_EB_FEM.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["L_m", "mode", "f_EB_analytical_Hz", "f_FEM_Hz", "f_ANSYS_Hz",
                    "FEM_vs_EB_pct", "ANSYS_vs_EB_pct", "ANSYS_vs_FEM_pct"])
        for row in rows:
            w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row])

    modes = np.arange(1, 7)
    LCOL = {4.0: "#2a78d6", 10.0: "#eb6834", 20.0: "#1baf7a"}   # one colour per length

    # ---- figure: frequencies, one panel per length (log y; 3 methods overlaid) ----
    # stage2_fig_ansys_vs_EB_FEM_freq.png. Theory = line, FEM = crosses (sit on the line),
    # ANSYS = hollow circles (drop below the line as mode number rises / length falls).
    # Each panel has its own y range because frequencies scale as 1/L^2. One shared legend below.
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.9))
    for ax, L in zip(axes, lengths):
        fa, fe, fn, _, _ = res[L]
        ax.semilogy(modes, fa, "-", color="#444441", lw=1.5, label="Euler–Bernoulli theory")
        ax.semilogy(modes, fe, "x", color="#2a78d6", ms=8, mew=2, label=f"FEM ({CMP_NE} el.)")
        ax.semilogy(modes, fn, "o", mfc="white", mec="#eb6834", ms=7, mew=1.8, label="ANSYS (256 el.)")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
        ax.set_xticks(modes); ax.set_xlabel("Mode")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(ansys_plain))
        ax.set_title(f"$L$ = {L:.0f} m", loc="left", fontsize=10)
    axes[0].set_ylabel("Natural Frequency (Hz)")
    fig.suptitle("Bending Frequencies: ANSYS vs Euler–Bernoulli Theory vs FEM", x=0.01, ha="left", fontsize=11)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=9, frameon=True, fancybox=False,
               edgecolor="#8a8984", framealpha=1, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(HERE / "stage2_fig_ansys_vs_EB_FEM_freq.png", dpi=160, bbox_inches="tight")

    # ---- figure: % difference from EB theory (ANSYS solid, FEM dashed ~0) ----
    # stage2_fig_ansys_vs_EB_FEM_pct.png. Linear y axis so the 0 % FEM lines are visible.
    # Colour = length (LCOL), line style = method.
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for L in lengths:
        _, _, _, dFE, dAN = res[L]
        ax.plot(modes, dAN, "-o", color=LCOL[L], lw=2, ms=5, label=f"ANSYS, {L:.0f} m")
        ax.plot(modes, dFE, "--x", color=LCOL[L], lw=1.2, ms=6, label=f"FEM, {L:.0f} m")
    ax.axhline(0, color="#8a8984", lw=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    ax.set_xticks(modes); ax.set_xlabel("Mode")
    ax.set_ylabel("Difference from Euler–Bernoulli Theory (%)")
    ax.set_title("ANSYS and FEM vs Euler–Bernoulli Theory", loc="left", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=8, frameon=True,
              fancybox=False, edgecolor="#8a8984", framealpha=1)
    fig.tight_layout()
    fig.savefig(HERE / "stage2_fig_ansys_vs_EB_FEM_pct.png", dpi=160, bbox_inches="tight")
    print("\nSaved stage2_table_ansys_vs_EB_FEM.csv and comparison figures.")


# ##########################################################################################
# MAIN: run everything in order
# ##########################################################################################
# main() picks the matplotlib backend (windows or --no-show), runs sections 1-6 with a banner
# for each, then checks every file in EXPECTED_OUTPUTS exists. Add any new output file here.
EXPECTED_OUTPUTS = [
    "stage2_EB_mode_shapes.png",
    "stage2_table_frequencies.csv",
    "stage2_table_comparison.csv",
    "stage2_table_mesh.csv",
    "stage2_fig_mesh_EB.png",
    "stage2_fig_refinement_EB.png",
    "stage2_fig_mesh_EB_zeroed.png",
    "stage2_fig_refinement_EB_zeroed.png",
    *[f"ansys_{L}_{kind}.png" for L in ("10m", "4m", "20m")
      for kind in ("frequencies", "change_Hz", "change_pct")],
    "ansys_bending_modes.csv",
    "ansys_mesh_changes.csv",
    "stage2_table_ansys_vs_EB_FEM.csv",
    "stage2_fig_ansys_vs_EB_FEM_freq.png",
    "stage2_fig_ansys_vs_EB_FEM_pct.png",
]


def banner(title):
    """Print a section heading between steps; flush so it appears before the step's output."""
    print("\n" + "=" * 78)
    print(f" {title}")
    print("=" * 78, flush=True)


def main():
    """Entry point: parse options, choose plot backend, run sections 1-6, report outputs."""
    # ---- command-line options ----
    parser = argparse.ArgumentParser(description="Run all Stage 2 calculations.")
    parser.add_argument("--no-show", action="store_true",
                        help="save outputs only, do not open figure windows")
    args = parser.parse_args()

    # ---- plot backend: must be chosen before pyplot is first imported (hence the local imports
    #      inside each section's function) ----
    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")              # headless: save figures, open no windows
    import matplotlib.pyplot as plt

    # ---- run every section in order: (banner title, function) ----
    steps = [
        ("1/6  Shared parameters", print_params),
        ("2/6  Euler-Bernoulli hand calculation", run_euler_bernoulli),
        ("3/6  Beam FEM checks (flow, mesh)", run_cantilever_checks),
        ("4/6  Comparison + mesh-study tables and figures", run_results),
        ("5/6  ANSYS mesh convergence", run_ansys),
        ("6/6  ANSYS vs EB theory vs FEM", run_comparison),
    ]
    for title, func in steps:
        banner(title)
        func()

    # ---- output check: confirm every expected file was written ----
    banner("Stage 2 complete - outputs in " + str(HERE))
    for f in EXPECTED_OUTPUTS:
        status = "OK     " if (HERE / f).exists() else "MISSING"
        print(f"  [{status}] {f}")
    sys.stdout.flush()

    if plt.get_backend().lower() != "agg":     # interactive run: show all figures at the end
        plt.show()


if __name__ == "__main__":
    main()
