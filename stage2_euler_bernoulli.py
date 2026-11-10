"""
Stage 2 - Euler-Bernoulli hand calculation for a liquid-filled cantilever
=========================================================================
Follows the hand-calc steps one-to-one:
  Step 4  roots of 1 + cos(bL) cosh(bL) = 0      (bisection, + asymptotic check)
  Step 6  section properties, mass per length, sqrt(EI/m)
  Step 7  worked example (L = 10 m, mode 1)
  Step 8  frequency table  f_n = C_L * (beta_n L)^2,  C_L = sqrt(EI/m) / (2 pi L^2)
  Step 9  mode shapes W_n(x) and node count
  Step 10 checks (1/L^2 scaling, mode ratios, fluid-mass effect)

References: Rao, Mechanical Vibrations (beam vibration chapter);
            Blevins, Formulas for Natural Frequency and Mode Shape (cantilever table).
Requires: numpy, matplotlib
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# INPUTS: shared group parameters (params.py, same folder)
# --------------------------------------------------------------------------
from params import D_OUT, D_IN, E_STEEL, RHO_STEEL, RHO_LNG
Do, Di, E, rho_s, rho_f = D_OUT, D_IN, E_STEEL, RHO_STEEL, RHO_LNG
LENGTHS = [4.0, 10.0, 20.0]   # m
N_MODES = 6

# --------------------------------------------------------------------------
# STEP 4 - roots of the characteristic equation  g(x) = 1 + cos(x) cosh(x) = 0
# --------------------------------------------------------------------------
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

# --------------------------------------------------------------------------
# STEP 6 - section properties
# --------------------------------------------------------------------------
I   = math.pi * (Do**4 - Di**4) / 64          # m^4
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

# --------------------------------------------------------------------------
# STEP 5/7/8 - frequencies
# --------------------------------------------------------------------------
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

# --------------------------------------------------------------------------
# STEP 9 - mode shapes
# --------------------------------------------------------------------------
def sigma(bL):
    return (math.cosh(bL) + math.cos(bL)) / (math.sinh(bL) + math.sin(bL))

def mode_shape(n, xi):
    """Normalised mode shape W_n at xi = x/L (tip value scaled to +1)."""
    b = BL[n - 1]
    s = sigma(b)
    z = b * xi
    W = np.cosh(z) - np.cos(z) - s * (np.sinh(z) - np.sin(z))
    return W / W[-1]

xi = np.linspace(0, 1, 2001)
print("\nSTEP 9 - Mode shape constants and nodes")
for n in range(1, N_MODES + 1):
    W = mode_shape(n, xi)
    nodes = np.sum(np.diff(np.sign(W[5:])) != 0)     # skip the clamped end itself
    print(f"  mode {n}: sigma = {sigma(BL[n-1]):.4f}, interior nodes = {nodes} (expected {n-1})")

fig, ax = plt.subplots(figsize=(8, 4.5))
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
for n in range(1, N_MODES + 1):
    ax.plot(xi, mode_shape(n, xi), lw=2, color=colors[n-1], label=f"mode {n}")
ax.axhline(0, color="#8a8984", lw=0.8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_xlabel("x / L  (0 = fixed end, 1 = free end)")
ax.set_ylabel("W_n (tip normalised to 1)")
ax.set_title("Euler-Bernoulli cantilever mode shapes (independent of L)", loc="left")
ax.legend(frameon=False, ncol=3, fontsize=8)
fig.tight_layout()
fig.savefig("stage2_EB_mode_shapes.png", dpi=150)

# --------------------------------------------------------------------------
# STEP 10 - checks
# --------------------------------------------------------------------------
print("\nSTEP 10 - Checks")
f10, f20, f4 = table[10][0], table[20][0], table[4][0]
print(f"  1/L^2 scaling: f1(10)/4 = {f10/4:.3f} vs f1(20) = {f20:.3f};"
      f"  f1(10)*6.25 = {f10*6.25:.2f} vs f1(4) = {f4:.2f}")
print("  Mode ratios f_n/f_1 :", " : ".join(f"{(b/BL[0])**2:.3f}" for b in BL))
f_steel = freqs(10, mass=m_s)[0]
print(f"  Steel only, L = 10 m: f1 = {f_steel:.3f} Hz;  ratio to filled = {f_steel/f10:.3f}"
      f"  (sqrt(m/m_s) = {math.sqrt(m/m_s):.3f})")
print("\nSaved stage2_EB_mode_shapes.png")
