"""
Stage 2 - Liquid-filled cantilever: analytical vs FEM natural bending frequencies.

Analytical : Euler-Bernoulli cantilever, f_n = (beta_n L)^2 / (2 pi L^2) * sqrt(EI/m)
             beta_n L roots of 1 + cos(bL) cosh(bL) = 0  (Blevins; Rao)
FEM        : 2-node beam elements, consistent mass.
             - Euler-Bernoulli (Hermite cubic) element
             - Timoshenko element (Friedman & Kosmatka 1993, Comput. Struct. 47:473-481)
Fluid      : internal LNG treated as mass moving with the pipe (added to m, not to EI).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh

# ---------------- inputs: shared group parameters (params.py, same folder) ----------------
from params import D_OUT, D_IN, RHO_STEEL, RHO_LNG, E_STEEL, U0, A_STEEL, A_BORE, I_SECTION, M_STEEL, M_LNG, M_TOTAL

Do, Di = D_OUT, D_IN
rho_s, rho_f = RHO_STEEL, RHO_LNG
E   = E_STEEL
U   = U0                 # m/s, flow velocity (for flow-effect check)

# Not yet in params.py - local assumption (propose adding to params.py):
nu  = 0.29               # ASSUMPTION: Poisson's ratio, 304 SS (check source)

As = A_STEEL                              # steel area
Af = A_BORE                               # fluid area
I  = I_SECTION                            # second moment of area (steel)
m_s, m_f = M_STEEL, M_LNG
m  = M_TOTAL                              # mass per unit length (pipe + contents)
G  = E / (2 * (1 + nu))
kappa = 2 * (1 + nu) / (4 + 3 * nu)       # Cowper (1966) shear coefficient, thin circular tube
EI = E * I

# ---------------- analytical (Euler-Bernoulli) ----------------
def bl_roots(n=6):
    g = lambda x: 1 + np.cos(x) * np.cosh(x)
    return np.array([brentq(g, (k - 0.5) * np.pi - 0.6, (k - 0.5) * np.pi + 0.6) for k in range(1, n + 1)])

BL = bl_roots()

def f_EB(L):
    return BL**2 / (2 * np.pi * L**2) * np.sqrt(EI / m)

# ---------------- FEM ----------------
def element(Le, timoshenko):
    phi = 12 * EI / (kappa * G * As * Le**2) if timoshenko else 0.0
    L = Le
    k = EI / ((1 + phi) * L**3) * np.array([
        [12,      6*L,            -12,     6*L],
        [6*L,     (4+phi)*L**2,   -6*L,    (2-phi)*L**2],
        [-12,    -6*L,             12,    -6*L],
        [6*L,     (2-phi)*L**2,   -6*L,    (4+phi)*L**2]])
    a1 = 13/35 + 7*phi/10 + phi**2/3
    a2 = (11/210 + 11*phi/120 + phi**2/24) * L
    a3 = 9/70 + 3*phi/10 + phi**2/6
    a4 = (13/420 + 3*phi/40 + phi**2/24) * L
    a5 = (1/105 + phi/60 + phi**2/120) * L**2
    a6 = (1/140 + phi/60 + phi**2/120) * L**2
    mt = m * L / (1 + phi)**2 * np.array([
        [a1,  a2,  a3, -a4],
        [a2,  a5,  a4, -a6],
        [a3,  a4,  a1, -a2],
        [-a4, -a6, -a2, a5]])
    if timoshenko:   # rotary inertia of the steel wall (ASSUMPTION: fluid rotary inertia neglected)
        b1 = 6/5
        b2 = (1/10 - phi/2) * L
        b3 = (2/15 + phi/6 + phi**2/3) * L**2
        b4 = (-1/30 - phi/6 + phi**2/6) * L**2
        mr = rho_s * I / ((1 + phi)**2 * L) * np.array([
            [b1,  b2, -b1,  b2],
            [b2,  b3, -b2,  b4],
            [-b1, -b2, b1, -b2],
            [b2,  b4, -b2,  b3]])
        mt = mt + mr
    return k, mt

def f_FEM(L, ne, timoshenko=False, nmodes=6):
    nd = 2 * (ne + 1)
    K, M = np.zeros((nd, nd)), np.zeros((nd, nd))
    ke, me = element(L / ne, timoshenko)
    for e in range(ne):
        s = slice(2*e, 2*e + 4)
        K[s, s] += ke; M[s, s] += me
    K, M = K[2:, 2:], M[2:, 2:]           # clamp node 0 (w = 0, theta = 0)
    w2 = eigh(K, M, eigvals_only=True)[:nmodes]
    f = np.full(nmodes, np.nan); f[:len(w2)] = np.sqrt(w2) / (2 * np.pi)
    return f

if __name__ == "__main__":
    print("---- Section properties ----")
    print(f"A_steel = {As:.5f} m2   A_fluid = {Af:.5f} m2")
    print(f"I       = {I:.5e} m4   EI = {EI:.4e} N m2")
    print(f"m_steel = {m_s:.2f} kg/m  m_fluid = {m_f:.2f} kg/m  m_total = {m:.2f} kg/m")
    print(f"radius of gyration r = sqrt(I/As) = {np.sqrt(I/As):.4f} m")
    print(f"G = {G/1e9:.2f} GPa, kappa = {kappa:.4f}")
    print("beta_n L =", np.round(BL, 4))

    print("\n---- Flow-velocity check (Paidoussis dimensionless velocity) ----")
    for L in (4, 10, 20):
        u = U * L * np.sqrt(m_f / EI)
        print(f"L = {L:2d} m : u* = U L sqrt(M_f/EI) = {u:.2e}")

    for L in (4, 10, 20):
        fa = f_EB(L)
        print(f"\n==== L = {L} m  (L/D = {L/Do:.1f}, L/r = {L/np.sqrt(I/As):.1f}) ====")
        print("Analytical EB (Hz):", np.round(fa, 3))
        print("Mesh study, EB FEM: % error vs analytical")
        for ne in (2, 4, 8, 16, 32, 64):
            ff = f_FEM(L, ne)
            err = (ff - fa) / fa * 100
            print(f"  ne={ne:3d}  " + "  ".join(("      n/a" if np.isnan(e) else f"{e:+9.4f}") for e in err))
        ft = f_FEM(L, 64, timoshenko=True)
        print("Timoshenko FEM, ne=64 (Hz):", np.round(ft, 3))
        print("Timoshenko vs EB (%):      ", np.round((ft - fa) / fa * 100, 2))
