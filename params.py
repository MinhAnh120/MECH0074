"""Shared parameters for the MECH0074 LNG expansion loop coursework.

Every model imports from here, so all stages use the same values:
    from params import *

Supplied values are taken from the brief. Anything marked TBC is still to be
agreed by the group; log any change and its reason in the assumptions log on Teams.
"""
import math

# ---------------------------------------------------------------- pipe geometry (brief)
D_OUT = 0.813            # outer diameter, m
D_IN = 0.79394           # inner diameter, m
T_WALL = (D_OUT - D_IN) / 2          # wall thickness, m (~9.53 mm)
BEND_RADIUS = 1.4        # radius of curvature of each of the four bends, m
N_BENDS = 4

# ---------------------------------------------------------------- material: 304 stainless (brief)
RHO_STEEL = 8000.0       # kg/m^3
E_STEEL = 213.233e9      # Young's modulus (NIST fit), Pa
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
# The brief gives the bulk sound speed. Wall elasticity lowers the speed in the pipe
# (Korteweg correction). Set USE_KORTEWEG and record the choice in the assumptions log.
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


if __name__ == "__main__":
    print(f"Wall thickness      {T_WALL*1e3:.2f} mm")
    print(f"Mass per length     {M_TOTAL:.1f} kg/m (steel {M_STEEL:.1f}, LNG {M_LNG:.1f})")
    print(f"Second moment       {I_SECTION:.4e} m^4")
    print(f"Wave speed used     {C_WAVE:.1f} m/s (Korteweg would give {C_KORTEWEG:.1f})")
    print(f"Joukowsky surge     {DP_JOUKOWSKY/1e6:.3f} MPa")
    print(f"2L/c, 4L/c          {T_RETURN:.3f} s, {T_PERIOD:.3f} s ({F_HYDRAULIC:.2f} Hz)")
    print(f"Contraction/100 m   {DL_PER_100M*1e3:.1f} mm")
    print(f"Restrained stress   {SIGMA_RESTRAINED/1e6:.1f} MPa")
