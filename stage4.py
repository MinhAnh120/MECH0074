"""
MECH0074 CW1 - Stage 4 (coupled response) toolkit  -- FULLY COMMENTED VERSION
=============================================================================
Big picture (read this first):

  valve closes -> pressure wave travels along the pipe (Part 1: MOC solver)
       -> wave reaches each bend of the expansion loop at a different time
       -> the pressure on each bend creates a force on the pipe (Part 4: bend forces)
       -> the pipe structure (a frame of beam elements, Part 2) responds to those forces
       -> we solve the vibration (Part 5), then get stress and support forces (Part 6)
       -> we repeat for a baseline and modified designs and compare (Part 10)

THIS IS A STARTING POINT. Replace the pieces your group already owns:
  * pressure histories  -> pass your Stage-1 arrays to `Hyd(...)` instead of calling moc()
  * mass/stiffness/modes -> pass your Stage-2/3 FE matrices to `System(...)`
Everything that is an ASSUMPTION is marked  # ASSUMPTION  - justify or change it.

Units: SI throughout (m, kg, s, Pa, N).  Run:  python stage4_commented.py
"""
# dataclass lets us define a neat "parameter bag"; replace() makes a modified copy of it
from dataclasses import dataclass, replace
import numpy as np                                   # arrays and maths
from scipy.linalg import eigh, lu_factor, lu_solve   # eigenvalues (modes) and fast linear solves
import matplotlib
matplotlib.use("Agg")                                # draw plots to files, no window needed
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------
# 0. Parameters (teaching case from the brief + assumptions)
# ----------------------------------------------------------------------------
@dataclass
class Par:
    """All inputs live here, so a 'modification' is just a copy with one value changed."""
    # --- values given in the coursework brief ---
    Do: float = 0.813              # pipe outer diameter (m)
    Di: float = 0.79394            # pipe inner diameter (m)
    rho_f: float = 455.0           # LNG density (kg/m^3)
    c: float = 1100.0              # sound speed in the LNG (m/s); see korteweg() for pipe elasticity
    u0: float = 1.5                # initial flow velocity (m/s)
    p_res: float = 0.7e6           # absolute pressure at the reservoir = initial pressure (Pa)
    p_ext: float = 101325.0        # external (atmospheric) absolute pressure (Pa)
    rho_s: float = 8000.0          # steel density (kg/m^3)
    E: float = 213.233e9           # Young's modulus of 304 stainless at cryogenic temperature (Pa)
    alpha: float = 1.4e-5          # thermal expansivity (1/K). CHECK the exponent in the PDF (garbled)
    T0: float = 293.0              # starting temperature (K)
    T1: float = 111.0              # final (LNG) temperature (K)
    R: float = 1.4                 # bend centre-line radius (m)
    t_close: float = 0.02          # valve closes linearly over this time (s)
    # --- loop layout. ASSUMPTION: the brief gives only the bends and radius ---
    H: float = 10.0                # length of the two vertical legs of the Pi-loop (m)
    W: float = 10.0                # length of the top (horizontal) leg of the loop (m)
    La: float = 30.0               # straight pipe between reservoir and loop (m)
    Lb: float = 30.0               # straight pipe between loop and valve (m)
    guides: tuple = ()             # lateral restraints: ((s_position_m, dof), ...) dof 0=x, 1=y
    # --- modelling choices ---
    f_darcy: float = 0.0           # ASSUMPTION: frictionless pipe (set ~0.01 to see waves decay)
    zeta: float = 0.01             # ASSUMPTION: 1 % damping ratio at modes 1 and 6 (Rayleigh damping)
    h_max: float = 1.0             # structural element size on straight pipes (m)
    n_arc: int = 8                 # chord (beam) elements used to approximate each 90 deg bend; keep even
    flex_bends: bool = False       # True: reduce bend bending stiffness by the Karman flexibility factor
    sif: bool = False              # True: multiply bend bending stress by the stress intensification factor
    momentum_term: bool = True     # include the rho*A*u^2 momentum term in the bend force
    dx: float = 0.5                # length of one hydraulic reach in the pressure solver (m)
    dt: float = 5e-4               # structural time step (s)
    T_end: float = 4.0             # length of the simulated record (s)
    n_modes: int = 0               # number of modes kept in modal solve; 0 = keep all

    # --- derived section properties (computed from the inputs above, never typed in) ---
    @property
    def As(self):                  # steel cross-section area (m^2)
        return np.pi / 4 * (self.Do**2 - self.Di**2)
    @property
    def Af(self):                  # fluid (bore) cross-section area (m^2)
        return np.pi / 4 * self.Di**2
    @property
    def I(self):                   # second moment of area of the pipe wall (m^4)
        return np.pi / 64 * (self.Do**4 - self.Di**4)
    @property
    def t_wall(self):              # wall thickness (m)
        return 0.5 * (self.Do - self.Di)
    @property
    def m_lin(self):               # mass per unit length = steel + the LNG inside (kg/m)
        return self.rho_s * self.As + self.rho_f * self.Af
    @property
    def dTc(self):                 # temperature drop during cooldown (K)
        return self.T0 - self.T1


def korteweg(p):
    """Wave speed corrected for pipe-wall elasticity (thin-wall form).
    A flexible pipe 'gives' under pressure, which slows the wave a little."""
    K = p.rho_f * p.c**2                                   # fluid bulk modulus from c = sqrt(K/rho)
    return p.c / np.sqrt(1 + K * p.Di / (p.E * p.t_wall))  # Korteweg formula


def bend_flex_factor(p):
    """Returns (k, i): Karman flexibility factor and stress intensification factor of a bend.
    Thin-walled bends ovalise, so they are MORE flexible and have HIGHER local stress than a straight beam."""
    rm = 0.5 * (p.Do - p.t_wall)            # mean radius of the pipe wall
    h = p.t_wall * p.R / rm**2              # pipe-bend characteristic h = tR/r^2
    return max(1.0, 1.65 / h), 0.9 / h ** (2 / 3)   # (flexibility factor k, SIF i), B31.3 App. D form


# ----------------------------------------------------------------------------
# 1. Geometry of the loop (polyline of straight beam elements; arcs by chords)
# ----------------------------------------------------------------------------
def build_loop(p):
    """Build the node coordinates, elements and bend data of a Pi-shaped expansion loop.
    Flow travels in +x, goes up, across, down, and continues in +x (four 90-degree bends)."""
    # path described as a list of segments: ("L", length) = straight, ("B", +1/-1) = left/right 90 deg bend
    segs = [("L", p.La), ("B", +1), ("L", p.H), ("B", -1), ("L", p.W),
            ("B", -1), ("L", p.H), ("B", +1), ("L", p.Lb)]
    pts, s, elems, bends = [(0.0, 0.0)], [0.0], [], []   # node (x,y) list; s = distance along pipe; elements; bend info
    heading = 0.0                                        # current travel direction angle (0 = +x)
    for kind, val in segs:                               # walk along the pipe segment by segment
        x, y = pts[-1]                                   # current end of the pipe
        s0 = s[-1]                                       # path length so far
        if kind == "L":                                  # ---- straight section ----
            n = max(1, int(np.ceil(val / p.h_max)))      # how many elements to use
            for k in range(1, n + 1):
                d = val * k / n                          # distance from start of this straight
                pts.append((x + d * np.cos(heading), y + d * np.sin(heading)))  # new node position
                s.append(s0 + d)                         # path length of the new node
                elems.append((len(pts) - 2, len(pts) - 1, False))   # element joining last two nodes (not a bend)
        else:                                            # ---- 90 degree bend ----
            sg, n, h0 = val, p.n_arc, heading            # sg = +1 left turn / -1 right turn; h0 = heading before bend
            cx, cy = x - sg * p.R * np.sin(h0), y + sg * p.R * np.cos(h0)   # centre of curvature of the arc
            start = len(pts) - 1                         # index of the node where the bend starts
            for k in range(1, n + 1):
                a = sg * k * (np.pi / 2) / n             # angle turned so far
                h = h0 + a                               # heading at this point of the arc
                pts.append((cx + sg * p.R * np.sin(h), cy - sg * p.R * np.cos(h)))  # point on the arc
                s.append(s0 + p.R * abs(a))              # path length = arc length so far
                elems.append((len(pts) - 2, len(pts) - 1, True))   # chord element, flagged as a bend
            heading = h0 + sg * np.pi / 2                # new travel direction after the bend
            n1 = np.array([np.cos(h0), np.sin(h0)])      # unit vector of the flow entering the bend
            n2 = np.array([np.cos(heading), np.sin(heading)])  # unit vector of the flow leaving the bend
            end = len(pts) - 1                           # index of the node where the bend ends
            bends.append(dict(start=start, end=end, mid=start + n // 2, n1=n1, n2=n2,
                              s_start=s[start], s_end=s[end],        # hydraulic stations of the two arc ends
                              Ps=np.array(pts[start]), Pe=np.array(pts[end]),  # coordinates of the arc ends
                              T=np.array(pts[start]) + p.R * n1))    # tangent intersection (distance = R for 90 deg)
    n_nodes = len(pts)
    # anchors: first and last node fully fixed (3 dofs each: ux, uy, rotation)
    fixed = [0, 1, 2, 3 * (n_nodes - 1), 3 * (n_nodes - 1) + 1, 3 * (n_nodes - 1) + 2]
    sa = np.array(s)
    for (sg_pos, dof) in p.guides:                       # optional lateral guides: fix one dof at the nearest node
        fixed.append(3 * int(np.argmin(abs(sa - sg_pos))) + dof)
    return dict(pts=np.array(pts), s=sa, elems=elems, bends=bends, fixed=sorted(set(fixed)),
                L_hyd=sa[-1])                            # L_hyd = total pipe length reservoir -> valve


def build_cantilever(p, L, n=40):
    """A straight clamped-free beam of length L, used only to verify the beam elements."""
    pts = np.column_stack([np.linspace(0, L, n + 1), np.zeros(n + 1)])   # nodes along x
    return dict(pts=pts, s=pts[:, 0], elems=[(i, i + 1, False) for i in range(n)], bends=[],
                fixed=[0, 1, 2], L_hyd=L)                # clamp = fix all 3 dofs of node 0


# ----------------------------------------------------------------------------
# 2. Frame FE (2D Euler-Bernoulli, consistent mass, fluid mass lumped into m_lin)
# ----------------------------------------------------------------------------
def _kloc(EA, EI, L):
    """Local 6x6 stiffness matrix of a 2D beam element (axial + bending).
    dof order per node: [axial u, transverse v, rotation theta]."""
    a, b, c_, d = EA / L, 12 * EI / L**3, 6 * EI / L**2, 4 * EI / L   # standard beam coefficients
    return np.array([[a, 0, 0, -a, 0, 0], [0, b, c_, 0, -b, c_], [0, c_, d, 0, -c_, d / 2],
                     [-a, 0, 0, a, 0, 0], [0, -b, -c_, 0, b, -c_], [0, c_, d / 2, 0, -c_, d]])


def _mloc(m, L):
    """Local 6x6 CONSISTENT mass matrix of a beam element with mass per length m."""
    return m * L / 420 * np.array([[140, 0, 0, 70, 0, 0], [0, 156, 22 * L, 0, 54, -13 * L],
                                   [0, 22 * L, 4 * L**2, 0, 13 * L, -3 * L**2],
                                   [70, 0, 0, 140, 0, 0], [0, 54, 13 * L, 0, 156, -22 * L],
                                   [0, -13 * L, -3 * L**2, 0, -22 * L, 4 * L**2]])


class System:
    """Assembles global stiffness K, mass M and thermal load vector for a frame geometry."""
    def __init__(self, p, geo, thermal=True):
        self.p, self.geo = p, geo
        pts = geo["pts"]
        n = len(pts)                                     # number of nodes
        nd = 3 * n                                       # total degrees of freedom
        K, M, fth = np.zeros((nd, nd)), np.zeros((nd, nd)), np.zeros(nd)   # empty global matrices
        kflex, _ = bend_flex_factor(p)                   # flexibility factor (used only if flex_bends)
        self.el = []                                     # per-element data kept for stress recovery
        eps0 = -p.alpha * p.dTc if thermal else 0.0      # free thermal strain (negative = contraction)
        for (i, j, bend) in geo["elems"]:                # loop over every element
            d = pts[j] - pts[i]                          # vector from node i to node j
            L = np.hypot(*d)                             # element length
            c, s_ = d / L                                # cos and sin of the element angle
            T = np.zeros((6, 6))                         # rotation matrix: global -> local axes
            T[:2, :2] = T[3:5, 3:5] = [[c, s_], [-s_, c]]
            T[2, 2] = T[5, 5] = 1                        # rotations are unaffected by the axis change
            EI = p.E * p.I / (kflex if (bend and p.flex_bends) else 1.0)   # softer bends if requested
            kl, ml = _kloc(p.E * p.As, EI, L), _mloc(p.m_lin, L)           # local stiffness and mass
            fl = p.E * p.As * eps0 * np.array([-1, 0, 0, 1, 0, 0])         # local thermal load (initial strain)
            dofs = np.r_[3 * i:3 * i + 3, 3 * j:3 * j + 3]                 # global dof numbers of this element
            K[np.ix_(dofs, dofs)] += T.T @ kl @ T        # add element stiffness (rotated to global) into K
            M[np.ix_(dofs, dofs)] += T.T @ ml @ T        # add element mass into M
            fth[dofs] += T.T @ fl                        # add element thermal load into the global vector
            self.el.append(dict(dofs=dofs, T=T, kl=kl, fl=fl, bend=bend))
        self.K, self.M, self.fth, self.nd = K, M, fth, nd
        self.fix = np.array(geo["fixed"])                # constrained dofs (anchors, guides)
        self.free = np.setdiff1d(np.arange(nd), self.fix)  # unconstrained dofs
        self.Kff = K[np.ix_(self.free, self.free)]       # stiffness of the free dofs only
        self.Mff = M[np.ix_(self.free, self.free)]       # mass of the free dofs only

    def modes(self, n=None):
        """Natural frequencies and mode shapes: solve K phi = w^2 M phi.
        Shapes are mass-normalised (phi^T M phi = 1), which simplifies the modal equations."""
        w2, Phi = eigh(self.Kff, self.Mff)
        w = np.sqrt(np.maximum(w2, 0))                   # angular frequencies (rad/s); max() guards rounding
        if n:
            w, Phi = w[:n], Phi[:, :n]
        return w, Phi

    def rayleigh(self, w, zeta, i=0, j=5):
        """Rayleigh damping C = alpha*M + beta*K, tuned to give damping ratio zeta at modes i and j."""
        wa, wb = w[i], w[min(j, len(w) - 1)]
        return 2 * zeta * wa * wb / (wa + wb), 2 * zeta / (wa + wb)   # alpha, beta


# ----------------------------------------------------------------------------
# 3. Hydraulics: method of characteristics (replace with your Stage-1 solver output)
# ----------------------------------------------------------------------------
class Hyd:
    """Holds pressure P and velocity U on a (time, station) grid and interpolates in both."""
    def __init__(self, t, s, P, U):
        self.t, self.s, self.P, self.U = t, s, P, U

    def at(self, s_query, tq, what="P"):
        """Time history of pressure (or velocity) at distance s_query from the reservoir,
        evaluated at times tq. Linear interpolation between stations, then in time."""
        A = self.P if what == "P" else self.U
        x = np.clip(s_query / (self.s[1] - self.s[0]), 0, len(self.s) - 1 - 1e-9)  # fractional station index
        i0 = int(x)                                       # station just before
        w = x - i0                                        # interpolation weight
        series = (1 - w) * A[:, i0] + w * A[:, i0 + 1]    # blend neighbouring stations
        return np.interp(tq, self.t, series)              # resample onto the requested times


def load_stage1(path, L_total=None, time_col=0):
    """Load your group's Stage-1 pressure histories into a Hyd object.
    Accepted formats
      * .npz with arrays t (nt,), s (ns,), P (nt, ns) [absolute Pa]  (optional U (nt, ns))
      * .csv  first row = header: 'time', then the STATION POSITION in metres for every pressure column
              e.g.  time, 0, 10, 25, 50, 100      (pressures in Pa, absolute)
              If your file is in MPa or gauge, convert BEFORE calling (see the checks below).
    Stations need not be evenly spaced; they are re-gridded onto an even grid by linear interpolation.
    NOTE: the more stations you export the better - interpolating between a few stations smears the
    wave front. Ideally export the full grid your solver already computes."""
    if path.endswith(".npz"):
        d = np.load(path)
        t, s, P = d["t"], d["s"], d["P"]
        U = d["U"] if "U" in d else np.full_like(P, 0.0)
    else:
        raw = np.genfromtxt(path, delimiter=",", names=None, skip_header=1)
        hdr = open(path).readline().strip().split(",")
        t = raw[:, time_col]
        cols = [i for i in range(raw.shape[1]) if i != time_col]
        s = np.array([float(hdr[i]) for i in cols])
        P = raw[:, cols]
        U = np.zeros_like(P)                    # velocity not exported -> momentum term will be ~0 (it is tiny)
    order = np.argsort(s)
    s, P, U = s[order], P[:, order], U[:, order]
    L = L_total or s[-1]
    sg = np.linspace(0, L, max(len(s), int(L / 0.5) + 1))     # even grid the Hyd class expects
    Pg = np.array([np.interp(sg, s, P[k]) for k in range(len(t))])
    Ug = np.array([np.interp(sg, s, U[k]) for k in range(len(t))])
    # sanity checks you should look at
    if P.max() < 1e5:
        print("WARNING: max pressure < 0.1 MPa - is your file in MPa? multiply by 1e6.")
    print(f"Loaded Stage 1: {len(t)} time steps, dt={t[1]-t[0]:.2e} s, {len(s)} stations over {s[-1]:.1f} m, "
          f"P range {P.min()/1e6:.3f} to {P.max()/1e6:.3f} MPa")
    return Hyd(t, sg, Pg, Ug)


class HydAtStations(Hyd):
    """Hydraulic results at a fixed list of stations (e.g. the 8 bend arc-end stations), as exported by
    the Stage-1 solver. at() only answers for those stations - it does NOT interpolate along the pipe."""
    def __init__(self, t, s, P, U):
        super().__init__(t, s, P, U)
    def at(self, s_query, tq, what="P"):
        A = self.P if what == "P" else self.U
        k = int(np.argmin(abs(self.s - s_query)))
        if abs(self.s[k] - s_query) > 1e-6:
            raise ValueError(f"station {s_query:.3f} m was not exported (have {self.s})")
        return np.interp(tq, self.t, A[:, k])


def hyd_from_stage1_solver(stage1_repo, geo, speed_mode="bulk", duration=4.0, reaches=None, station_tol=1e-6):
    """Run the GROUP'S Stage-1 solver (MinhAnh120/MECH0074, branch stage1-pressure-waves) on the actual loop
    route and return pressures/velocities at every bend's inlet and outlet station.
      stage1_repo : path of the cloned Stage-1 repository
      geo         : output of build_loop(p)  (gives total hydraulic length and bend arc-end stations)
      speed_mode  : 'bulk' (1100 m/s) or 'elastic' (997.9 m/s) - must match the group decision
    Uses the Stage-1 nominal Darcy friction, fixed reservoir, linear flow closure. Records the phase screen."""
    import sys, math
    from dataclasses import replace as _rep
    sys.path.insert(0, stage1_repo)
    from complete_stage1 import shared_case
    from fluid_assumptions import friction_estimate
    from stage1 import simulate
    base = shared_case(speed_mode, 0)
    f = friction_estimate(base.density, base.velocity, base.diameter)["darcy_factor"]
    case = _rep(base, friction=f, length=geo["L_hyd"])
    stations = sorted({round(float(x), 9) for b in geo["bends"] for x in (b["s_start"], b["s_end"])})
    # Stage-1 requires stations strictly inside [0, L]; choose N so every station is (nearly) a node
    if reaches is None:
        reaches = int(round(geo["L_hyd"] / 0.3125))
    r = simulate(case, reaches, duration, stations=stations)
    A = case.area
    print(f"Stage-1 solver: N={reaches}, dt={r['dt']:.3e} s, dx={r['dx']:.4f} m, interpolated stations={r['interpolated_stations']}, "
          f"min p_abs={r['minimum_absolute_pa']/1e6:.3f} MPa, phase status={r['phase_status']}")
    return HydAtStations(r["time"], np.array(stations), r["pressure"], r["flow"] / A), r


def moc(p, L, p_floor=None):
    """1D water-hammer solver. Reservoir at x=0, valve at x=L, valve closes linearly over p.t_close.
    Method of characteristics: information travels along lines dx/dt = +/-c, so with dt = dx/c each new
    node value follows exactly from its two neighbours at the previous step."""
    N = max(2, int(round(L / p.dx)))                     # number of reaches
    dx = L / N                                           # actual reach length
    dt = dx / p.c                                        # hydraulic time step (Courant number = 1)
    nt = int(p.T_end / dt) + 1                           # number of time steps
    rho, c, D = p.rho_f, p.c, p.Di
    s = np.linspace(0, L, N + 1)                         # station positions
    F = lambda u: p.f_darcy * rho / (2 * D) * u * np.abs(u)   # Darcy friction term (zero if f_darcy=0)
    P = np.empty((nt, N + 1))                            # pressure[time, station]
    U = np.empty((nt, N + 1))                            # velocity[time, station]
    U[0] = p.u0                                          # initial condition: steady flow...
    P[0] = p.p_res - F(p.u0) * s                         # ...with pressure = reservoir pressure (minus friction drop)
    for n in range(nt - 1):                              # march forward in time
        pa, ua = P[n, :-2], U[n, :-2]                    # values at node i-1 (feeds the C+ characteristic)
        pb, ub = P[n, 2:], U[n, 2:]                      # values at node i+1 (feeds the C- characteristic)
        Fa, Fb = F(ua), F(ub)
        Un, Pn = np.empty(N + 1), np.empty(N + 1)        # new time level
        # interior nodes: intersect the C+ and C- characteristic equations
        Un[1:-1] = (pa - pb + rho * c * (ua + ub) - c * dt * (Fa + Fb)) / (2 * rho * c)
        Pn[1:-1] = 0.5 * (pa + pb + rho * c * (ua - ub) - c * dt * (Fa - Fb))
        Pn[0] = p.p_res                                  # reservoir boundary: pressure held constant
        Un[0] = U[n, 1] + (Pn[0] - P[n, 1] - c * dt * F(U[n, 1])) / (rho * c)   # velocity from the C- line
        t1 = (n + 1) * dt                                # time of the new level
        Un[-1] = p.u0 * max(0.0, 1 - t1 / p.t_close)     # valve boundary: velocity ramps linearly to zero
        Pn[-1] = P[n, -2] - rho * c * (Un[-1] - U[n, -2]) - c * dt * F(U[n, -2])   # pressure from the C+ line
        if p_floor is not None:                          # crude vapour clip (NOT a proper cavitation model)
            Pn = np.maximum(Pn, p_floor)
        P[n + 1], U[n + 1] = Pn, Un                      # store the new level
    return Hyd(np.arange(nt) * dt, s, P, U)


# ----------------------------------------------------------------------------
# 4. Bend forces
# ----------------------------------------------------------------------------
def bend_force_vector(dp_in, dp_out, u_in, u_out, n1, n2, A, rho, momentum=True):
    """Force of the fluid ON the pipe at a bend (from a momentum balance on the fluid):
         F = [(dp1)A + rho A u1^2] n1  -  [(dp2)A + rho A u2^2] n2
    dp = gauge pressure (p - p_ext); n1/n2 = unit vectors of flow in/out.
    Returned as two separate forces so each can be applied on its own leg line."""
    m1 = rho * A * u_in**2 if momentum else 0.0          # momentum-flux term at the inlet (small here)
    m2 = rho * A * u_out**2 if momentum else 0.0         # ... and at the outlet
    return np.outer(n1, dp_in * A + m1), -np.outer(n2, dp_out * A + m2)   # shape (2, time) each


def check_uniform_bend(p, bend, dp=1.0e6):
    """Sanity check with UNIFORM pressure: the resultant force must be sqrt(2)*dp*A, point along the
    bisector (outwards), and have zero moment about the tangent intersection."""
    A = p.Af
    F1, F2 = dp * A * bend["n1"], -dp * A * bend["n2"]   # inlet and outlet forces
    Fr = F1 + F2                                         # resultant force
    bis = (bend["n1"] - bend["n2"]) / np.linalg.norm(bend["n1"] - bend["n2"])   # bisector direction
    cross = lambda r, f: r[0] * f[1] - r[1] * f[0]       # 2D moment r x f
    Mt = cross(bend["Ps"] - bend["T"], F1) + cross(bend["Pe"] - bend["T"], F2)  # moment about tangent intersection
    return dict(F_resultant=np.linalg.norm(Fr), F_expected=np.sqrt(2) * dp * A,
                angle_to_bisector_deg=np.degrees(np.arccos(np.clip(Fr @ bis / np.linalg.norm(Fr), -1, 1))),
                moment_about_tangent_intersection=Mt)


def assemble_forces(p, sysm, hyd, t):
    """Build the force history on every structural dof, F[dof, time].
    Each bend's inlet and outlet forces use the pressure at THEIR OWN hydraulic stations, so wave travel
    time between bends (and across each arc) is retained automatically."""
    A, rho = p.Af, p.rho_f
    F = np.zeros((sysm.nd, len(t)))
    for b in sysm.geo["bends"]:
        pa, pe = hyd.at(b["s_start"], t) - p.p_ext, hyd.at(b["s_end"], t) - p.p_ext   # gauge pressures
        ua, ue = hyd.at(b["s_start"], t, "U"), hyd.at(b["s_end"], t, "U")             # velocities
        Fi, Fo = bend_force_vector(pa, pe, ua, ue, b["n1"], b["n2"], A, rho, p.momentum_term)
        F[3 * b["start"]:3 * b["start"] + 2] += Fi       # apply inlet force (x,y) at the arc-start node
        F[3 * b["end"]:3 * b["end"] + 2] += Fo           # apply outlet force (x,y) at the arc-end node
    return F


# ----------------------------------------------------------------------------
# 5. Dynamic solvers
# ----------------------------------------------------------------------------
def newmark_direct(M, C, K, F, dt, beta=0.25, gamma=0.5):
    """Newmark time integration of M a + C v + K u = F(t) for the full model.
    beta=1/4, gamma=1/2 is the 'average acceleration' scheme: unconditionally stable, no numerical damping."""
    n, nt = F.shape
    u, v, a = np.zeros((n, nt)), np.zeros((n, nt)), np.zeros((n, nt))   # displacement, velocity, acceleration
    a0, a1, a2 = 1 / (beta * dt**2), gamma / (beta * dt), 1 / (beta * dt)   # Newmark constants
    a3, a4, a5 = 1 / (2 * beta) - 1, gamma / beta - 1, dt / 2 * (gamma / beta - 2)
    lu = lu_factor(K + a1 * C + a0 * M)                  # factorise the effective stiffness ONCE (it never changes)
    a[:, 0] = np.linalg.solve(M, F[:, 0] - C @ v[:, 0] - K @ u[:, 0])   # initial acceleration from equilibrium
    for k in range(nt - 1):
        rhs = (F[:, k + 1] + M @ (a0 * u[:, k] + a2 * v[:, k] + a3 * a[:, k])
               + C @ (a1 * u[:, k] + a4 * v[:, k] + a5 * a[:, k]))      # effective load for the next step
        u[:, k + 1] = lu_solve(lu, rhs)                  # new displacement
        a[:, k + 1] = a0 * (u[:, k + 1] - u[:, k]) - a2 * v[:, k] - a3 * a[:, k]   # new acceleration
        v[:, k + 1] = v[:, k] + dt * ((1 - gamma) * a[:, k] + gamma * a[:, k + 1]) # new velocity
    return u, v, a


def newmark_modal(w, zeta_n, f_modal, dt, beta=0.25, gamma=0.5):
    """Same integration, but on the decoupled MODAL equations (one scalar equation per mode):
         q'' + 2 zeta_n w_n q' + w_n^2 q = phi_n^T F(t)       (mass-normalised modes)
    Much cheaper, and it shows which modes are driven. f_modal[mode, time] = Phi^T F."""
    nm, nt = f_modal.shape
    q, qd, qdd = (np.zeros((nm, nt)) for _ in range(3))  # modal coordinate, velocity, acceleration
    a0, a1, a2 = 1 / (beta * dt**2), gamma / (beta * dt), 1 / (beta * dt)
    a3, a4, a5 = 1 / (2 * beta) - 1, gamma / beta - 1, dt / 2 * (gamma / beta - 2)
    cn = 2 * zeta_n * w                                  # modal damping coefficients
    keff = w**2 + a1 * cn + a0                           # effective stiffness: a scalar per mode
    qdd[:, 0] = f_modal[:, 0]                            # starts at rest, so q'' = modal force
    for k in range(nt - 1):
        rhs = (f_modal[:, k + 1] + (a0 * q[:, k] + a2 * qd[:, k] + a3 * qdd[:, k])
               + cn * (a1 * q[:, k] + a4 * qd[:, k] + a5 * qdd[:, k]))
        q[:, k + 1] = rhs / keff                         # divide, not solve: modes are decoupled
        qdd[:, k + 1] = a0 * (q[:, k + 1] - q[:, k]) - a2 * qd[:, k] - a3 * qdd[:, k]
        qd[:, k + 1] = qd[:, k] + dt * ((1 - gamma) * qdd[:, k] + gamma * qdd[:, k + 1])
    return q, qd, qdd


# ----------------------------------------------------------------------------
# 6. Post-processing: stress and reactions
# ----------------------------------------------------------------------------
def member_stress(p, sysm, u_tot, p_gauge_peak):
    """Nominal stress envelope over all elements and all times.
    Conservative longitudinal stress = |axial| + |bending|, combined with hoop stress by von Mises."""
    _, sif = bend_flex_factor(p)
    ro = 0.5 * p.Do                                      # outer radius (where bending stress is largest)
    sig_h = p_gauge_peak * p.Di / (2 * p.t_wall)         # hoop stress from the peak pressure: p D / 2t
    peak, where = 0.0, -1
    sig_ax_max = sig_b_max = 0.0
    for e, el in enumerate(sysm.el):
        ul = el["T"] @ u_tot[el["dofs"], :]              # element displacements in local axes
        f = el["kl"] @ ul - el["fl"][:, None]            # internal forces (minus the free thermal strain part)
        N = f[3]                                         # axial force
        Mb = np.maximum(abs(f[2]), abs(f[5]))            # larger of the two end bending moments
        i_f = sif if (el["bend"] and p.sif) else 1.0     # apply SIF only on bends and only if switched on
        sl = abs(N) / p.As + i_f * Mb * ro / p.I         # longitudinal stress: N/A + M r / I
        vm = np.sqrt(sl**2 - sl * sig_h + sig_h**2)      # von Mises for a biaxial (longitudinal + hoop) state
        if vm.max() > peak:                              # remember the worst element
            peak, where = vm.max(), e
            sig_ax_max, sig_b_max = (abs(N) / p.As).max(), (i_f * Mb * ro / p.I).max()
    return dict(vm_peak=peak, element=where, axial=sig_ax_max, bending=sig_b_max, hoop=sig_h)


def run_case(p, label="case", method="modal", p_floor=None, hyd=None, verbose=True):
    """Run one complete design: geometry -> FE -> pressure -> forces -> response -> metrics."""
    geo = build_loop(p)                                  # nodes, elements, bends
    sysm = System(p, geo)                                # K, M, thermal load
    fr, fx = sysm.free, sysm.fix                         # free / fixed dof lists
    if hyd is None:
        hyd = moc(p, geo["L_hyd"], p_floor)              # pressure solution (or use the one you pass in)
    t = np.arange(0, p.T_end, p.dt)                      # structural time grid
    F = assemble_forces(p, sysm, hyd, t)                 # pressure -> bend forces on every dof
    F0 = F[:, 0].copy()                                  # steady thrust before the valve moves
    dF = F - F0[:, None]                                 # surge INCREMENT (the steady part is static, not vibration)
    # --- static state: thermal contraction + steady pressure thrust ---
    u_th = np.zeros(sysm.nd)
    u_th[fr] = np.linalg.solve(sysm.Kff, sysm.fth[fr])             # thermal only (gives 'thermal movement')
    u_st = np.zeros(sysm.nd)
    u_st[fr] = np.linalg.solve(sysm.Kff, sysm.fth[fr] + F0[fr])    # thermal + steady pressure
    # --- dynamics: modes, damping, then time integration ---
    w, Phi = sysm.modes()
    al, be = sysm.rayleigh(w, p.zeta)                    # Rayleigh coefficients
    zn = 0.5 * (al / w + be * w)                         # damping ratio of every mode implied by them
    if method == "modal":
        nm = p.n_modes or len(w)                         # how many modes to keep
        fm = Phi[:, :nm].T @ dF[fr]                      # MODAL FORCES: project force onto each mode shape
        q, qd, qdd = newmark_modal(w[:nm], zn[:nm], fm, p.dt)
        ud, vd, ad = Phi[:, :nm] @ q, Phi[:, :nm] @ qd, Phi[:, :nm] @ qdd   # back to physical dofs
    else:
        C = al * sysm.Mff + be * sysm.Kff                # full damping matrix
        ud, vd, ad = newmark_direct(sysm.Mff, C, sysm.Kff, dF[fr], p.dt)
    u_dyn = np.zeros((sysm.nd, len(t)))                  # expand free-dof results to all dofs
    u_dyn[fr] = ud
    v_full, a_full = np.zeros_like(u_dyn), np.zeros_like(u_dyn)
    v_full[fr], a_full[fr] = vd, ad
    u_tot = u_st[:, None] + u_dyn                        # total = static + surge vibration
    # --- peak surge displacement: biggest in-plane displacement magnitude over all nodes and times ---
    mag = np.hypot(u_dyn[0::3], u_dyn[1::3])             # |(ux,uy)| at every node and time
    k_node, k_t = np.unravel_index(np.argmax(mag), mag.shape)
    # --- support reactions at anchors: R = M a + C v + K u - F - f_thermal, on the fixed rows ---
    Cfull = al * sysm.M + be * sysm.K
    R = (sysm.M[fx] @ a_full + Cfull[fx] @ v_full + sysm.K[fx] @ u_tot
         - (F[fx] + sysm.fth[fx][:, None]))
    Rmag = np.hypot(R[0], R[1])                          # force magnitude at anchor 1 (node 0)
    Rmag2 = np.hypot(R[3], R[4])                         # force magnitude at anchor 2 (last node)
    Rpk = max(Rmag.max(), Rmag2.max())                   # peak anchor force over time
    Rst = max(np.hypot(*(R[0:2, 0])), np.hypot(*(R[3:5, 0])))   # static anchor force (t=0)
    pg_peak = max(0.0, hyd.P.max() - p.p_ext)            # peak gauge pressure, for hoop stress
    st = member_stress(p, sysm, u_tot, pg_peak)          # nominal stress envelope
    res = dict(label=label, p=p, geo=geo, sys=sysm, hyd=hyd, t=t, w=w, Phi=Phi,
               F=F, dF=dF, u_dyn=u_dyn, u_th=u_th, k_node=k_node, k_t=k_t,
               therm_move_mm=1e3 * np.abs(np.hypot(u_th[0::3], u_th[1::3])).max(),   # largest thermal movement
               contraction_mm=1e3 * p.alpha * p.dTc * geo["pts"][-1, 0],             # free contraction between anchors
               surge_mm=1e3 * mag.max(), react_peak_kN=Rpk / 1e3, react_static_kN=Rst / 1e3,
               stress=st, dp_joukowsky=p.rho_f * p.c * p.u0,
               dp_peak_loop=max(hyd.at(b["s_end"], t).max() - p.p_res for b in geo["bends"]))
    if verbose:
        f = w[:6] / 2 / np.pi                            # first six natural frequencies in Hz
        print(f"[{label}] f1-6 (Hz): " + " ".join(f"{x:6.2f}" for x in f))
    return res


# ----------------------------------------------------------------------------
# 7. Verification helpers
# ----------------------------------------------------------------------------
def cantilever_check(p, L, n_el=40, n_modes=6):
    """Compare FE (fluid-filled) cantilever bending frequencies with Euler-Bernoulli theory:
         f_n = (beta_n L)^2 / (2 pi L^2) * sqrt(EI / m)      (m includes the fluid mass)"""
    sysm = System(p, build_cantilever(p, L, n_el), thermal=False)
    w, Phi = sysm.modes()
    ux = Phi[0::3]                                       # axial component of every mode shape
    axial = (ux**2).sum(0) / (Phi**2).sum(0) > 0.5       # a mode is 'axial' if most of its motion is axial
    fe = (w[~axial] / 2 / np.pi)[:n_modes]               # keep bending modes only, in Hz
    bl = np.array([1.87510407, 4.69409113, 7.85475744, 10.99554073, 14.13716839, 17.27875953])  # roots of cantilever equation
    th = bl[:len(fe)] ** 2 / (2 * np.pi * L**2) * np.sqrt(p.E * p.I / p.m_lin)
    return fe, th


def joukowsky_numbers(p, L=100.0):
    """Hand-calculation numbers for Stage 1: surge pressure, periods and the 'rapid closure' limit."""
    return dict(dp_MPa=p.rho_f * p.c * p.u0 / 1e6, dp_korteweg_MPa=p.rho_f * korteweg(p) * p.u0 / 1e6,
                period_4L_over_c_s=4 * L / p.c, travel_time_L_over_c_s=L / p.c,
                rapid_limit_2L_over_c_s=2 * L / p.c, c_korteweg=korteweg(p))


# ----------------------------------------------------------------------------
# 8. Plots
# ----------------------------------------------------------------------------
def plot_pressure(res, fn):
    """Pressure histories at each bend inlet and at the valve, with the Joukowsky level marked."""
    p, hyd, geo, t = res["p"], res["hyd"], res["geo"], res["t"]
    fig, ax = plt.subplots(figsize=(8, 4))
    for i, b in enumerate(geo["bends"]):
        ax.plot(t, hyd.at(b["s_start"], t) / 1e6, label=f"bend {i+1} inlet (s={b['s_start']:.1f} m)")
    try:                                       # valve trace only if the hydraulics include that station
        ax.plot(t, hyd.at(geo["L_hyd"], t) / 1e6, "k--", lw=.8, label="valve")
    except ValueError:
        pass
    ax.axhline((p.p_res + res["dp_joukowsky"]) / 1e6, color="r", ls=":", label=r"p0+$\rho c \Delta u$")
    ax.set(xlim=(0, 1.0), xlabel="t (s)", ylabel="absolute pressure (MPa)")
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def plot_spectra(res, fn, n_show=6):
    """Frequency content of the modal forces vs the loop natural frequencies and the pressure-wave
    harmonics (2k-1)c/4L. Peaks sitting on the dashed lines = resonance risk."""
    p, t, w, Phi, sysm = res["p"], res["t"], res["w"], res["Phi"], res["sys"]
    fr = sysm.free
    fm = Phi[:, :n_show].T @ res["dF"][fr]               # modal force histories of the first modes
    f = np.fft.rfftfreq(len(t), p.dt)                    # frequency axis (Hz)
    fig, ax = plt.subplots(figsize=(8, 4))
    for k in range(n_show):
        ax.semilogy(f, np.abs(np.fft.rfft(fm[k])) * p.dt, lw=.8, label=f"modal force {k+1}")
    for k in range(n_show):
        ax.axvline(w[k] / 2 / np.pi, color="k", ls="--", lw=.6)       # natural frequencies
    L = res["geo"]["L_hyd"]
    for k in range(1, 8):
        ax.axvline((2 * k - 1) * p.c / (4 * L), color="r", ls=":", lw=.8)   # wave harmonics
    ax.set(xlim=(0, 40), xlabel="frequency (Hz)", ylabel=r"modal force amplitude spectrum $|\Phi_n^T \Delta F|$ (arb.)",
           title="dashed black: loop modes 1-6   dotted red: (2k-1)c/4L")
    ax.legend(fontsize=6, ncol=3)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def plot_response(results, fn):
    """Surge displacement time history of the most-excited node, for each design on one plot."""
    fig, ax = plt.subplots(figsize=(8, 4))
    for r in results:
        k = r["k_node"]
        u = np.hypot(r["u_dyn"][3 * k], r["u_dyn"][3 * k + 1]) * 1e3   # displacement magnitude in mm
        ax.plot(r["t"], u, lw=.8, label=f"{r['label']} (peak {r['surge_mm']:.1f} mm)")
    ax.set(xlabel="t (s)", ylabel="surge displacement at most-excited node (mm)", xlim=(0, 4))
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def plot_modes(res, fn, n=6):
    """Draw the first n mode shapes (grey = undeformed loop, blue = exaggerated mode shape)."""
    sysm, w, Phi, pts = res["sys"], res["w"], res["Phi"], res["geo"]["pts"]
    fig, axs = plt.subplots(2, 3, figsize=(11, 5))
    for k, ax in enumerate(axs.ravel()[:n]):
        ph = np.zeros(sysm.nd)
        ph[sysm.free] = Phi[:, k]                        # expand mode to all dofs
        d = np.column_stack([ph[0::3], ph[1::3]])        # (ux, uy) at each node
        d = d / np.abs(d).max() * 2.5                    # scale so the largest motion is 2.5 m on the plot
        ax.plot(*pts.T, "0.7", lw=1)
        ax.plot(*(pts + d).T, "C0", lw=1.5)
        ax.set_aspect("equal")
        ax.set_title(f"mode {k+1}: {w[k]/2/np.pi:.2f} Hz", fontsize=9)
    fig.tight_layout()
    fig.savefig(fn, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
# 9. Studies
# ----------------------------------------------------------------------------
def sensitivity(base, which):
    """Repeat the baseline while varying ONE numerical setting; if the results barely move the
    setting is 'converged' (this is your evidence for mesh / resolution / damping checks)."""
    out = []
    if which == "damping":                               # how strongly does the assumed damping matter?
        for z in (0.005, 0.01, 0.02, 0.05):
            r = run_case(replace(base, zeta=z), f"zeta={z}", verbose=False)
            out.append((f"zeta={z}", r["surge_mm"], r["stress"]["vm_peak"] / 1e6, r["react_peak_kN"]))
    if which == "hydraulic_dx":                          # pressure-solver resolution
        for dx in (2.0, 1.0, 0.5, 0.25):
            r = run_case(replace(base, dx=dx), f"dx={dx}", verbose=False)
            out.append((f"dx={dx} m", r["surge_mm"], r["stress"]["vm_peak"] / 1e6, r["react_peak_kN"]))
    if which == "structural_mesh":                       # structural mesh refinement
        for h, n in ((2.0, 4), (1.0, 8), (0.5, 16)):
            r = run_case(replace(base, h_max=h, n_arc=n), f"h={h}", verbose=False)
            out.append((f"h_max={h} m, {n} chords/bend", r["surge_mm"], r["stress"]["vm_peak"] / 1e6, r["react_peak_kN"]))
    if which == "time_step":                             # time-integration step
        for dt in (2e-3, 1e-3, 5e-4, 2.5e-4):
            r = run_case(replace(base, dt=dt), f"dt={dt}", verbose=False)
            out.append((f"dt={dt*1e3:.2f} ms", r["surge_mm"], r["stress"]["vm_peak"] / 1e6, r["react_peak_kN"]))
    print(f"\nSensitivity: {which}\n  {'setting':32s} surge(mm)  VM peak(MPa)  reaction(kN)")
    for row in out:
        print(f"  {row[0]:32s} {row[1]:9.2f} {row[2]:12.1f} {row[3]:13.0f}")
    return out


# ----------------------------------------------------------------------------
# 10. Demo  (runs only when you execute this file directly)
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    import csv
    base = Par()                                         # baseline design = the teaching case
    # hand-calculation numbers to quote on your slides
    print("Theory:", {k: round(float(v), 4) for k, v in joukowsky_numbers(base).items()})
    k, i = bend_flex_factor(base)
    print(f"Bend flexibility factor k={k:.1f}, SIF i={i:.2f} (set flex_bends/sif True to use)")
    print(f"Fully restrained thermal stress E*alpha*dT = {base.E*base.alpha*base.dTc/1e6:.0f} MPa")

    # verification 1: beam elements vs Euler-Bernoulli for 4, 10, 20 m fluid-filled cantilevers
    print("\nBeam-element verification, fluid-filled cantilevers (FE vs Euler-Bernoulli, Hz):")
    for L in (4.0, 10.0, 20.0):
        fe, th = cantilever_check(base, L)
        print(f"  L={L:4.0f} m  FE: " + " ".join(f"{x:7.3f}" for x in fe))
        print(f"           EB: " + " ".join(f"{x:7.3f}" for x in th))

    # verification 2: bend force direction and moment balance with uniform pressure
    print("\nUniform-pressure bend check:", check_uniform_bend(base, build_loop(base)["bends"][0]))

    # the designs to compare: baseline + three modifications (edit these!)
    cases = [("baseline", base),
             ("slow closure 0.3 s", replace(base, t_close=0.3)),    # modification 1: closure time
             ("larger radius 2.5 m", replace(base, R=2.5)),         # modification 2: bend curvature
             ("larger loop 14x14 m", replace(base, H=14.0, W=14.0))]  # modification 3: loop size
    results = [run_case(p, lab) for lab, p in cases]

    # comparison table: the four quantities the brief asks for
    print(f"\n{'case':22s}{'thermal mv(mm)':>15s}{'surge pk(mm)':>14s}{'anchor pk(kN)':>15s}{'static(kN)':>12s}{'VM pk(MPa)':>12s}")
    rows = []
    for r in results:
        print(f"{r['label']:22s}{r['therm_move_mm']:15.1f}{r['surge_mm']:14.2f}{r['react_peak_kN']:15.0f}"
              f"{r['react_static_kN']:12.0f}{r['stress']['vm_peak']/1e6:12.1f}")
        rows.append([r["label"], r["therm_move_mm"], r["surge_mm"], r["react_peak_kN"],
                     r["react_static_kN"], r["stress"]["vm_peak"] / 1e6])
    with open("stage4_comparison.csv", "w", newline="") as fh:     # save the table for your report
        w_ = csv.writer(fh)
        w_.writerow(["case", "thermal_move_mm", "peak_surge_mm", "anchor_peak_kN", "anchor_static_kN", "VM_peak_MPa"])
        w_.writerows(rows)

    # cross-check: the two independent time integrators should agree
    rd = run_case(base, "baseline-direct", method="direct", verbose=False)
    print(f"\nModal vs direct Newmark peak surge: {results[0]['surge_mm']:.3f} vs {rd['surge_mm']:.3f} mm")

    # hydraulic benchmark: solver pressure rise vs Joukowsky, and a vapour-pressure warning
    h = results[0]["hyd"]
    pv = h.P[:, -1]                                      # pressure history at the valve
    print(f"\nValve rise (MOC) = {(pv.max()-base.p_res)/1e6:.4f} MPa vs rho*c*du = {results[0]['dp_joukowsky']/1e6:.4f} MPa")
    print(f"Min absolute pressure at valve = {pv.min()/1e6:.3f} MPa  (<~0.1 MPa means LNG would vapourise: linear theory invalid)")

    # figures for your slides
    plot_pressure(results[0], "fig_pressure.png")
    plot_modes(results[0], "fig_modes.png")
    plot_spectra(results[0], "fig_spectra.png")
    plot_response(results, "fig_response.png")

    # convergence / sensitivity evidence required by the brief
    for wh in ("damping", "hydraulic_dx", "structural_mesh", "time_step"):
        sensitivity(base, wh)
