"""Literature-based methods with explicit, provisional coursework assumptions.

See docs/FRICTION_AND_PHASE.md: no actual LNG composition or pipe finish is known.
"""

import math

TEMPERATURE_K = 111.0
MU_PA_S = 1.20e-4
MU_BOUNDS_PA_S = (0.96e-4, 1.44e-4)
ROUGHNESS_M = 15e-6
ROUGHNESS_BOUNDS_M = (0.0, 50e-6)
PHASE_BASIS = "Pure methane surrogate at 111 K; actual LNG composition unknown"


def methane_saturation_pa(temperature=TEMPERATURE_K):
    """Setzmann & Wagner (1991), Eq. 3.2, Table 26; absolute pressure."""
    if not math.isfinite(temperature) or not 90.6941 <= temperature < 190.564:
        raise ValueError("Methane correlation requires 90.6941 <= T < 190.564 K.")
    theta = 1-temperature/190.564
    series = (-6.036219*theta + 1.409353*theta**1.5
              - 0.4945199*theta**2 - 1.443048*theta**4.5)
    return 4.5992e6*math.exp(190.564/temperature*series)


def colebrook(reynolds, relative_roughness):
    """Darcy factor: solve the turbulent Colebrook equation by bisection."""
    if (not math.isfinite(reynolds) or reynolds < 4000
            or not math.isfinite(relative_roughness)
            or not 0 <= relative_roughness < 0.05):
        raise ValueError("Require turbulent Re >= 4000 and 0 <= epsilon/D < 0.05.")
    lower, upper = 1.0, 30.0  # y = 1/sqrt(f_D)
    for _ in range(80):
        y = (lower+upper)/2
        residual = y+2*math.log10(relative_roughness/3.7+2.51*y/reynolds)
        if residual > 0:
            upper = y
        else:
            lower = y
    return ((lower+upper)/2)**-2


def friction_estimate(density=455.0, velocity=1.5, diameter=0.79394):
    """Bounds propagate the stated scenario envelope, not a confidence interval."""
    if any(not math.isfinite(v) or v <= 0 for v in (density, velocity, diameter)):
        raise ValueError("Positive density, velocity and diameter are required.")
    numerator = density*velocity*diameter
    corners = [{"mu_pa_s": mu, "roughness_m": rough,
                "reynolds": numerator/mu,
                "darcy_factor": colebrook(numerator/mu, rough/diameter)}
               for mu in MU_BOUNDS_PA_S for rough in ROUGHNESS_BOUNDS_M]
    return {"mu_pa_s": MU_PA_S, "roughness_m": ROUGHNESS_M,
            "reynolds": numerator/MU_PA_S,
            "darcy_factor": colebrook(numerator/MU_PA_S, ROUGHNESS_M/diameter),
            "lower": min(c["darcy_factor"] for c in corners),
            "upper": max(c["darcy_factor"] for c in corners), "corners": corners,
            "bounds_basis": "Assumed viscosity +/-20%, equivalent roughness 0-50 micrometres"}
