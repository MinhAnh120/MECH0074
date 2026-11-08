"""Head/flow MOC adapted from the supplied MECH0074 notebook, cells 6-7.

SI units throughout. This is a single-phase teaching model, not a design approval.
"""

from dataclasses import dataclass
import math
import numpy as np

G = 9.81


@dataclass(frozen=True)
class Case:
    length: float = 100.0
    diameter: float = 0.79394
    outer_diameter: float = 0.813
    density: float = 455.0
    bulk_speed: float = 1100.0
    young: float = 213.233e9
    velocity: float = 1.5
    pressure_abs: float = 700000.0  # Assigned to the reservoir boundary.
    external_abs: float = 101325.0
    closure_time: float = 0.020
    friction: float = 0.0  # Darcy factor; 0.020 is only a teaching sensitivity.
    speed_mode: str = "elastic"
    phase_limit: float | None = None

    def __post_init__(self):
        positive = (self.length, self.diameter, self.outer_diameter, self.density,
                    self.bulk_speed, self.young, self.pressure_abs,
                    self.external_abs, self.closure_time)
        if any(not math.isfinite(v) or v <= 0 for v in positive):
            raise ValueError("Dimensions, properties, pressures and closure must be positive.")
        if self.outer_diameter <= self.diameter:
            raise ValueError("Outer diameter must exceed inner diameter.")
        if any(not math.isfinite(v) or v < 0 for v in (self.velocity, self.friction)):
            raise ValueError("Initial velocity and Darcy factor must be nonnegative.")
        if self.speed_mode not in ("elastic", "bulk"):
            raise ValueError("speed_mode must be 'elastic' or 'bulk'.")
        if self.phase_limit is not None:
            if not math.isfinite(self.phase_limit) or self.phase_limit <= 0:
                raise ValueError("An entered phase limit must be positive absolute pressure.")

    @property
    def area(self):
        return math.pi * self.diameter**2 / 4

    @property
    def speed(self):
        if self.speed_mode == "bulk":
            return self.bulk_speed
        wall = (self.outer_diameter - self.diameter) / 2
        return (self.bulk_speed**-2 + self.density*self.diameter/(wall*self.young))**-0.5


def theory(case):
    a, length = case.speed, case.length
    return {"speed_m_s": a, "area_m2": case.area,
            "initial_flow_m3_s": case.area*case.velocity,
            "surge_pa": case.density*a*case.velocity,
            "one_way_s": length/a, "round_trip_s": 2*length/a,
            "period_s": 4*length/a, "fundamental_hz": a/(4*length),
            "closure_ratio": case.closure_time/(2*length/a)}


def simulate(case, reaches, duration, stations=(10.0, 25.0, 50.0, 100.0),
             hold_open=False):
    """Keep current pipe fields and requested histories; inspect all nodes for minima.

    Off-grid stations are explicitly sampled by linear interpolation. hold_open
    suppresses closure solely for the steady-state verification case.
    """
    if isinstance(reaches, bool) or not isinstance(reaches, (int, np.integer)) or reaches < 1:
        raise ValueError("reaches must be a positive integer.")
    if not math.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be positive and finite.")
    stations = np.asarray(stations, dtype=float)
    if (stations.ndim != 1 or stations.size == 0 or not np.isfinite(stations).all()
            or np.any(stations < 0) or np.any(stations > case.length)
            or np.any(np.diff(stations) <= 0)):
        raise ValueError("Stations must be ordered, distinct and within the pipe.")
    a, area = case.speed, case.area
    dx = case.length/reaches
    dt = dx/a
    steps = math.ceil(duration/dt)
    time = np.arange(steps + 1)*dt
    nodes = np.linspace(0, case.length, reaches + 1)
    q0 = area*case.velocity
    reservoir_head = (case.pressure_abs-case.external_abs)/(case.density*G)
    h = reservoir_head-case.friction*nodes*case.velocity**2/(2*G*case.diameter)
    q = np.full(reaches + 1, q0)
    b = a/(G*area)
    r = case.friction*dx/(2*G*case.diameter*area**2)
    position = stations/dx
    left = np.minimum(np.floor(position).astype(int), reaches-1)
    weight = position-left
    interpolated = bool(np.any(np.minimum(abs(weight), abs(1-weight)) > 1e-10))
    pressure = np.empty((steps + 1, stations.size))
    flow = np.empty_like(pressure)
    minimum = math.inf
    first_negative = first_phase = None

    for n, instant in enumerate(time):
        pressure[n] = case.external_abs+case.density*G*((1-weight)*h[left]+weight*h[left+1])
        flow[n] = (1-weight)*q[left]+weight*q[left+1]
        local_min = float(case.external_abs+case.density*G*h.min())
        minimum = min(minimum, local_min)
        if local_min <= 0 and first_negative is None:
            first_negative = float(instant)
        if case.phase_limit is not None and local_min <= case.phase_limit and first_phase is None:
            first_phase = float(instant)
        if n == steps:
            break

        # All characteristics come from the previous time level.
        cplus = h[:-2]+b*q[:-2]-r*q[:-2]*np.abs(q[:-2])
        cminus = h[2:]-b*q[2:]+r*q[2:]*np.abs(q[2:])
        hn, qn = np.empty_like(h), np.empty_like(q)
        hn[1:-1] = (cplus+cminus)/2
        qn[1:-1] = (cplus-cminus)/(2*b)
        incoming = h[1]-b*q[1]+r*q[1]*abs(q[1])
        hn[0] = reservoir_head
        qn[0] = (reservoir_head-incoming)/b
        incoming = h[-2]+b*q[-2]-r*q[-2]*abs(q[-2])
        fraction = 1.0 if hold_open else max(1-time[n+1]/case.closure_time, 0.0)
        qn[-1] = q0*fraction
        hn[-1] = incoming-b*qn[-1]
        h, q = hn, qn

    if first_negative is not None:
        status = "nonpositive_absolute_pressure"
    elif case.phase_limit is None:
        status = "LNG_phase_limit_not_supplied"
    elif first_phase is not None:
        status = "phase_limit_crossed"
    else:
        status = "above_entered_phase_limit_only"
    return {"time": time, "stations": stations, "pressure": pressure, "flow": flow,
            "dx": dx, "dt": dt, "interpolated_stations": interpolated,
            "minimum_absolute_pa": minimum, "first_nonpositive_s": first_negative,
            "first_phase_crossing_s": first_phase, "phase_status": status}


def estimate_period(time, pressure, start):
    """Median interval between equivalent upward zero crossings after startup."""
    indices = np.flatnonzero((pressure[:-1] <= 0) & (pressure[1:] > 0))
    crossings = time[indices]-pressure[indices]*np.diff(time)[indices]/np.diff(pressure)[indices]
    crossings = crossings[crossings > start]
    return float(np.median(np.diff(crossings))) if crossings.size >= 2 else None


def spectrum(time, pressure):
    """Mean-removed, Hann-windowed single-sided amplitude spectrum in Pa."""
    time, pressure = np.asarray(time), np.asarray(pressure)
    if (time.ndim != 1 or time.size < 4 or pressure.shape != time.shape
            or not np.isfinite(time).all() or not np.isfinite(pressure).all()):
        raise ValueError("Spectrum requires matching finite one-dimensional records.")
    dt = float(time[1]-time[0])
    if dt <= 0 or not np.allclose(np.diff(time), dt, rtol=1e-9, atol=1e-12):
        raise ValueError("FFT samples must be uniformly spaced in increasing time.")
    signal = pressure-pressure.mean()
    if np.max(abs(signal)) < 1e-10:
        raise ValueError("No oscillation is available for frequency estimation.")
    window = np.hanning(time.size)
    amplitude = 2*np.abs(np.fft.rfft(signal*window))/window.sum()
    amplitude[0] /= 2
    if time.size % 2 == 0:
        amplitude[-1] /= 2
    frequency = np.fft.rfftfreq(time.size, dt)
    peak = int(np.argmax(amplitude[1:])+1)
    return frequency, amplitude, float(frequency[peak])
