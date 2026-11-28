"""Standalone Stage 1 column-separation sensitivity. No loop loads are exported."""
from dataclasses import replace
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import fluid_assumptions
import stage1
from fluid_assumptions import friction_estimate, methane_saturation_pa
from stage1 import Case, G, simulate, theory

STATIONS = (10., 25., 50., 100.)
GRIDS = (80, 160, 320, 640, 1280, 2560)
DURATION = 1.
REFERENCE = "be7f68ab1348566e677c8b58c22f372f017e7f39"

def cavity_run(case, reaches, phase_pressure, weight=.5, duration=DURATION):
    """Discrete vapour cavities at MOC nodes, with separate face discharges.

    Integrate dV/dt = Q_down - Q_up using weighted adjacent time samples.
    A collapsed cavity rejoins at the current grid time. Reopening within
    the same time step cannot be resolved by this scheme.
    """
    dx = case.length/reaches
    dt = dx/case.speed
    area = case.area
    impedance = case.speed/(G*area)
    resistance = case.friction*dx/(2*G*case.diameter*area**2)
    reservoir = (case.pressure_abs - case.external_abs)/(case.density*G)
    vapour = (phase_pressure - case.external_abs)/(case.density*G)
    nodes = dx*np.arange(reaches + 1)
    time = dt*np.arange(math.floor(duration/dt + 1e-9) + 1)
    head = reservoir - case.friction*nodes*case.velocity**2/(2*G*case.diameter)
    if np.min(head) <= vapour:
        raise ValueError("Initial liquid pressure must exceed the selected phase threshold.")
    initial_q = area*case.velocity
    incoming_q = np.full(reaches+1, initial_q)
    outgoing_q = incoming_q.copy()
    maximum_balance_error = 0.
    volume = np.zeros(reaches+1)
    gaps = np.zeros_like(volume)
    pressure = np.empty((len(time),len(STATIONS)))
    valve_volume = np.empty(len(time))
    maximum_ratio = 0.
    first_activation = None
    valve_activation = None
    valve_collapse = None
    interior_activation = None
    largest_pressure = case.pressure_abs
    reopens = 0
    for k,t in enumerate(time):
        absolute = case.external_abs + case.density*G*head
        pressure[k] = np.interp(STATIONS, nodes, absolute)
        valve_volume[k] = volume[-1]
        maximum_ratio = max(maximum_ratio, float(volume.max())/(area*dx))
        largest_pressure = max(largest_pressure,float(absolute.max()))
        if k == len(time)-1:
            break
        following = time[k+1]
        # A positive characteristic leaves the downstream face of its old node.
        plus = head[:-1] + impedance*outgoing_q[:-1] - resistance*outgoing_q[:-1]*np.abs(outgoing_q[:-1])
        minus = head[1:] - impedance*incoming_q[1:] + resistance*incoming_q[1:]*np.abs(incoming_q[1:])
        valve_q = initial_q*max(0., 1.-following/case.closure_time)
        head[0] = reservoir
        head[1:-1] = (plus[:-1]+minus[1:])/2
        head[-1] = plus[-1]-impedance*valve_q

        upstream = (plus-vapour)/impedance
        downstream = np.append((vapour-minus[1:])/impedance,valve_q)
        new_gaps = downstream-upstream
        previously_active = volume[1:] > 0.
        liquid_below = head[1:] < vapour
        old_gap = np.where(previously_active,gaps[1:],0.)
        proposed = volume[1:] + dt*((1.-weight)*old_gap+weight*new_gaps)
        unresolved_reopening = liquid_below & (proposed <= 0.)
        reopens += int(np.count_nonzero(unresolved_reopening))
        proposed[unresolved_reopening] = dt*weight*new_gaps[unresolved_reopening]
        active = (previously_active | liquid_below) & (proposed > 0.)
        if active.any() and first_activation is None:
            first_activation = float(following)
        if active[:-1].any() and interior_activation is None:
            interior_activation = float(following)
        if active[-1] and valve_activation is None:
            valve_activation = float(following)
        if previously_active[-1] and not active[-1] and valve_collapse is None:
            valve_collapse = float(following)

        # Check continuity at nodes with a cavity at both time levels.
        persistent = previously_active & active & ~unresolved_reopening
        if persistent.any():
            implied = (proposed[persistent]-volume[1:][persistent])/dt
            target = (1.-weight)*old_gap[persistent] + weight*new_gaps[persistent]
            maximum_balance_error = max(maximum_balance_error,float(np.max(np.abs(implied-target))))
        volume[1:] = np.where(active,proposed,0.)
        gaps[1:] = np.where(active,new_gaps,0.)
        incoming_q[0] = outgoing_q[0] = (reservoir-minus[0])/impedance
        shared = (plus[:-1]-minus[1:])/(2*impedance)
        incoming_q[1:-1] = np.where(active[:-1],upstream[:-1],shared)
        outgoing_q[1:-1] = np.where(active[:-1],downstream[:-1],shared)
        outgoing_q[-1] = valve_q
        incoming_q[-1] = upstream[-1] if active[-1] else valve_q
        head[1:][active] = vapour
    return dict(time=time,pressure=pressure,valve_volume=valve_volume,
                max_cavity_to_reach_volume=maximum_ratio,first_activation_s=first_activation,
                valve_activation_s=valve_activation,interior_activation_s=interior_activation,
                valve_first_collapse_s=valve_collapse,largest_pressure_pa=largest_pressure,
                substep_reopening_count=reopens,dt_s=dt,max_balance_error_m3_s=maximum_balance_error)

def compare(case, baseline, run, grid, phase, weight, label):
    time = run["time"]
    old = baseline["pressure"][:len(time)]
    np.testing.assert_allclose(time, baseline["time"][:len(time)], rtol=0, atol=1e-13)
    pressure = run["pressure"]
    active_start = run["first_activation_s"]
    prior = time < active_start if active_start is not None else np.ones(len(time),dtype=bool)
    pre_error = np.max(np.abs(pressure[prior]-old[prior]))
    if pre_error > 1e-5:
        raise RuntimeError(f"Pre-activation mismatch: {pre_error:g} Pa")
    if (np.min(pressure) < phase-1e-6 or np.min(run["valve_volume"]) < 0.
            or run["max_balance_error_m3_s"] > 1e-9):
        raise RuntimeError("Cavity pressure bound or cavity-volume check failed.")
    crest = time < theory(case)["round_trip_s"]
    valve = pressure[:, -1]
    baseline_valve = old[:, -1]
    first_collapse = run["valve_first_collapse_s"]
    post = time >= first_collapse if first_collapse is not None else np.zeros(len(time),dtype=bool)
    return dict(case=label,reaches=grid,dt_ms=run["dt_s"]*1e3,phase_pressure_mpa=phase/1e6,
                volume_weight=weight,baseline_first_rise_mpa=np.max(baseline_valve[crest]-baseline_valve[0])/1e6,
                cavity_first_rise_mpa=np.max(valve[crest]-valve[0])/1e6,
                baseline_min_valve_mpa=np.min(baseline_valve)/1e6,cavity_min_valve_mpa=np.min(valve)/1e6,
                first_cavity_ms=None if active_start is None else active_start*1e3,
                valve_first_collapse_ms=None if first_collapse is None else first_collapse*1e3,
                interior_first_cavity_ms=None if run["interior_activation_s"] is None else run["interior_activation_s"]*1e3,
                peak_valve_from_first_collapse_mpa=np.max(valve[post])/1e6 if post.any() else None,
                peak_anywhere_mpa=run["largest_pressure_pa"]/1e6,max_valve_cavity_l=np.max(run["valve_volume"])*1e3,
                max_cavity_to_reach_volume_percent=run["max_cavity_to_reach_volume"]*100,
                pre_activation_difference_pa=pre_error,substep_reopenings=run["substep_reopening_count"])

def plot(output, phase, base, cases, rows):
    colour = "#432260"
    fig, (pressure, volume) = plt.subplots(2,1,figsize=(9.2,6.1),sharex=True,layout="constrained",
                                          gridspec_kw={"height_ratios":[2,1]})
    t = cases[640]["time"]
    pressure.plot(t,base["pressure"][:len(t),-1]/1e6,"--",color="#64616a",lw=1.4,label="Original single-phase model")
    pressure.plot(t,cases[640]["pressure"][:,-1]/1e6,color=colour,lw=1.6,label="Cavity sensitivity, N = 640")
    pressure.axhline(phase/1e6,color="#ac6f09",ls=":",lw=1.5,label="Pure-methane proxy at 111 K")
    pressure.axvspan(cases[640]["valve_first_collapse_s"],1.,color="#bdb5c3",alpha=.12,lw=0,
                     label="After first modelled collapse")
    pressure.set(ylabel="Valve pressure [MPa abs.]")
    pressure.legend(loc="upper left",fontsize=8,frameon=True)
    for grid,colour in [(320,"#578697"),(640,"#432260"),(1280,"#c87821")]:
        run=cases[grid]
        volume.plot(run["time"],run["valve_volume"]*1e3,lw=1.2,label=f"N = {grid}",color=colour)
    volume.set(xlabel="Time after valve closure starts [s]",ylabel="Valve cavity volume [L]",xlim=(0,1))
    volume.legend(loc="upper left",fontsize=8)
    for axis in (pressure,volume):
        axis.axvline(cases[640]["valve_activation_s"],color="#ac6f09",lw=1.,alpha=.5)
        axis.grid(alpha=.22)
    fig.suptitle("Cavitation sensitivity of the 100 m benchmark")
    fig.savefig(output/"cavitation_valve_comparison.png",dpi=220)
    plt.close(fig)

    fig, (timing, rebound) = plt.subplots(1,2,figsize=(9.2,3.5),layout="constrained")
    resolution=[r for r in rows if r["case"]=="methane_proxy"]
    ns=[r["reaches"] for r in resolution]
    timing.plot(ns,[r["first_cavity_ms"] for r in resolution],"o-",label="First cavity")
    timing.plot(ns,[r["valve_first_collapse_ms"] for r in resolution],"s-",label="First valve collapse")
    timing.set(ylabel="Time [ms]")
    timing.legend(fontsize=8)
    rebound.plot(ns,[r["peak_valve_from_first_collapse_mpa"] for r in resolution],"o-",color="#7b459d")
    rebound.set(ylabel="Later valve peak [MPa abs.]")
    for axis in (timing,rebound):
        axis.set(xlabel="Number of reaches")
        axis.set_xscale("log",base=2)
        axis.set_xticks(ns,labels=list(map(str,ns)))
        axis.grid(alpha=.22)
    fig.suptitle("Sensitivity of cavity predictions to mesh refinement")
    fig.savefig(output/"cavitation_mesh_sensitivity.png",dpi=220)
    plt.close(fig)

def run(output):
    output.mkdir(parents=True,exist_ok=True)
    case = Case(speed_mode="bulk",friction=friction_estimate()["darcy_factor"])
    phase = methane_saturation_pa(111.)
    rows, cases = [], {}
    for grid in GRIDS:
        baseline = simulate(case,grid,DURATION)
        result = cavity_run(case,grid,phase)
        cases[grid] = result
        rows.append(compare(case,baseline,result,grid,phase,.5,"methane_proxy"))
    for grid in (320,640,1280,2560):
        baseline = simulate(case,grid,DURATION)
        for limit,weight,label in [(100000.,.5,"rounded_100_kpa"),(phase,1.,"implicit_weight")]:
            result = cavity_run(case,grid,limit,weight)
            rows.append(compare(case,baseline,result,grid,limit,weight,label))
    no_cavity_case = replace(case,closure_time=.25)
    slow = simulate(no_cavity_case,320,DURATION)
    check = cavity_run(no_cavity_case,320,phase)
    if check["first_activation_s"] is not None:
        raise RuntimeError("Unexpected cavity in slow-closure control.")
    np.testing.assert_allclose(check["pressure"],slow["pressure"][:len(check["time"])],rtol=0,atol=1e-5)
    ideal = replace(case,friction=0.)
    early = cavity_run(ideal,640,phase,duration=.22)
    # The first frictionless valve reflection ramps down by 2*Joukowsky rise.
    target = theory(ideal)
    onset = target["round_trip_s"] + ideal.closure_time/2*(1. + (ideal.pressure_abs-phase)/target["surge_pa"])
    detected = early["valve_activation_s"]
    if not 0. <= detected-onset <= early["dt_s"]:
        raise RuntimeError("Frictionless valve activation differs from the analytical crossing.")
    full = simulate(case,640,DURATION)
    with (output/"cavitation_summary.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savetxt(output/"valve_history_N640.csv",
               np.column_stack((cases[640]["time"],full["pressure"][:len(cases[640]["time"]),-1],
                                cases[640]["pressure"][:,-1],cases[640]["valve_volume"])),
               delimiter=",",fmt="%.12g",comments="",
               header="time_s,original_single_phase_valve_pressure_pa_abs,dvcm_sensitivity_valve_pressure_pa_abs,valve_cavity_volume_m3")
    metadata=dict(source_commit=REFERENCE,python=platform.python_version(),numpy=np.__version__,
                  matplotlib=matplotlib.__version__,phase_pressure_pa=phase,
                  max_pre_activation_difference_pa=max(row["pre_activation_difference_pa"] for row in rows),
                  slow_closure_control_max_pressure_difference_pa=float(np.max(np.abs(
                      check["pressure"]-slow["pressure"][:len(check["time"])]))),
                  frictionless_analytical_activation_s=onset,
                  frictionless_detected_activation_s=detected,
                  methane_runs_max_cavity_balance_residual_m3_s=max(
                      item["max_balance_error_m3_s"] for item in cases.values()),
                  single_phase_inputs=case.__dict__,record_s=DURATION,
                  source_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in map(Path, (stage1.__file__, fluid_assumptions.__file__, __file__))})
    (output/"run_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    plot(output,phase,full,cases,rows)
    for r in rows:
        print(r["case"],r["reaches"],*[f'{r[k]:.5g}' if r[k] is not None else '--' for k in
              ("first_cavity_ms","valve_first_collapse_ms","peak_valve_from_first_collapse_mpa",
               "max_cavity_to_reach_volume_percent","max_valve_cavity_l","substep_reopenings")])

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("studies/stage1_cavitation/results"))
    run(parser.parse_args().out)
