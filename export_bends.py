"""Export four real hydraulic stations after the group supplies a route JSON."""

import argparse
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import numpy as np
from complete_stage1 import provenance, shared_case
from fluid_assumptions import PHASE_BASIS, friction_estimate
from stage1 import simulate, theory


def export_route(config, output, code_version="uncommitted; use source_sha256"):
    length, bends = config.get("length_m"), config.get("bend_stations_m")
    if (not isinstance(length, (int, float)) or isinstance(length, bool)
            or not math.isfinite(length) or length <= 0
            or not isinstance(bends, list) or len(bends) != 4):
        raise ValueError("Supply a positive total hydraulic length and four real bend stations.")
    if (any(not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x)
            for x in bends) or not 0 < bends[0] < bends[1] < bends[2] < bends[3] < length):
        raise ValueError("Four bend stations must increase strictly between reservoir and valve.")
    if config.get("speed_mode") not in ("bulk", "elastic") or not config.get("geometry_source"):
        raise ValueError("Supply speed_mode and a traceable geometry_source.")
    base = shared_case(config["speed_mode"], 0)
    friction = friction_estimate(base.density, base.velocity, base.diameter)["darcy_factor"]
    case = replace(base, friction=friction, length=length)
    if "phase_limit_pa" in config:
        if not config.get("phase_basis"):
            raise ValueError("An overridden phase limit needs phase_basis.")
        limit = config["phase_limit_pa"]
        if not isinstance(limit, (int, float)) or isinstance(limit, bool) or not math.isfinite(limit) or limit <= 0:
            raise ValueError("An overridden phase limit must be positive absolute Pa.")
        case = replace(case, phase_limit=limit)
    result = simulate(case, config["reaches"], config["duration_s"], stations=bends)
    output.mkdir(parents=True, exist_ok=True)
    crossing = result["first_phase_crossing_s"]
    pre_crossing = np.ones(len(result["time"])) if crossing is None else (result["time"] < crossing)
    headings = (["time_s"]+[f"p_abs_bend{i}_pa" for i in range(1, 5)]
                +[f"Q_bend{i}_m3_s" for i in range(1, 5)]+["before_global_phase_crossing"])
    np.savetxt(output/"bend_histories.csv",
               np.column_stack((result["time"], result["pressure"], result["flow"], pre_crossing)),
               delimiter=",", fmt="%.12g", header=",".join(headings), comments="")
    metadata = {**provenance(code_version), "route": config, "parameters_SI": asdict(case),
                "phase_basis": config.get("phase_basis", PHASE_BASIS),
                "actual_LNG_composition_known": False, "analytical": theory(case),
                "pressure_convention": "Absolute Pa; structural pressure = p_abs - external_abs",
                "flow_convention": "Signed volume flow, m3/s; positive reservoir to valve",
                "closure_law": "Q_valve(t) = A*U0*max(1-t/T_close,0)",
                "station_definition": "Arc length from reservoir; compact-bend approximation",
                "hydraulic_model": "Uniform-area horizontal pipe; no bend losses, reflections or structural feedback",
                "dt_s": result["dt"], "dx_m": result["dx"],
                "interpolated_stations": result["interpolated_stations"],
                "minimum_absolute_pa": result["minimum_absolute_pa"],
                "first_phase_crossing_s": crossing, "phase_status": result["phase_status"],
                "approved_physical_load": False,
                "interpretation": "Surrogate/supplied-threshold screening only; no cavitation model"}
    (output/"metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--out", type=Path, default=Path("results/latest/bends"))
    parser.add_argument("--code-version", default="uncommitted; use source_sha256")
    args = parser.parse_args()
    export_route(json.loads(args.config.read_text()), args.out, args.code_version)
