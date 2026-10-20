#!/usr/bin/env python3
"""Run the S04 fractional-cache cyclicAMI field and GAMG integration screen."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import sys
import uuid

sys.dont_write_bytecode = True
import runS03ApplicationScreens as screens
import runTests


TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]
B0_PROJECT_DIR = Path("/usr/lib/openfoam/openfoam2606")


def load_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-id", required=True, choices=("C0",))
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--library-root", required=True, type=Path)
    parser.add_argument("--solid-foam", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=300)
    return parser.parse_args()


def main() -> int:
    args = load_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError("--timeout must be finite and positive")

    project_dir = args.project_dir.resolve(strict=True)
    library_root = args.library_root.resolve(strict=True)
    solid_foam = args.solid_foam.resolve(strict=True)
    run_root = args.run_root.resolve()
    results_dir = args.results_dir.resolve()
    report_path = args.report.resolve()
    protected_roots = {REPO_ROOT.resolve(), B0_PROJECT_DIR.resolve(strict=True)}
    for path in (run_root, results_dir, report_path):
        if any(path.is_relative_to(root) for root in protected_roots):
            raise ValueError(
                f"S04 artifacts cannot be inside an OpenFOAM tree: {path}"
            )

    env = os.environ.copy()
    if Path(env.get("FOAM_LIBBIN", "")).resolve() != library_root:
        raise ValueError("FOAM_LIBBIN does not match --library-root")
    if not solid_foam.is_file() or not library_root.is_dir():
        raise ValueError("S04 library root or solidFoam is missing")
    env["PATH"] = os.pathsep.join(
        (
        str(solid_foam.parent),
        env.get("PATH", ""),
        )
    )
    target = screens.target_identity(
        args.target_id,
        project_dir,
        env,
        external_tool_roots=(solid_foam.parent,),
    )
    run_root.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite {report_path}")
    identity_dir = run_root / "target-identity"
    identity_dir.mkdir(parents=True, exist_ok=False)
    solid_identity = runTests.inspect_target(
        args.target_id,
        project_dir,
        Path(target["tools"]["solidFoam"]["path"]),
        env["WM_OPTIONS"],
        env,
        identity_dir,
        allowed_library_roots=(library_root,),
    )

    variants_path = (
        TEST_DIR / "cases" / "solidFoam" / "variants.json"
    )
    variants = json.loads(variants_path.read_text(encoding="utf-8"))
    values = variants["T18"].copy()
    values.update(
        {
            "endTime": 0.002,
            "deltaT": 0.001,
            "hotT": 320,
            "coldT": 280,
            "motionFunction": values["cacheMotionFunction"],
            "motionAmplitude": [0, 0, 5],
            "motionOmega": math.pi/(2*0.0013),
            "cacheSize": 360,
            "nThetaStencilMax": 8,
            "requireGAMG": True,
            "requireFractionalCache": True,
        }
    )

    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id += "-" + uuid.uuid4().hex[:8]
    case_dir = run_root / f"T18-fractional-cache-{run_id}" / "case"
    case_dir.parent.mkdir(parents=True, exist_ok=False)
    record = screens.run_solid_case(
        "T18",
        values,
        case_dir,
        env,
        args.timeout,
        args.target_id,
        dynamic_ami=True,
        dynamic_topology=False,
        cache_ami=True,
        record_id="T18-fractional-cache-field-gamg",
    )

    fractional_records = record["fractional_cache_records"]
    fractional_records_valid = bool(fractional_records) and all(
        item["index0"] >= 0
        and item["index1"] >= 0
        and 0 < item["weight"] < 1
        for item in fractional_records
    )
    field_metrics_valid = (
        isinstance(record["field_metrics"], dict)
        and "error" not in record["field_metrics"]
    )
    checks = {
        "runtime_completed": record["status"] == "runtime_pass",
        "pure_zone_mixture_selected": record["runtime_mixture_selected"],
        "cache_reuse_verified": record["cache_reuse_verified"],
        "fractional_cache_verified": (
            record["cache_fraction_verified"]
            and fractional_records_valid
        ),
        "fraction_used_by_gamg": (
            record["cache_fraction_gamg_verified"]
            and record["gamg_solver_events"] > 0
        ),
        "field_metrics_available": field_metrics_valid,
    }
    status = "pass" if all(checks.values()) else "failure"
    report = {
        "schema_version": 1,
        "stage": "S04",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_id": run_id,
        "status": status,
        "target": target,
        "solid_foam_identity": solid_identity,
        "case_parameters": {
            "end_time": values["endTime"],
            "delta_t": values["deltaT"],
            "hot_temperature": values["hotT"],
            "cold_temperature": values["coldT"],
            "motion_function": values["motionFunction"],
            "motion_amplitude": values["motionAmplitude"],
            "motion_omega": values["motionOmega"],
            "cache_size": values["cacheSize"],
            "n_theta_stencil_max": values["nThetaStencilMax"],
            "library_root": str(library_root),
            "solid_foam": str(solid_foam),
            "cache_owner_patch": "AMI2",
            "dynamic_topology": False,
            "material_model": "pureZoneMixture",
        },
        "checks": checks,
        "record": record,
        "input_manifest": runTests.input_manifest_files(),
        "report_path": str(report_path),
    }
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"S04 integration report: {report_path}")
    print(f"S04 integration status: {status}")
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"S04 cache integration error: {error}", file=sys.stderr)
        raise SystemExit(2)
