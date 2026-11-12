#!/usr/bin/env python3
"""Build and run the cyclicAMI / faceAreaWeightAMI regression suite."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import uuid

sys.dont_write_bytecode = True
from solidMetrics import SOLID_ACCEPTANCE_LIMITS, solid_reference_metrics
from t03Variants import write_case as write_t03_variant


TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]
B0_PROJECT_DIR = Path("/usr/lib/openfoam/openfoam2606")
ALL_TEST_IDS = [f"T{i:02d}" for i in range(20)]

TEST_SPECS = {
    "T01": {
        "name": "planar transfer control",
        "variant": "planar-nonmatching-1-to-2",
        "case": "analytical",
        "arguments": ["-test", "T01"],
        "assertions": {
            "T01_CALCULATED",
            "T01_SOURCE_COVERAGE",
            "T01_CONSTANT_SOURCE",
            "T01_TARGET_COVERAGE",
            "T01_CONSTANT_TARGET",
            "T01_NONNEGATIVE_WEIGHTS",
            "T01_WEIGHTED_ADJOINT",
        },
    },
    "T02": {
        "name": "cached interpolation fraction",
        "variant": "cache-midpoint-fraction",
        "case": "analytical",
        "arguments": ["-test", "T02"],
        "assertions": {
            "T02_CACHE_RESTORED",
            "T02_CACHE_WEIGHT_TYPE",
            "T02_CACHE_FRACTION",
            "T02_CACHE_MIDPOINT",
        },
    },
    "T03": {
        "name": "rotational area vector after topology update",
        "variant": "rotation-90deg-topology-update",
        "case": "topology",
        "prepare_mesh": True,
        "arguments": ["-test", "T03"],
        "assertions": {
            "T03_AMI_PATCHES_PRESENT",
            "T03_OWNER_FACE_EXISTS",
            "T03_NEIGHBOUR_FACE_EXISTS",
            "T03_OWNER_AREA_GEOMETRY",
            "T03_OWNER_AREA_FRAME",
            "T03_NEIGHBOUR_AREA_GEOMETRY",
            "T03_NEIGHBOUR_AREA_FRAME",
            "T03_NEIGHBOUR_MAGNITUDE",
            "T03_NEIGHBOUR_CENTRE_GEOMETRY",
            "T03_CELL_VOLUMES_POSITIVE",
            "T03_CELL_CLOSURE",
            "T03_AMI_DISTRIBUTED",
        },
    },
    "T04": {
        "name": "rotational low-weight default frame",
        "variant": "low-weight-rotation",
        "case": "rotation",
        "prepare_mesh": True,
        "arguments": ["-test", "T04"],
        "assertions": {
            "T04_AMI_PATCHES_PRESENT",
            "T04_FACE_FRAME_U_AMI1",
            "T04_FACE_FRAME_U_AMI2",
            "T04_FACES_CHECKED_U_AMI1",
            "T04_FACES_CHECKED_U_AMI2",
            "T04_FACE_FRAME_S_AMI1",
            "T04_FACE_FRAME_S_AMI2",
            "T04_FACES_CHECKED_S_AMI1",
            "T04_FACES_CHECKED_S_AMI2",
            "T04_FACE_FRAME_T_AMI1",
            "T04_FACE_FRAME_T_AMI2",
            "T04_FACES_CHECKED_T_AMI1",
            "T04_FACES_CHECKED_T_AMI2",
            "T04_FACE_FRAME_ST_AMI1",
            "T04_FACE_FRAME_ST_AMI2",
            "T04_FACES_CHECKED_ST_AMI1",
            "T04_FACES_CHECKED_ST_AMI2",
        },
    },
    "T05": {
        "name": "seed selection and face ordering",
        "variant": "empty-first-seed-no-restart",
        "case": "analytical",
        "arguments": ["-test", "T05", "-no-restart"],
        "assertions": {
            "T05_CALCULATED",
            "T05_SOURCE_COVERAGE",
            "T05_TARGET_COVERAGE",
            "T05_PHYSICAL_MAPPING",
        },
    },
    "T06": {
        "name": "fine-face overlap cutoff",
        "variant": "fully-covered-2.5e-7-face",
        "case": "analytical",
        "arguments": ["-test", "T06"],
        "assertions": {
            "T06_CALCULATED",
            "T06_SOURCE_COVERAGE",
            "T06_FINE_FACE_COVERAGE",
        },
        "expected_failures": {
            "T06_SOURCE_COVERAGE",
            "T06_FINE_FACE_COVERAGE",
        },
    },
    "T07": {
        "name": "MPI discovery versus accepted geometric pairs",
        "variant": "two-rank-gap-0.03-max-distance2-0.01",
        "case": "analytical",
        "ranks": 2,
        "arguments": [
            "-parallel",
            "-test",
            "T07",
            "-gap",
            "0.03",
            "-maxDistance2",
            "0.01",
        ],
        # With restart on and requireMatch true, the empty target box on one
        # rank (C0 discovery defect) aborts in faceAreaWeightAMI::calculate()
        # before any T07 assertion is reported. Every rank must print the
        # fatal.
        "assertions": set(),
        "expected_fatal": (
            "1 source faces have no overlapping target face"
            " (requireMatch is true)"
        ),
    },
    "T08": {
        "name": "MPI topology and local communicator matrix",
        "variant": "rank-topology-matrix",
        "case": "analytical",
        "arguments": ["-parallel", "-test", "T08"],
        "assertions": {
            "T08_CALCULATED",
            "T08_SOURCE_COVERAGE",
            "T08_TARGET_COVERAGE",
            "T08_CONSTANT_SOURCE",
            "T08_CONSTANT_TARGET",
            "T08_DISTRIBUTED",
            "T08_SINGLE_PATCH_PROC",
        },
    },
    "T09": {
        "name": "AMI state transitions",
        "variant": "distributed-local-distributed",
        "case": "analytical",
        "arguments": ["-parallel", "-test", "T09", "-state-cycle"],
        "assertions": {
            f"T09_{state}_{assertion}"
            for state in (
                "INITIAL_DISTRIBUTED",
                "LOCAL_AFTER_DISTRIBUTED",
                "DISTRIBUTED_AFTER_LOCAL",
            )
            for assertion in (
                "CALCULATED",
                "DISTRIBUTED",
                "SINGLE_PATCH_PROC",
                "MAPS_PRESENT",
                "SOURCE_COVERAGE",
                "TARGET_COVERAGE",
                "CONSTANT_SOURCE",
                "CONSTANT_TARGET",
            )
        },
    },
    "T10": {
        "name": "partial coverage diagnostics",
        "variant": "source-0.9-target-1.0-and-0.8",
        "case": "analytical",
        "arguments": ["-test", "T10"],
        "assertions": {
            "T10_CALCULATED",
            "T10_SOURCE_COVERAGE",
            "T10_TARGET_COVERAGE_0",
            "T10_TARGET_COVERAGE_1",
        },
    },
    "T11": {
        "name": "Scale and coordinate-offset matrix",
        "variant": "scale-offset",
        "case": "analytical",
        "arguments": ["-test", "T11"],
        "assertions": {
            "T11_CALCULATED",
            "T11_SOURCE_COVERAGE",
            "T11_TARGET_COVERAGE",
            "T11_CONSTANT_SOURCE",
            "T11_CONSTANT_TARGET",
            "T11_SOURCE_AREA_SCALE_RATIO",
            "T11_TARGET_AREA_SCALE_RATIO",
            "T11_SOURCE_OVERLAP_AREA_RATIO",
            "T11_TARGET_OVERLAP_AREA_RATIO",
            "T11_CENTROID_NORMALIZED_ERROR",
        },
    },
    "T12": {
        "name": "empty target input",
        "variant": "empty-target",
        "case": "analytical",
        "arguments": ["-test", "T12"],
        "assertions": {
            "T12_CALCULATED",
            "T12_EMPTY_TARGET_COVERAGE",
        },
    },
    "T13": {
        "name": "piecewise-constant area average",
        "variant": "two-half-target-average",
        "case": "analytical",
        "arguments": ["-test", "T13"],
        "assertions": {
            "T13_CALCULATED",
            "T13_AREA_AVERAGE",
        },
    },
}


T02_S04_QUERY_IDS = (
    "Q000",
    "Q0250",
    "Q0500",
    "Q0750",
    "Q1000",
    "W000",
    "W0500",
    "W1000",
)
T02_S04_QUERY_ASSERTIONS = (
    "RESTORED",
    "INDEX_STATE",
    "WEIGHT",
    "WEIGHT_BOUNDS",
    "RESULT_PRESENT",
    "RESULT",
    "CONSTANT",
    "FINITE",
)
T02_S04_ASSERTIONS = {"T02_CACHE_WEIGHT_TYPE"}
for _direction in ("SOURCE", "TARGET"):
    T02_S04_ASSERTIONS.update(
        {
            f"T02_{_direction}_CALCULATED",
            f"T02_{_direction}_OPERATORS_SET",
            f"T02_{_direction}_ENDPOINT_SEPARATION",
            f"T02_{_direction}_WRAP_OPERATORS_SET",
        }
    )
    T02_S04_ASSERTIONS.update(
        f"T02_{_direction}_{query}_{assertion}"
        for query in T02_S04_QUERY_IDS
        for assertion in T02_S04_QUERY_ASSERTIONS
    )


DEFERRED_TESTS = {}


def topology_scenarios() -> list[dict]:
    scenarios = [
        {
            "variant": "serial",
            "arguments": [
                "-test", "T08",
                "-src-rank", "0",
                "-tgt-rank", "1",
                "-opt-switch", "localAMIComm=1",
            ],
        }
    ]

    layouts = (
        ("two-rank-split", 2, 0, 1),
        ("four-rank-split-empty-ranks", 4, 0, 1),
        ("four-rank-local-pair", 4, 0, 0),
        ("four-rank-noncontiguous-1-3", 4, 1, 3),
    )
    for name, ranks, source_rank, target_rank in layouts:
        for local_comm in (0, 1, 2):
            scenarios.append(
                {
                    "variant": f"{name}-localAMIComm-{local_comm}",
                    "ranks": ranks,
                    "arguments": [
                        "-parallel",
                        "-test", "T08",
                        "-src-rank", str(source_rank),
                        "-tgt-rank", str(target_rank),
                        "-opt-switch", f"localAMIComm={local_comm}",
                    ],
                }
            )
    return scenarios


def cache_fraction_scenarios() -> list[dict]:
    scenarios = [
        {
            "variant": "serial",
            "arguments": ["-test", "T02", "-s04-fractions"],
        },
        {
            "variant": "two-rank-split",
            "ranks": 2,
            "arguments": [
                "-parallel",
                "-test", "T02",
                "-s04-fractions",
                "-src-rank", "0",
                "-tgt-rank", "1",
                "-opt-switch", "localAMIComm=1",
            ],
        },
    ]
    for local_comm in (0, 1, 2):
        scenarios.append(
            {
                "variant": f"four-rank-noncontiguous-localAMIComm-{local_comm}",
                "ranks": 4,
                "arguments": [
                    "-parallel",
                    "-test", "T02",
                    "-s04-fractions",
                    "-src-rank", "1",
                    "-tgt-rank", "3",
                    "-opt-switch", f"localAMIComm={local_comm}",
                ],
            }
        )
    return scenarios


def state_transition_scenarios() -> list[dict]:
    scenarios = [
        {
            "variant": "serial",
            "arguments": [
                "-test", "T09",
                "-state-cycle",
                "-opt-switch", "localAMIComm=1",
            ],
        },
        {
            "variant": "two-rank-localAMIComm-1",
            "ranks": 2,
            "arguments": [
                "-parallel",
                "-test", "T09",
                "-state-cycle",
                "-opt-switch", "localAMIComm=1",
            ],
        },
    ]
    for local_comm in (0, 1, 2):
        scenarios.append(
            {
                "variant": f"four-rank-localAMIComm-{local_comm}",
                "ranks": 4,
                "arguments": [
                    "-parallel",
                    "-test", "T09",
                    "-state-cycle",
                    "-opt-switch", f"localAMIComm={local_comm}",
                ],
            }
        )
    return scenarios


def scale_offset_scenarios() -> list[dict]:
    cases = (
        ("unit-control", "1", "0"),
        ("minimum-scale", "1e-9", "0"),
        ("small-scale-offset", "1e-6", "1e-6"),
        ("medium-scale-offset", "1e-3", "1e6"),
        ("large-scale-large-offset", "1", "1e12"),
    )
    scenarios = []
    for name, scale, offset in cases:
        for ranks in (None, 2):
            rank_name = "serial" if ranks is None else "two-rank"
            arguments = [
                "-test", "T11",
                "-scale", scale,
                "-offset", offset,
            ]
            if ranks:
                arguments.extend(
                    (
                        "-parallel",
                        "-src-rank", "0",
                        "-tgt-rank", "1",
                        "-opt-switch", "localAMIComm=1",
                    )
                )
            scenarios.append(
                {
                    "variant": f"{name}-{rank_name}",
                    "ranks": ranks,
                    "arguments": arguments,
                }
            )
    return scenarios


S03_SCENARIOS = {
    "T08": topology_scenarios(),
    "T09": state_transition_scenarios(),
    "T11": scale_offset_scenarios(),
}

S04_SCENARIOS = {
    "T02": cache_fraction_scenarios(),
    "T08": topology_scenarios(),
    "T09": state_transition_scenarios(),
}

T03_ANGLE_53_13 = 53.13010235415598

T03_SERIAL_VARIANTS = (
    ("rot90", "rotational", (0, 0, 1), (0, 0, 0), 90.0, (0, 0, 0)),
    ("rot0", "rotational", (0, 0, 1), (0, 0, 0), 0.0, (0, 0, 0)),
    ("trans0", "translational", (0, 0, 1), (0, 0, 0), 0.0, (0, 0, 0)),
    ("rot180", "rotational", (0, 0, 1), (0, 0, 0), 180.0, (0, 0, 0)),
    ("rot90off", "rotational", (0, 0, 1), (0.5, 0.5, 0), 90.0, (0, 0, 0)),
    ("rotpos_z", "rotational", (0, 0, 1), (0, 0, 0), T03_ANGLE_53_13, (0, 0, 0)),
    ("rotneg_z", "rotational", (0, 0, 1), (0, 0, 0), -T03_ANGLE_53_13, (0, 0, 0)),
    ("rot_ax34", "rotational", (3, 4, 0), (0.5, 0.5, 0), T03_ANGLE_53_13, (0, 0, 0)),
    ("rotpos_z30", "rotational", (0, 0, 1), (0, 0, 0), 30.0, (0, 0, 0)),
    ("rotneg_z30", "rotational", (0, 0, 1), (0, 0, 0), -30.0, (0, 0, 0)),
    ("trans_y2", "translational", (0, 0, 1), (0, 0, 0), 0.0, (0, 2, 0)),
)

T03_PARALLEL_VARIANTS = (
    ("rot90", "rotational", (0, 0, 1), (0, 0, 0), 90.0, (0, 0, 0)),
    ("rotpos_z", "rotational", (0, 0, 1), (0, 0, 0), T03_ANGLE_53_13, (0, 0, 0)),
)


def t03_s05_scenarios() -> list[dict]:
    scenarios = []
    for name, transform, axis, centre, angle, sep in T03_SERIAL_VARIANTS:
        scenarios.append(
            {
                "variant": f"{name}-serial",
                "variant_config": {
                    "transform": transform,
                    "axis": axis,
                    "centre": centre,
                    "angle": angle,
                    "sep": sep,
                },
                "arguments": ["-test", "T03"],
            }
        )
    for name, transform, axis, centre, angle, sep in T03_PARALLEL_VARIANTS:
        scenarios.append(
            {
                "variant": f"{name}-two-rank-distributed",
                "ranks": 2,
                "decompose": True,
                "variant_config": {
                    "transform": transform,
                    "axis": axis,
                    "centre": centre,
                    "angle": angle,
                    "sep": sep,
                },
                "arguments": [
                    "-parallel",
                    "-test", "T03",
                    "-opt-switch", "localAMIComm=0",
                ],
            }
        )
    return scenarios


S05_SCENARIOS = {
    "T03": t03_s05_scenarios(),
    "T08": topology_scenarios(),
    "T09": state_transition_scenarios(),
}

T04_S06_SCENARIOS = [
    {
        "variant": "low-weight-rotation-serial",
        "arguments": ["-test", "T04"],
    },
    {
        "variant": "low-weight-rotation-two-rank-distributed",
        "ranks": 2,
        "decompose": True,
        "rank_system_files": ("controlDict",),
        "arguments": ["-parallel", "-test", "T04"],
    },
]

S06_SCENARIOS = {
    "T04": T04_S06_SCENARIOS,
    "T08": topology_scenarios(),
    "T09": state_transition_scenarios(),
}

T05_S07_ORDERINGS = ("as-built", "src-reversed", "tgt-reversed", "both-reversed")


def t05_s07_scenarios() -> list[dict]:
    base_assertions = {
        "T05_CALCULATED",
        "T05_SOURCE_COVERAGE",
        "T05_TARGET_COVERAGE",
    }
    all_assertions = base_assertions | {"T05_PHYSICAL_MAPPING"}
    swapped_ownership = ["-src-rank", "1", "-tgt-rank", "0"]
    scenarios = []
    for order in T05_S07_ORDERINGS:
        if order != "as-built":
            scenarios.append(
                {
                    "variant": f"{order}-serial",
                    "assertions": set(all_assertions),
                    "arguments": [
                        "-test", "T05",
                        "-no-restart",
                        "-t05-order", order,
                    ],
                }
            )
        scenarios.append(
            {
                "variant": f"{order}-two-rank-distributed",
                "ranks": 2,
                "assertions": set(all_assertions),
                "arguments": [
                    "-parallel",
                    "-test", "T05",
                    "-no-restart",
                    "-t05-order", order,
                ],
            }
        )
        scenarios.append(
            {
                "variant": f"restart-{order}-serial",
                "assertions": set(all_assertions),
                "arguments": ["-test", "T05", "-t05-order", order],
            }
        )
        scenarios.append(
            {
                "variant": f"restart-{order}-two-rank-distributed",
                "ranks": 2,
                "assertions": set(all_assertions),
                "arguments": [
                    "-parallel",
                    "-test", "T05",
                    "-t05-order", order,
                ],
            }
        )
        scenarios.append(
            {
                "variant": f"restart-{order}-two-rank-swapped-ownership",
                "ranks": 2,
                "assertions": set(all_assertions),
                "arguments": [
                    "-parallel",
                    "-test", "T05",
                    "-t05-order", order,
                    *swapped_ownership,
                ],
            }
        )
    scenarios.append(
        {
            "variant": "multi-target-restart-serial",
            "assertions": set(all_assertions),
            "arguments": ["-test", "T05", "-t05-multi"],
        }
    )
    scenarios.append(
        {
            "variant": "multi-target-restart-two-rank-distributed",
            "ranks": 2,
            "assertions": set(all_assertions),
            "arguments": ["-parallel", "-test", "T05", "-t05-multi"],
        }
    )
    scenarios.append(
        {
            "variant": "multi-target-restart-two-rank-swapped-ownership",
            "ranks": 2,
            "assertions": set(all_assertions),
            "arguments": [
                "-parallel",
                "-test", "T05",
                "-t05-multi",
                *swapped_ownership,
            ],
        }
    )
    return scenarios


S04_REQUIRED_TESTS = {"T01", "T02", "T08", "T09", "T13"}

S03_REQUIRED_TESTS = {
    "T01", "T02", "T03", "T04", "T05", "T06", "T07",
    "T08", "T09", "T10", "T11", "T12", "T13",
}

TIMING_PHASES = {
    "T01": {
        "constructor",
        "mapping_rebuild",
        "interpolate_to_source",
        "interpolate_to_target",
        "checksum",
    },
    "T02": {"cache_restore", "cached_interpolation"},
    "T03": {"topology_rebuild"},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_process(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
    timeout: float,
) -> dict:
    started = time.monotonic()
    process_env = env.copy()
    process_env["PWD"] = str(cwd.resolve())
    process_env["PYTHONDONTWRITEBYTECODE"] = "1"
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=process_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    timed_out = False
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            output, _ = process.communicate(timeout=2)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = process.communicate()

    return {
        "command": command,
        "cwd": str(cwd),
        "return_code": process.returncode,
        "timed_out": timed_out,
        "duration_seconds": time.monotonic() - started,
        "output": output,
    }


def run_profiled_process(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
    timeout: float,
    metrics_path: Path,
) -> dict:
    time_command = Path("/usr/bin/time")
    if not time_command.is_file():
        raise ValueError("/usr/bin/time is required for performance measurements")

    wrapped_command = [
        str(time_command),
        "-v",
        "-o",
        str(metrics_path),
        "--",
        *command,
    ]
    result = run_process(wrapped_command, cwd, env, timeout)
    result["executed_command"] = result["command"]
    result["command"] = command

    metrics = {}
    if metrics_path.is_file():
        text = metrics_path.read_text(encoding="utf-8")
        for key, pattern in {
            "peak_rss_kb": r"Maximum resident set size \(kbytes\):\s*(\d+)",
            "user_seconds": r"User time \(seconds\):\s*([0-9.]+)",
            "system_seconds": r"System time \(seconds\):\s*([0-9.]+)",
            "elapsed": r"Elapsed \(wall clock\) time .*?:\s*([0-9:.]+)",
        }.items():
            match = re.search(pattern, text)
            if match:
                metrics[key] = (
                    int(match.group(1))
                    if key == "peak_rss_kb"
                    else match.group(1)
                    if key == "elapsed"
                    else float(match.group(1))
                )
    result["performance"] = {
        "metrics": metrics,
        "metrics_log": str(metrics_path),
        "peak_rss_scope": (
            "direct process; MPI launcher peak is not the sum of rank memory"
        ),
    }
    return result


def write_log(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(result["output"], encoding="utf-8")


def parse_assertions(output: str) -> list[dict]:
    assertions = []
    for line in output.splitlines():
        marker = line.find("AMI_ASSERT ")
        if marker < 0:
            continue
        values = {}
        for token in shlex.split(line[marker + len("AMI_ASSERT ") :]):
            if "=" in token:
                key, value = token.split("=", 1)
                values[key] = value
        if {"test", "id", "status"} <= values.keys():
            assertions.append(values)
    return assertions


def parse_timing_records(
    output: str,
    test_id: str,
    performance_iterations: int,
) -> list[dict]:
    records = []
    for line in output.splitlines():
        marker = line.find("AMI_TIMING ")
        if marker < 0:
            continue
        fields = {}
        for token in shlex.split(line[marker + len("AMI_TIMING ") :]):
            if "=" in token:
                key, value = token.split("=", 1)
                fields[key] = value
        if fields.get("test") != test_id:
            continue
        if not {"phase", "iterations"} <= fields.keys():
            raise ValueError(f"incomplete timing record for {test_id}")
        phase = fields["phase"]
        try:
            iterations = int(fields["iterations"])
        except ValueError as error:
            raise ValueError(
                f"invalid iteration count in {test_id}/{phase}"
            ) from error
        expected_iterations = (
            performance_iterations if test_id == "T01" else 1
        )
        if iterations != expected_iterations:
            raise ValueError(
                f"unexpected iteration count in {test_id}/{phase}"
            )

        record = {
            "test": test_id,
            "phase": phase,
            "iterations": iterations,
        }
        if {"total_seconds", "mean_seconds"} <= fields.keys():
            try:
                total_seconds = float(fields["total_seconds"])
                mean_seconds = float(fields["mean_seconds"])
            except ValueError as error:
                raise ValueError(
                    f"invalid elapsed time in {test_id}/{phase}"
                ) from error
            if (
                not math.isfinite(total_seconds)
                or not math.isfinite(mean_seconds)
                or total_seconds < 0
                or mean_seconds < 0
            ):
                raise ValueError(
                    f"non-finite or negative time in {test_id}/{phase}"
                )
            record["total_seconds"] = total_seconds
            record["mean_seconds"] = mean_seconds
        elif phase == "checksum" and "value" in fields:
            try:
                value = float(fields["value"])
            except ValueError as error:
                raise ValueError(
                    f"invalid checksum in {test_id}/{phase}"
                ) from error
            if not math.isfinite(value):
                raise ValueError(f"non-finite checksum in {test_id}/{phase}")
            record["value"] = value
        else:
            raise ValueError(f"missing timing values in {test_id}/{phase}")
        records.append(record)

    phases = [record["phase"] for record in records]
    if set(phases) != TIMING_PHASES[test_id] or len(phases) != len(set(phases)):
        raise ValueError(f"missing or duplicate timing phases for {test_id}")
    return records


def terminate_self_check(
    run_dir: Path,
    env: dict[str, str],
) -> dict:
    child = "import time; time.sleep(30)"
    parent = (
        "import subprocess,sys,time; "
        f"child=subprocess.Popen([sys.executable, '-c', {child!r}]); "
        "print(child.pid, flush=True); time.sleep(30)"
    )
    result = run_process(
        [sys.executable, "-c", parent],
        run_dir,
        env,
        0.25,
    )
    log_path = run_dir / "T00-timeout-self-check.log"
    write_log(log_path, result)

    child_pid = None
    for line in result["output"].splitlines():
        try:
            child_pid = int(line.strip())
        except ValueError:
            continue
        break

    child_stopped = child_pid is not None
    if child_pid is not None:
        for _ in range(20):
            try:
                stat = Path(f"/proc/{child_pid}/stat").read_text(encoding="utf-8")
            except (FileNotFoundError, ProcessLookupError):
                break
            if len(stat.split()) > 2 and stat.split()[2] == "Z":
                break
            child_stopped = False
            time.sleep(0.05)

    passed = (
        result["timed_out"]
        and result["return_code"] != 0
        and child_stopped
    )
    return {
        "id": "T00_TIMEOUT",
        "status": "pass" if passed else "failure",
        "timed_out": result["timed_out"],
        "return_code": result["return_code"],
        "child_pid": child_pid,
        "child_stopped": child_stopped,
        "log": str(log_path),
    }


def expected_status(
    test_id: str,
    spec: dict,
    result: dict,
    assertions: list[dict],
    observation_mode: bool = False,
) -> str:
    required_fields = {
        "test",
        "id",
        "status",
        "observed",
        "expected",
        "error",
        "threshold",
    }
    for assertion in assertions:
        if (
            not required_fields <= assertion.keys()
            or assertion["test"] != test_id
        ):
            return "failure"
        if assertion["status"] not in {"pass", "fail"}:
            return "failure"
        try:
            values = [
                float(assertion[key])
                for key in ("observed", "expected", "error", "threshold")
            ]
        except ValueError:
            return "failure"
        if not all(math.isfinite(value) for value in values):
            return "failure"
        if values[-1] < 0:
            return "failure"

    assertion_ids = [assertion["id"] for assertion in assertions]
    if len(assertion_ids) != len(set(assertion_ids)):
        return "failure"
    expected_ids = set(spec.get("assertions", set()))
    observed_ids = set(assertion_ids)

    failed_ids = {
        assertion["id"]
        for assertion in assertions
        if assertion["status"] == "fail"
    }
    unexpected_ids = {
        assertion["id"]
        for assertion in assertions
        if assertion["status"] not in {"pass", "fail"}
    }
    expected_failures = set(spec.get("expected_failures", set()))

    if result["timed_out"] or unexpected_ids:
        return "failure"
    expected_fatal = spec.get("expected_fatal")
    if expected_fatal:
        if (
            result["return_code"] != 0
            and not assertions
            and result["output"].count(expected_fatal)
            >= (spec.get("ranks") or 1)
        ):
            return "observed_failure" if observation_mode else "expected_failure"
        return "failure"
    if observation_mode:
        if (
            result["return_code"] == 0
            and observed_ids == expected_ids
            and not failed_ids
            and assertions
        ):
            return "observed_pass"
        if (
            result["return_code"] == 1
            and failed_ids
            and assertions
            and (
                observed_ids == expected_ids
                or (
                    observed_ids <= expected_ids
                    and "FOAM FATAL ERROR" in result["output"]
                )
            )
        ):
            return "observed_failure"
        return "failure"
    if observed_ids != expected_ids:
        return "failure"
    if expected_failures:
        if (
            result["return_code"] != 0
            and failed_ids == expected_failures
            and assertions
        ):
            return "expected_failure"
        return "failure"
    if result["return_code"] == 0 and not failed_ids and assertions:
        return "pass"
    return "failure"


def prepare_case(
    case_name: str,
    run_dir: Path,
) -> Path:
    source = TEST_DIR / "cases" / case_name
    destination = run_dir / "case"
    shutil.copytree(source, destination)
    return destination


def build_driver(
    build_root: Path,
    run_id: str,
    env: dict[str, str],
    timeout: float,
) -> tuple[Path | None, dict]:
    build_dir = build_root / run_id
    source_dir = build_dir / "source"
    app_bin = build_dir / "bin"
    source_dir.parent.mkdir(parents=True, exist_ok=False)
    shutil.copytree(TEST_DIR, source_dir)
    app_bin.mkdir(parents=True, exist_ok=False)

    build_env = env.copy()
    build_env["FOAM_USER_APPBIN"] = str(app_bin)
    result = run_process(["wmake"], source_dir, build_env, timeout)
    log_path = build_dir / "build.log"
    write_log(log_path, result)
    binary = app_bin / "Test-AMIInterpolation"
    if result["return_code"] != 0 or not binary.is_file():
        return None, {
            "status": "failure",
            "command": result["command"],
            "return_code": result["return_code"],
            "timed_out": result["timed_out"],
            "duration_seconds": result["duration_seconds"],
            "log": str(log_path),
        }
    return binary, {
        "status": "pass",
        "command": result["command"],
        "return_code": result["return_code"],
        "duration_seconds": result["duration_seconds"],
        "log": str(log_path),
        "binary": str(binary),
        "sha256": sha256_file(binary),
    }


def inspect_target(
    target_id: str,
    project_dir: Path,
    binary: Path,
    expected_options: str,
    env: dict[str, str],
    work_dir: Path,
    allowed_library_roots: tuple[Path, ...] = (),
) -> dict:
    expected_project = {
        "B0": B0_PROJECT_DIR,
        "C0": REPO_ROOT,
    }[target_id].resolve()
    if project_dir.resolve() != expected_project:
        raise ValueError(
            f"{target_id} project directory must be {expected_project}"
        )
    actual_project = Path(env.get("WM_PROJECT_DIR", "")).resolve()
    if actual_project != project_dir.resolve():
        raise ValueError(
            f"WM_PROJECT_DIR={actual_project} does not match "
            f"--project-dir={project_dir.resolve()}"
        )
    if env.get("WM_PROJECT_VERSION") != "v2606":
        raise ValueError(f"{target_id} requires a v2606 environment")
    if env.get("WM_OPTIONS") != expected_options:
        raise ValueError(
            f"WM_OPTIONS={env.get('WM_OPTIONS')} does not match "
            f"{expected_options}"
        )

    ldd = run_process(
        ["ldd", str(binary.resolve())],
        work_dir,
        env,
        20,
    )
    if ldd["return_code"] != 0 or ldd["timed_out"]:
        raise ValueError(f"ldd failed for {binary}: {ldd['output']}")

    resolved = {}
    library_roots = (project_dir.resolve(),) + tuple(
        root.resolve() for root in allowed_library_roots
    )
    for library in (
        "libOpenFOAM.so",
        "libfiniteVolume.so",
        "libmeshTools.so",
        "libPstream.so",
    ):
        match = re.search(
            rf"{re.escape(library)}\s+=>\s+(\S+)",
            ldd["output"],
        )
        if not match:
            raise ValueError(f"ldd did not resolve required library {library}")
        library_path = Path(match.group(1)).resolve()
        if not any(
            library_path.is_relative_to(root)
            for root in library_roots
        ):
            raise ValueError(
                f"{library} resolves outside the approved {target_id} "
                f"library roots: {library_path}"
            )
        resolved[library] = str(library_path)

    tool_paths = {}
    for tool in ("wmake", "blockMesh", "checkMesh", "foamDictionary", "mpirun"):
            tool_path = shutil.which(tool, path=env.get("PATH"))
            if not tool_path:
                raise ValueError(f"required tool {tool} is not available")
            resolved_tool = Path(tool_path).resolve()
            if tool != "mpirun" and not resolved_tool.is_relative_to(
                project_dir.resolve()
            ):
                raise ValueError(
                    f"{tool} resolves outside {target_id}: {resolved_tool}"
                )
            tool_paths[tool] = str(resolved_tool)

    wmake_version = run_process(
        [tool_paths["wmake"], "-version"],
        work_dir,
        env,
        20,
    )
    if wmake_version["return_code"] != 0 or wmake_version["timed_out"]:
        raise ValueError(f"wmake -version failed: {wmake_version['output']}")

    commit = None
    worktree_dirty = None
    if target_id == "C0":
        commit = subprocess.run(
            ["git", "-C", str(project_dir), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(project_dir), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        worktree_dirty = bool(status.strip())

    return {
        "target_id": target_id,
        "project_dir": str(project_dir.resolve()),
        "project_version": env["WM_PROJECT_VERSION"],
        "wmake_version": wmake_version["output"].strip(),
        "wm_options": env["WM_OPTIONS"],
        "environment": {
            key: env.get(key)
            for key in (
                "WM_PRECISION_OPTION",
                "WM_LABEL_SIZE",
                "WM_MPLIB",
                "WM_COMPILER",
                "WM_COMPILE_OPTION",
            )
        },
        "commit": commit,
        "worktree_dirty": worktree_dirty,
        "binary": str(binary.resolve()),
        "binary_sha256": sha256_file(binary),
        "ldd": ldd["output"],
        "openfoam_libraries": resolved,
        "tools": tool_paths,
    }


def run_test(
    test_id: str,
    spec: dict,
    binary: Path,
    run_root: Path,
    run_id: str,
    env: dict[str, str],
    timeout: float,
    observation_mode: bool = False,
    measure_performance: bool = False,
    performance_iterations: int = 100,
    scenario_id: str | None = None,
) -> dict:
    record_id = scenario_id or test_id
    run_dir = run_root / run_id / record_id
    run_dir.mkdir(parents=True, exist_ok=False)
    case_dir = prepare_case(spec["case"], run_dir)
    if spec.get("variant_config"):
        write_t03_variant(case_dir, **spec["variant_config"])

    prep_result = None
    if spec.get("prepare_mesh"):
        block_mesh = shutil.which("blockMesh", path=env.get("PATH"))
        if not block_mesh:
            raise ValueError("blockMesh is not available in the selected target")
        prep_command = [block_mesh, "-case", str(case_dir)]
        prep_result = (
            run_profiled_process(
                prep_command,
                run_dir,
                env,
                timeout,
                run_dir / "time.blockMesh.txt",
            )
            if measure_performance
            else run_process(prep_command, run_dir, env, timeout)
        )
        write_log(run_dir / "log.blockMesh", prep_result)
        if prep_result["return_code"] != 0 or prep_result["timed_out"]:
            return {
                "id": record_id,
                "test_id": test_id,
                "name": spec["name"],
                "variant": spec.get("variant", "default"),
                "case": spec["case"],
                "rank_count": spec.get("ranks") or 1,
                "status": "failure",
                "reason": "blockMesh fixture preparation failed",
                "prepare": {
                    "return_code": prep_result["return_code"],
                    "timed_out": prep_result["timed_out"],
                    "log": str(run_dir / "log.blockMesh"),
                },
            }

    if spec.get("decompose"):
        decompose_command = shutil.which(
            "decomposePar", path=env.get("PATH")
        )
        if not decompose_command:
            raise ValueError("decomposePar is not available in the selected target")
        decompose_result = run_process(
            [decompose_command, "-force", "-case", str(case_dir)],
            run_dir,
            env,
            timeout,
        )
        write_log(run_dir / "log.decomposePar", decompose_result)
        if decompose_result["return_code"] != 0 or decompose_result["timed_out"]:
            return {
                "id": record_id,
                "test_id": test_id,
                "name": spec["name"],
                "variant": spec.get("variant", "default"),
                "case": spec["case"],
                "rank_count": spec.get("ranks") or 1,
                "status": "failure",
                "reason": "decomposePar fixture preparation failed",
                "prepare": {
                    "return_code": decompose_result["return_code"],
                    "timed_out": decompose_result["timed_out"],
                    "log": str(run_dir / "log.decomposePar"),
                },
            }
        for rank_dir in sorted(case_dir.glob("processor*")):
            (rank_dir / "system").mkdir(exist_ok=True)
            for name in spec.get("rank_system_files", ("t03Variant", "controlDict")):
                shutil.copy2(
                    case_dir / "system" / name,
                    rank_dir / "system" / name,
                )

    command = [str(binary.resolve())]
    if spec.get("ranks"):
        mpirun = shutil.which("mpirun", path=env.get("PATH"))
        if not mpirun:
            raise ValueError("mpirun is not available in the selected target")
        command = [
            mpirun,
            "-np",
            str(spec["ranks"]),
            *command,
        ]
    command.extend(["-case", str(case_dir), *spec["arguments"]])

    if measure_performance and test_id in {"T01", "T02", "T03"}:
        command.append("-timing")
        if test_id == "T01":
            command.extend(
                ["-benchmark-iterations", str(performance_iterations)]
            )
    result = (
        run_profiled_process(
            command,
            run_dir,
            env,
            timeout,
            run_dir / "time.test.txt",
        )
        if measure_performance
        else run_process(command, run_dir, env, timeout)
    )
    log_path = run_dir / "log.test"
    write_log(log_path, result)
    assertions = parse_assertions(result["output"])
    status = expected_status(
        test_id,
        spec,
        result,
        assertions,
        observation_mode,
    )
    timing_records = None
    timing_error = None
    if measure_performance and test_id in TIMING_PHASES:
        try:
            timing_records = parse_timing_records(
                result["output"],
                test_id,
                performance_iterations,
            )
        except ValueError as error:
            timing_error = str(error)
            status = "failure"

    record = {
        "id": record_id,
        "test_id": test_id,
        "name": spec["name"],
        "variant": spec.get("variant", "default"),
        "case": spec["case"],
        "rank_count": spec.get("ranks") or 1,
        "status": status,
        "command": command,
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "assertions": assertions,
        "measured_values": [
            {
                "assertion_id": assertion["id"],
                "observed": assertion.get("observed"),
                "expected": assertion.get("expected"),
            }
            for assertion in assertions
        ],
        "thresholds": [
            {
                "assertion_id": assertion["id"],
                "absolute_threshold": assertion.get("threshold"),
            }
            for assertion in assertions
        ],
        "failed_assertion_ids": sorted(
            {
                assertion["id"]
                for assertion in assertions
                if assertion["status"] == "fail"
            }
        ),
        "expected_failure_ids": (
            []
            if observation_mode
            else sorted(spec.get("expected_failures", set()))
        ),
        "log": str(log_path),
    }
    if spec.get("expected_fatal"):
        record["expected_fatal"] = spec["expected_fatal"]
    if measure_performance:
        record["performance"] = result["performance"]
        if test_id in TIMING_PHASES:
            record["timings"] = timing_records
            if timing_error:
                record["performance_timing_error"] = timing_error
    if observation_mode:
        record["c0_known_failure_ids"] = sorted(
            spec.get("expected_failures", set())
        )
    if prep_result:
        record["prepare"] = {
            "return_code": prep_result["return_code"],
            "timed_out": prep_result["timed_out"],
            "log": str(run_dir / "log.blockMesh"),
        }
        if measure_performance:
            record["prepare"]["performance"] = prep_result["performance"]
    if result["timed_out"]:
        record["reason"] = "test exceeded its configured timeout"
    elif timing_error:
        record["reason"] = (
            f"performance timing output is invalid: {timing_error}"
        )
    elif status == "failure":
        record["reason"] = "exit status or assertion IDs did not match the test contract"
    return record


def run_cache_topology_rejection(
    run_root: Path,
    run_id: str,
    env: dict[str, str],
    timeout: float,
) -> dict:
    record_id = "T09-cache-plus-face-creation"
    run_dir = run_root / run_id / record_id
    run_dir.mkdir(parents=True, exist_ok=False)
    case_dir = prepare_case("topology", run_dir)
    block_mesh_dict = case_dir / "system/blockMeshDict"
    contents = block_mesh_dict.read_text(encoding="utf-8")

    for patch_name in ("AMI1", "AMI2"):
        pattern = re.compile(
            rf"(?m)^([ \t]*{patch_name}[ \t]*\n[ \t]*\{{\n)"
        )
        entries = "        cacheSize 2;\n        createAMIFaces true;\n"
        contents, count = pattern.subn(
            lambda match: match.group(1) + entries,
            contents,
            count=1,
        )
        if count != 1:
            raise ValueError(
                f"could not configure {patch_name} for cache/topology test"
            )
    block_mesh_dict.write_text(contents, encoding="utf-8")

    block_mesh = shutil.which("blockMesh", path=env.get("PATH"))
    if not block_mesh:
        raise ValueError("blockMesh is not available in the selected target")
    command = [block_mesh, "-case", str(case_dir)]
    result = run_process(command, run_dir, env, timeout)
    log_path = run_dir / "log.blockMesh"
    write_log(log_path, result)
    expected_message = (
        "createAMIFaces and AMI caching cannot be used together."
    )
    rejected = (
        result["return_code"] != 0
        and not result["timed_out"]
        and expected_message in result["output"]
    )

    return {
        "id": record_id,
        "test_id": "T09",
        "name": "unsupported cache plus face-creation combination",
        "variant": "explicit-fatal-io-rejection",
        "case": "topology",
        "rank_count": 1,
        "status": "observed_pass" if rejected else "failure",
        "expected_behavior": expected_message,
        "observed_rejection": rejected,
        "command": command,
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "log": str(log_path),
        "reason": (
            ""
            if rejected
            else "unsupported combination did not produce the expected explicit rejection"
        ),
    }


def run_t05_unmatched_rejection(
    run_root: Path,
    run_id: str,
    env: dict[str, str],
    timeout: float,
    binary: Path,
    ranks: int,
    restart: bool = False,
) -> dict:
    layout = "two-rank-distributed" if ranks > 1 else "serial"
    restart_part = "restart-" if restart else ""
    record_id = f"T05-unmatched-source-{restart_part}{layout}"
    run_dir = run_root / run_id / record_id
    run_dir.mkdir(parents=True, exist_ok=False)
    case_dir = prepare_case("analytical", run_dir)

    command = [str(binary.resolve())]
    if ranks > 1:
        mpirun = shutil.which("mpirun", path=env.get("PATH"))
        if not mpirun:
            raise ValueError("mpirun is not available in the selected target")
        command = [mpirun, "-np", str(ranks), *command]
    command.extend(["-case", str(case_dir)])
    if ranks > 1:
        command.append("-parallel")
    command.extend(["-test", "T05", "-unmatched-source"])
    if not restart:
        command.append("-no-restart")

    result = run_process(command, run_dir, env, timeout)
    log_path = run_dir / "log.test"
    write_log(log_path, result)
    expected_message = "source faces have no overlapping target face"
    rejected = (
        result["return_code"] != 0
        and not result["timed_out"]
        and expected_message in result["output"]
    )

    return {
        "id": record_id,
        "test_id": "T05",
        "name": (
            "unmatched source face with required match and restart"
            if restart
            else "unmatched source face with required match and no restart"
        ),
        "variant": (
            "restart-explicit-fatal-required-match"
            if restart
            else "explicit-fatal-required-match"
        ),
        "case": "analytical",
        "rank_count": ranks,
        "status": "pass" if rejected else "failure",
        "expected_behavior": expected_message,
        "observed_rejection": rejected,
        "command": command,
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "log": str(log_path),
        "reason": (
            ""
            if rejected
            else "unmatched source face did not produce the expected collective fatal"
        ),
    }


def run_s04_integration_screen(
    run_root: Path,
    results_dir: Path,
    run_id: str,
    project_dir: Path,
    library_root: Path,
    solid_foam: Path,
    env: dict[str, str],
    target_id: str,
    timeout: float,
) -> dict:
    integration_results = results_dir / "integration"
    integration_results.mkdir(parents=True, exist_ok=True)
    integration_report = (
        integration_results / f"S04-cache-integration-{run_id}.json"
    )
    if integration_report.exists():
        raise FileExistsError(
            f"refusing to overwrite {integration_report}"
        )

    runner_dir = run_root / run_id / "S04-cache-integration-runner"
    runner_dir.mkdir(parents=True, exist_ok=False)
    command = [
        sys.executable,
        str(TEST_DIR / "runS04CacheIntegration.py"),
        "--target-id",
        target_id,
        "--project-dir",
        str(project_dir),
        "--library-root",
        str(library_root),
        "--solid-foam",
        str(solid_foam),
        "--run-root",
        str(run_root / f"{run_id}-integration"),
        "--results-dir",
        str(integration_results),
        "--report",
        str(integration_report),
        "--timeout",
        str(timeout),
    ]
    result = run_process(
        command,
        TEST_DIR,
        env,
        max(timeout*10, 600),
    )
    log_path = runner_dir / "log.integration-runner"
    write_log(log_path, result)

    payload = None
    reason = None
    if integration_report.is_file():
        try:
            payload = json.loads(
                integration_report.read_text(encoding="utf-8")
            )
        except json.JSONDecodeError as error:
            reason = f"invalid integration report: {error}"
    else:
        reason = "integration report was not created"

    passed = (
        result["return_code"] == 0
        and not result["timed_out"]
        and isinstance(payload, dict)
        and payload.get("status") == "pass"
    )
    if reason is None and not passed:
        reason = "fractional field/GAMG integration screen failed"

    return {
        "id": "T18-fractional-cache-field-gamg",
        "test_id": "T18",
        "status": "pass" if passed else "failure",
        "command": command,
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "timeout_seconds": max(timeout*10, 600),
        "log": str(log_path),
        "report": str(integration_report),
        "result": payload,
        "reason": reason or "",
    }


def fmt_number(value: float | int) -> str:
    return f"{value:.16g}"


def fmt_vector(values: list[float | int]) -> str:
    return "(" + " ".join(fmt_number(value) for value in values) + ")"


def render_solid_case(
    template_dir: Path,
    case_dir: Path,
    values: dict,
) -> None:
    shutil.copytree(template_dir, case_dir)
    replacements = {
        "@LEFT_NY@": str(values["leftNy"]),
        "@RIGHT_NY@": str(values["rightNy"]),
        "@LEFT_RHO@": fmt_number(values["leftRho"]),
        "@RIGHT_RHO@": fmt_number(values["rightRho"]),
        "@LEFT_CP@": fmt_number(values["leftCp"]),
        "@RIGHT_CP@": fmt_number(values["rightCp"]),
        "@LEFT_KAPPA@": fmt_vector(values["leftKappa"]),
        "@RIGHT_KAPPA@": fmt_vector(values["rightKappa"]),
        "@INITIAL_T@": fmt_number(values["initialT"]),
        "@HOT_T@": fmt_number(values["hotT"]),
        "@COLD_T@": fmt_number(values["coldT"]),
        "@END_TIME@": fmt_number(values["endTime"]),
        "@DELTA_T@": fmt_number(values["deltaT"]),
        "@E1@": fmt_vector(values["e1"]),
        "@E2@": fmt_vector(values["e2"]),
    }

    for path in case_dir.rglob("*"):
        if not path.is_file():
            continue
        contents = path.read_text(encoding="utf-8")
        for token, replacement in replacements.items():
            contents = contents.replace(token, replacement)
        if re.search(r"@[A-Z0-9_]+@", contents):
            raise ValueError(f"unrendered token in {path}")
        path.write_text(contents, encoding="utf-8")


def parse_foam_scalar(output: str, entry: str) -> float:
    try:
        value = float(output.strip())
    except ValueError as error:
        raise ValueError(f"foamDictionary returned invalid scalar for {entry}") from error
    if not math.isfinite(value):
        raise ValueError(f"foamDictionary returned non-finite scalar for {entry}")
    return value


def parse_foam_vector(output: str, entry: str) -> list[float]:
    values = [
        float(value)
        for value in re.findall
        (
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
            output,
        )
    ]
    if len(values) != 3 or not all(math.isfinite(value) for value in values):
        raise ValueError(f"foamDictionary returned invalid vector for {entry}")
    return values


def validate_material_properties(
    foam_dictionary: str,
    case_dir: Path,
    values: dict,
    env: dict[str, str],
    timeout: float,
    log_path: Path,
) -> tuple[bool, dict, list[str]]:
    observed = {}
    errors = []
    log_parts = []

    for zone in ("left", "right"):
        expected = {
            "rho": values[f"{zone}Rho"],
            "Cp": values[f"{zone}Cp"],
            "kappa": values[f"{zone}Kappa"],
        }
        actual = {}
        entries = {
            "rho": f"mixture.{zone}.equationOfState.rho",
            "Cp": f"mixture.{zone}.thermodynamics.Cp",
            "kappa": f"mixture.{zone}.transport.kappa",
        }
        for key, entry in entries.items():
            result = run_process(
                [
                    foam_dictionary,
                    "-entry",
                    entry,
                    "-value",
                    str(case_dir / "constant" / "thermophysicalProperties"),
                ],
                case_dir,
                env,
                timeout,
            )
            log_parts.append(
                f"$ {shlex.join(result['command'])}\n"
                f"return_code={result['return_code']} "
                f"timed_out={result['timed_out']}\n"
                f"{result['output']}"
            )
            if result["return_code"] != 0 or result["timed_out"]:
                errors.append(f"could not read {entry}")
                continue
            try:
                actual[key] = (
                    parse_foam_vector(result["output"], entry)
                    if key == "kappa"
                    else parse_foam_scalar(result["output"], entry)
                )
            except ValueError as error:
                errors.append(str(error))

        if actual.keys() != expected.keys():
            errors.append(f"incomplete property values for zone {zone}")
            matches = False
        else:
            matches = all(
                math.isclose(actual[key], expected[key], rel_tol=1e-12, abs_tol=1e-12)
                if key != "kappa"
                else all(
                    math.isclose(
                        actual_component,
                        expected_component,
                        rel_tol=1e-12,
                        abs_tol=1e-12,
                    )
                    for actual_component, expected_component
                    in zip(actual[key], expected[key])
                )
                for key in expected
            )
            if not matches:
                errors.append(f"property values differ from the manifest for {zone}")

        observed[zone] = {
            "expected": expected,
            "observed": actual,
            "matches": matches,
        }

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(log_parts), encoding="utf-8")
    return not errors, observed, errors


def parse_cell_zones(text: str) -> dict[str, list[int]]:
    pattern = re.compile(
        r"(?ms)^\s*(left|right)\s*\{\s*type\s+cellZone;\s*"
        r"cellLabels\s+List<label>\s+(\d+)\(([^)]*)\);"
    )
    zones = {}
    for match in pattern.finditer(text):
        zone_name = match.group(1)
        declared_count = int(match.group(2))
        labels = [
            int(value)
            for value in re.findall(r"\d+", match.group(3))
        ]
        if declared_count != len(labels):
            continue
        zones[zone_name] = labels
    return zones


def validate_solid_fixtures(
    run_root: Path,
    run_id: str,
    env: dict[str, str],
    timeout: float,
) -> list[dict]:
    template = TEST_DIR / "cases" / "solidFoam" / "template"
    variants = json.loads(
        (
            TEST_DIR / "cases" / "solidFoam" / "variants.json"
        ).read_text(encoding="utf-8")
    )
    block_mesh = shutil.which("blockMesh", path=env.get("PATH"))
    check_mesh = shutil.which("checkMesh", path=env.get("PATH"))
    foam_dictionary = shutil.which("foamDictionary", path=env.get("PATH"))
    if not block_mesh or not check_mesh or not foam_dictionary:
        raise ValueError(
            "blockMesh, checkMesh, and foamDictionary are required "
            "for solid fixture validation"
        )

    records = []
    for test_id, values in variants.items():
        case_dir = run_root / run_id / "solid-fixtures" / test_id
        case_dir.parent.mkdir(parents=True, exist_ok=True)
        render_solid_case(template, case_dir, values)
        logs = {}

        mixture = run_process(
            [
                foam_dictionary,
                "-entry",
                "thermoType.mixture",
                "-value",
                str(case_dir / "constant" / "thermophysicalProperties"),
            ],
            case_dir,
            env,
            timeout,
        )
        write_log(case_dir / "log.foamDictionary", mixture)
        logs["foamDictionary"] = str(case_dir / "log.foamDictionary")
        parsed_mixture = mixture["output"].strip()
        mixture_ok = (
            mixture["return_code"] == 0
            and not mixture["timed_out"]
            and parsed_mixture == "pureZoneMixture"
        )

        properties_ok, properties, property_errors = validate_material_properties(
            foam_dictionary,
            case_dir,
            values,
            env,
            timeout,
            case_dir / "log.materialProperties",
        )
        logs["material_properties"] = str(case_dir / "log.materialProperties")

        mesh_result = run_process(
            [block_mesh, "-case", str(case_dir)],
            case_dir,
            env,
            timeout,
        )
        write_log(case_dir / "log.blockMesh", mesh_result)
        logs["blockMesh"] = str(case_dir / "log.blockMesh")

        check_result = None
        zones_ok = False
        zone_counts = {}
        total_cells = 0
        if mesh_result["return_code"] == 0 and not mesh_result["timed_out"]:
            check_result = run_process(
                [check_mesh, "-case", str(case_dir)],
                case_dir,
                env,
                timeout,
            )
            write_log(case_dir / "log.checkMesh", check_result)
            logs["checkMesh"] = str(case_dir / "log.checkMesh")

            zones_file = case_dir / "constant" / "polyMesh" / "cellZones"
            zones_text = (
                zones_file.read_text(encoding="utf-8")
                if zones_file.is_file()
                else ""
            )
            zones = parse_cell_zones(zones_text)
            zone_counts = {zone: len(labels) for zone, labels in zones.items()}
            expected_counts = {
                "left": 2*values["leftNy"],
                "right": 2*values["rightNy"],
            }
            total_cells = sum(expected_counts.values())
            all_labels = [
                label
                for labels in zones.values()
                for label in labels
            ]
            zones_ok = (
                set(zones) == {"left", "right"}
                and zone_counts == expected_counts
                and len(all_labels) == total_cells
                and len(set(all_labels)) == total_cells
                and set(all_labels) == set(range(total_cells))
            )

        try:
            thermal_reference = solid_reference_metrics(test_id, values)
            reference_ok = True
            reference_error = None
        except (KeyError, TypeError, ValueError) as error:
            thermal_reference = {}
            reference_ok = False
            reference_error = str(error)

        static_ok = (
            mixture_ok
            and properties_ok
            and mesh_result["return_code"] == 0
            and not mesh_result["timed_out"]
            and check_result is not None
            and check_result["return_code"] == 0
            and not check_result["timed_out"]
            and zones_ok
            and reference_ok
        )

        record = {
            "id": test_id,
            "variant": values.get("variant", "default"),
            "rank_count": 1,
            "intended_rank_count": values.get("decompose", 1),
            "status": "fixture_validated" if static_ok else "failure",
            "solver_status": "not_run",
            "solver_reason": "solidFoam runtime characterization is S03",
            "purpose": values["purpose"],
            "mixture": parsed_mixture,
            "mixture_valid": mixture_ok,
            "cell_zones_valid": zones_ok,
            "cell_zone_counts": zone_counts,
            "cell_count": total_cells,
            "material_properties_valid": properties_ok,
            "material_properties": properties,
            "material_property_errors": property_errors,
            "thermal_reference": thermal_reference,
            "reference_valid": reference_ok,
            "reference_error": reference_error,
            "acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
            "logs": logs,
        }
        records.append(record)

    return records


def coverage_status(
    run_records: list[dict],
    fixture_records: list[dict],
    timeout_record: dict,
    reference_test_record: dict,
) -> dict:
    runs_by_test: dict[str, list[dict]] = {}
    for record in run_records:
        test_id = record.get("test_id", record["id"])
        runs_by_test.setdefault(test_id, []).append(record)

    runner_status = (
        "pass"
        if timeout_record["status"] == "pass"
        and reference_test_record["status"] == "pass"
        else "failure"
    )
    fixtures_by_id = {
        record["id"]: record for record in fixture_records
    }

    coverage = {}
    for test_id in ALL_TEST_IDS:
        if test_id == "T00":
            coverage[test_id] = {
                "status": runner_status,
                "reason": (
                    ""
                    if runner_status == "pass"
                    else "timeout or solid reference self-check failed"
                ),
            }
        elif test_id in runs_by_test:
            records = runs_by_test[test_id]
            statuses = [record["status"] for record in records]
            if len(records) == 1:
                status = statuses[0]
            else:
                status = (
                    "complete"
                    if all(value != "failure" for value in statuses)
                    else "failure"
                )
            details = ", ".join(
                f"{record['id']}={record['status']}"
                for record in records
            )
            coverage[test_id] = {
                "status": status,
                "reason": details,
            }
        elif test_id in fixtures_by_id:
            record = fixtures_by_id[test_id]
            coverage[test_id] = {
                "status": record["status"],
                "reason": record.get("reason", ""),
            }
        elif test_id in DEFERRED_TESTS:
            coverage[test_id] = {
                "status": "not_run",
                "reason": DEFERRED_TESTS[test_id],
            }
        else:
            coverage[test_id] = {
                "status": "not_run",
                "reason": "not included in the targeted S02 C0 gate",
            }
    return coverage


def run_reference_tests(
    run_dir: Path,
    env: dict[str, str],
    timeout: float,
) -> dict:
    result = run_process(
        [sys.executable, str(TEST_DIR / "test_solidMetrics.py")],
        run_dir,
        env,
        timeout,
    )
    log_path = run_dir / "T00-solid-reference-tests.log"
    write_log(log_path, result)
    passed = (
        result["return_code"] == 0
        and not result["timed_out"]
        and "OK" in result["output"]
    )
    return {
        "id": "T00_SOLID_REFERENCE",
        "status": "pass" if passed else "failure",
        "command": result["command"],
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "log": str(log_path),
    }


def input_manifest_files() -> list[dict]:
    files = []
    paths = set(TEST_DIR.rglob("*"))
    paths.add(TEST_DIR.parent / "README.md")
    for path in sorted(paths):
        if not path.is_file():
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        files.append(
            {
                "path": str(path.relative_to(REPO_ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return files


def write_input_manifest(
    results_dir: Path,
    run_id: str,
    report_path: Path,
    args: argparse.Namespace,
    identity: dict,
    build_record: dict | None,
) -> Path:
    manifest_path = results_dir / f"{args.stage}-input-manifest-{run_id}.json"
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest = {
        "schema_version": 1,
        "stage": args.stage,
        "run_id": run_id,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "target": identity,
        "command_line": [sys.executable, *sys.argv],
        "selected_tests": args.tests,
        "solid_fixture_validation": args.validate_solid_fixtures,
        "performance_measurement": args.measure_performance,
        "performance_iterations": args.performance_iterations,
        "timeout_seconds": args.timeout,
        "build": build_record,
        "suite_directory": str(TEST_DIR),
        "inputs": input_manifest_files(),
        "report": str(report_path),
    }
    if args.stage == "S04":
        implementation_path = (
            REPO_ROOT
            / "src/meshTools/AMIInterpolation/AMIInterpolation/AMICache.H"
        )
        manifest["library_root"] = str(args.library_root.resolve())
        manifest["solid_foam"] = str(args.solid_foam.resolve())
        manifest["implementation_sources"] = [
            {
                "path": str(implementation_path),
                "sha256": sha256_file(implementation_path),
            }
        ]
    elif args.stage == "S05":
        implementation_path = (
            REPO_ROOT
            / "src/meshTools/AMIInterpolation/patches/cyclicAMI"
            / "cyclicAMIPolyPatch"
            / "cyclicAMIPolyPatchTopologyChange.C"
        )
        library_path = args.library_root.resolve() / "libmeshTools.so"
        if not implementation_path.is_file():
            raise FileNotFoundError(
                f"S05 implementation source is missing: {implementation_path}"
            )
        if not library_path.is_file():
            raise FileNotFoundError(
                f"S05 implementation library is missing: {library_path}"
            )
        manifest["library_root"] = str(args.library_root.resolve())
        manifest["implementation_sources"] = [
            {
                "path": str(implementation_path),
                "sha256": sha256_file(implementation_path),
            }
        ]
        manifest["implementation_libraries"] = [
            {
                "soname": "libmeshTools.so",
                "path": str(library_path),
                "resolved_path": str(library_path.resolve()),
                "sha256": sha256_file(library_path),
            }
        ]
    manifest_path.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest_path


def load_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("S02", "S03", "S04", "S05", "S06"), default="S02"
    )
    parser.add_argument("--target-id", choices=("B0", "C0"), default="C0")
    parser.add_argument(
        "--observation-mode",
        action="store_true",
        help="record baseline outcomes without applying C0 expected-failure rules",
    )
    parser.add_argument("--project-dir", type=Path)
    parser.add_argument("--build-root", type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--library-root", type=Path)
    parser.add_argument("--solid-foam", type=Path)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--validate-solid-fixtures", action="store_true")
    parser.add_argument("--measure-performance", action="store_true")
    parser.add_argument("--performance-iterations", type=int, default=100)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--expected-options", default="linux64GccDPInt32Opt")
    parser.add_argument(
        "--tests",
        nargs="+",
        default=None,
    )
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.tests is None:
        args.tests = [
            "T01", "T02", "T03", "T04", "T05", "T06", "T07",
            "T10", "T12", "T13",
        ]
        if args.stage == "S03":
            args.tests = [
                "T01", "T02", "T03", "T04", "T05", "T06", "T07",
                "T08", "T09", "T10", "T11", "T12", "T13",
            ]
        elif args.stage == "S04":
            args.tests = ["T01", "T02", "T08", "T09", "T13"]
        elif args.stage == "S05":
            args.tests = ["T03", "T08", "T09"]
        elif args.stage == "S06":
            args.tests = ["T04", "T08", "T09"]
    return args


def main() -> int:
    args = load_args()
    if args.list:
        print(json.dumps(
            {
                "tests": {
                    test_id: {
                        "name": spec["name"],
                        "variant": spec.get("variant", "default"),
                        "assertions": sorted(spec.get("assertions", set())),
                        "expected_failures": sorted(
                            spec.get("expected_failures", set())
                        ),
                    }
                    for test_id, spec in TEST_SPECS.items()
                },
                "deferred": DEFERRED_TESTS,
                "scenarios": {
                    test_id: [
                        scenario["variant"]
                        for scenario in scenarios
                    ]
                    for test_id, scenarios in S03_SCENARIOS.items()
                },
                "s04_scenarios": {
                    test_id: [
                        scenario["variant"]
                        for scenario in scenarios
                    ]
                    for test_id, scenarios in S04_SCENARIOS.items()
                },
                "s04_t02_assertions": sorted(T02_S04_ASSERTIONS),
                "solid_fixtures": sorted(
                    json.loads(
                        (
                            TEST_DIR
                            / "cases"
                            / "solidFoam"
                            / "variants.json"
                        ).read_text(encoding="utf-8")
                    )
                ),
            },
            indent=2,
        ))
        return 0

    env = os.environ.copy()
    if not all(
        (
            args.project_dir,
            args.build_root,
            args.run_root,
            args.results_dir,
        )
    ):
        raise ValueError(
            "--project-dir, --build-root, --run-root, and --results-dir are required"
        )

    if args.stage == "S02" and (
        args.target_id != "C0" or args.observation_mode
    ):
        raise ValueError("S02 supports only C0 expected-failure mode")
    if args.stage == "S03" and not args.observation_mode:
        raise ValueError("S03 requires --observation-mode")
    if args.stage == "S04" and (
        args.target_id != "C0" or args.observation_mode
    ):
        raise ValueError("S04 supports only C0 post-fix mode")
    if args.stage == "S05" and (
        args.target_id != "C0" or args.observation_mode
    ):
        raise ValueError("S05 supports only C0 post-fix mode")
    if args.stage == "S06" and (
        args.target_id != "C0" or args.observation_mode
    ):
        raise ValueError("S06 supports only C0 post-fix mode")
    if args.stage == "S04":
        if not args.library_root or not args.solid_foam:
            raise ValueError(
                "S04 requires --library-root and --solid-foam"
            )
        build_root = args.build_root.resolve()
        library_root = args.library_root.resolve(strict=True)
        solid_foam = args.solid_foam.resolve(strict=True)
        protected_roots = {
            REPO_ROOT.resolve(),
            B0_PROJECT_DIR.resolve(strict=True),
        }
        if any(
            path.is_relative_to(root)
            for path in (library_root, solid_foam)
            for root in protected_roots
        ):
            raise ValueError(
                "S04 libraries and solidFoam must be outside OpenFOAM trees"
            )
        if (
            not library_root.is_relative_to(build_root)
            or not solid_foam.is_relative_to(build_root)
        ):
            raise ValueError(
                "S04 libraries and solidFoam must be under --build-root"
            )
        if not library_root.is_dir() or not solid_foam.is_file():
            raise ValueError("S04 library root or solidFoam is missing")
        configured_libbin = Path(env.get("FOAM_LIBBIN", "")).resolve()
        if configured_libbin != library_root:
            raise ValueError(
                f"FOAM_LIBBIN={configured_libbin} does not match "
                f"--library-root={library_root}"
            )
    if args.stage == "S05":
        if not args.library_root:
            raise ValueError("S05 requires --library-root")
        build_root = args.build_root.resolve()
        library_root = args.library_root.resolve(strict=True)
        protected_roots = {REPO_ROOT.resolve(), B0_PROJECT_DIR.resolve(strict=True)}
        if any(library_root.is_relative_to(root) for root in protected_roots):
            raise ValueError(
                "S05 libraries must be outside OpenFOAM trees"
            )
        if not library_root.is_relative_to(build_root):
            raise ValueError("S05 libraries must be under --build-root")
        if not library_root.is_dir():
            raise ValueError("S05 library root is missing")
        configured_libbin = Path(env.get("FOAM_LIBBIN", "")).resolve()
        if configured_libbin != library_root:
            raise ValueError(
                f"FOAM_LIBBIN={configured_libbin} does not match "
                f"--library-root={library_root}"
            )
    if args.measure_performance and args.stage != "S03":
        raise ValueError("--measure-performance is available only in S03")
    if args.measure_performance and args.performance_iterations < 1:
        raise ValueError("--performance-iterations must be positive")

    project_dir = args.project_dir.resolve(strict=True)
    expected_project = {
        "B0": B0_PROJECT_DIR,
        "C0": REPO_ROOT,
    }[args.target_id].resolve(strict=True)
    if project_dir != expected_project:
        raise ValueError(
            f"{args.target_id} project directory must be "
            f"{expected_project}, not {project_dir}"
        )
    actual_project = Path(env.get("WM_PROJECT_DIR", "")).resolve()
    if actual_project != project_dir:
        raise ValueError(
            f"WM_PROJECT_DIR={actual_project} does not match "
            f"--project-dir={project_dir}"
        )
    if env.get("WM_PROJECT_VERSION") != "v2606":
        raise ValueError(f"{args.target_id} requires a v2606 environment")
    if env.get("WM_OPTIONS") != args.expected_options:
        raise ValueError(
            f"WM_OPTIONS={env.get('WM_OPTIONS')} does not match "
            f"{args.expected_options}"
        )
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError("--timeout must be a finite positive number")
    if len(args.tests) != len(set(args.tests)):
        raise ValueError("--tests must not contain duplicate IDs")
    for test_id in args.tests:
        if test_id not in TEST_SPECS:
            raise ValueError(f"unknown runnable test ID: {test_id}")
        if args.stage == "S02" and test_id in S03_SCENARIOS:
            raise ValueError(f"{test_id} is available only in S03")

    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id += "-" + uuid.uuid4().hex[:8]
    build_root = args.build_root.resolve()
    run_root = args.run_root.resolve()
    results_dir = args.results_dir.resolve()
    protected_roots = {REPO_ROOT.resolve()}
    if B0_PROJECT_DIR.exists():
        protected_roots.add(B0_PROJECT_DIR.resolve())
    for artifact_root in (build_root, run_root, results_dir):
        if any(
            artifact_root.is_relative_to(root)
            for root in protected_roots
        ):
            raise ValueError(
                f"test artifacts must stay outside both OpenFOAM trees: "
                f"{artifact_root}"
            )
    build_root.mkdir(parents=True, exist_ok=True)
    run_root.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    run_dir = run_root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    timeout_record = terminate_self_check(run_dir, env)
    reference_test_record = run_reference_tests(run_dir, env, args.timeout)
    binary = args.binary
    build_record = None
    if args.build:
        binary, build_record = build_driver(
            build_root,
            run_id,
            env,
            args.timeout,
        )
        if binary is None:
            report = {
                "schema_version": 3,
                "stage": args.stage,
                "run_id": run_id,
                "target": args.target_id,
                "build": build_record,
                "timeout_self_check": timeout_record,
                "solid_reference_tests": reference_test_record,
                "suite_status": "failure",
            }
            report_path = (
                results_dir
                / f"{args.stage}-{args.target_id}-{run_id}.json"
            )
            report["report_path"] = str(report_path)
            report_path.write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            print(f"Build failed; see {build_record['log']}", file=sys.stderr)
            return 1
    if binary is None:
        raise ValueError("provide --binary or use --build")
    binary = binary.resolve(strict=True)

    identity = inspect_target(
        args.target_id,
        project_dir,
        binary,
        args.expected_options,
        env,
        run_dir,
        (
            (args.library_root.resolve(),)
            if args.stage in {"S04", "S05"}
            else ()
        ),
    )

    run_records = []
    for test_id in args.tests:
        scenario_matrix = (
            S03_SCENARIOS if args.stage == "S03"
            else S04_SCENARIOS if args.stage == "S04"
            else S05_SCENARIOS if args.stage == "S05"
            else S06_SCENARIOS if args.stage == "S06"
            else {}
        )
        if test_id in scenario_matrix:
            for scenario in scenario_matrix[test_id]:
                spec = {**TEST_SPECS[test_id], **scenario}
                if args.stage in {"S05", "S06"}:
                    spec["expected_failures"] = set()
                if args.stage == "S04" and test_id == "T02":
                    spec["assertions"] = T02_S04_ASSERTIONS
                    spec["expected_failures"] = set()
                run_records.append(
                    run_test(
                        test_id,
                        spec,
                        binary,
                        run_root,
                        run_id,
                        env,
                        args.timeout,
                        args.observation_mode,
                        args.measure_performance,
                        args.performance_iterations,
                        scenario_id=f"{test_id}-{scenario['variant']}",
                    )
                )
            if test_id == "T09":
                rejection_record = run_cache_topology_rejection(
                    run_root,
                    run_id,
                    env,
                    args.timeout,
                )
                if args.stage in {"S04", "S05", "S06"}:
                    rejection_record["status"] = (
                        "pass"
                        if rejection_record["status"] == "observed_pass"
                        else "failure"
                    )
                run_records.append(rejection_record)
        else:
            run_records.append(
                run_test(
                    test_id,
                    TEST_SPECS[test_id],
                    binary,
                    run_root,
                    run_id,
                    env,
                    args.timeout,
                    args.observation_mode,
                    args.measure_performance,
                    args.performance_iterations,
                )
            )
            if test_id == "T05" and args.stage == "S02":
                for scenario in t05_s07_scenarios():
                    spec = {**TEST_SPECS[test_id], **scenario}
                    run_records.append(
                        run_test(
                            test_id,
                            spec,
                            binary,
                            run_root,
                            run_id,
                            env,
                            args.timeout,
                            args.observation_mode,
                            args.measure_performance,
                            args.performance_iterations,
                            scenario_id=f"T05-{scenario['variant']}",
                        )
                    )
                for ranks in (1, 2):
                    for restart in (False, True):
                        run_records.append(
                            run_t05_unmatched_rejection(
                                run_root,
                                run_id,
                                env,
                                args.timeout,
                                binary,
                                ranks,
                                restart,
                            )
                        )

    s04_integration_record = None
    if args.stage == "S04":
        s04_integration_record = run_s04_integration_screen(
            run_root,
            results_dir,
            run_id,
            project_dir,
            args.library_root,
            args.solid_foam,
            env,
            args.target_id,
            args.timeout,
        )

    fixture_records = []
    if args.validate_solid_fixtures:
        fixture_records = validate_solid_fixtures(
            run_root,
            run_id,
            env,
            args.timeout,
        )

    coverage = coverage_status(
        run_records,
        fixture_records,
        timeout_record,
        reference_test_record,
    )
    if args.stage == "S02":
        gate_ids = {"T01", "T02", "T03", "T04", "T05", "T06", "T07"}
        gate_records = {
            record["id"]: record
            for record in run_records
            if record["id"] in gate_ids
        }
        gate_pass = gate_ids <= gate_records.keys()
        gate_pass = gate_pass and all(
            gate_records[test_id]["status"] in {"pass", "expected_failure"}
            for test_id in gate_ids
        )
        gate_pass = gate_pass and any(
            record["id"] == "T01" and record["status"] == "pass"
            for record in run_records
        )
        gate_pass = gate_pass and timeout_record["status"] == "pass"
        gate_pass = gate_pass and reference_test_record["status"] == "pass"
        required_solid_ids = {"T14", "T15", "T16", "T17", "T18", "T19"}
        fixture_by_id = {record["id"]: record for record in fixture_records}
        gate_pass = gate_pass and args.validate_solid_fixtures
        gate_pass = gate_pass and required_solid_ids <= fixture_by_id.keys()
        gate_pass = gate_pass and all(
            fixture_by_id[test_id]["status"] == "fixture_validated"
            for test_id in required_solid_ids
        )
        suite_pass = (
            gate_pass
            and reference_test_record["status"] == "pass"
            and all(
                record["status"] in {"pass", "expected_failure"}
                for record in run_records
            )
        )
    elif args.stage == "S03":
        records_by_test: dict[str, list[dict]] = {}
        for record in run_records:
            test_id = record.get("test_id", record["id"])
            records_by_test.setdefault(test_id, []).append(record)

        required_record_ids = {
            test_id
            for test_id in S03_REQUIRED_TESTS
            if test_id not in S03_SCENARIOS
        }
        required_record_ids |= {
            f"{test_id}-{scenario['variant']}"
            for test_id, scenarios in S03_SCENARIOS.items()
            for scenario in scenarios
        }
        required_record_ids.add("T09-cache-plus-face-creation")
        executed_record_ids = {record["id"] for record in run_records}
        valid_observation_statuses = {"observed_pass", "observed_failure"}
        gate_pass = (
            S03_REQUIRED_TESTS <= records_by_test.keys()
            and required_record_ids <= executed_record_ids
            and records_by_test["T01"][0]["status"] == "observed_pass"
            and all(
                record["status"] in valid_observation_statuses
                for record in run_records
            )
            and timeout_record["status"] == "pass"
            and reference_test_record["status"] == "pass"
        )
        suite_pass = gate_pass
    elif args.stage == "S05":
        required_record_ids = {
            f"{test_id}-{scenario['variant']}"
            for test_id, scenarios in S05_SCENARIOS.items()
            for scenario in scenarios
        }
        required_record_ids.add("T09-cache-plus-face-creation")
        records_by_id = {record["id"]: record for record in run_records}
        gate_pass = (
            required_record_ids <= records_by_id.keys()
            and all(
                records_by_id[record_id]["status"] == "pass"
                for record_id in required_record_ids
            )
            and timeout_record["status"] == "pass"
            and reference_test_record["status"] == "pass"
        )
        suite_pass = (
            gate_pass
            and all(record["status"] == "pass" for record in run_records)
        )
    elif args.stage == "S06":
        required_record_ids = {
            f"{test_id}-{scenario['variant']}"
            for test_id, scenarios in S06_SCENARIOS.items()
            for scenario in scenarios
        }
        required_record_ids.add("T09-cache-plus-face-creation")
        records_by_id = {record["id"]: record for record in run_records}
        gate_pass = (
            required_record_ids <= records_by_id.keys()
            and all(
                records_by_id[record_id]["status"] == "pass"
                for record_id in required_record_ids
            )
            and timeout_record["status"] == "pass"
            and reference_test_record["status"] == "pass"
        )
        suite_pass = (
            gate_pass
            and all(record["status"] == "pass" for record in run_records)
        )
    else:
        required_record_ids = {"T01", "T13", "T09-cache-plus-face-creation"}
        required_record_ids |= {
            f"{test_id}-{scenario['variant']}"
            for test_id, scenarios in S04_SCENARIOS.items()
            for scenario in scenarios
        }
        records_by_id = {record["id"]: record for record in run_records}
        gate_pass = (
            required_record_ids <= records_by_id.keys()
            and all(
                records_by_id[record_id]["status"] == "pass"
                for record_id in required_record_ids
            )
            and timeout_record["status"] == "pass"
            and reference_test_record["status"] == "pass"
            and s04_integration_record is not None
            and s04_integration_record["status"] == "pass"
        )
        suite_pass = (
            gate_pass
            and all(record["status"] == "pass" for record in run_records)
        )

    report_path = (
        results_dir / f"{args.stage}-{args.target_id}-{run_id}.json"
    )
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite {report_path}")
    manifest_path = write_input_manifest(
        results_dir,
        run_id,
        report_path,
        args,
        identity,
        build_record,
    )
    report = {
        "schema_version": (
            2 if args.stage == "S02"
            else 3 if args.stage == "S03"
            else 5 if args.stage == "S05"
            else 6 if args.stage == "S06"
            else 4
        ),
        "stage": args.stage,
        "run_id": run_id,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "report_path": str(report_path),
        "input_manifest": str(manifest_path),
        "target": identity,
        "build": build_record,
        "timeout_self_check": timeout_record,
        "solid_reference_tests": reference_test_record,
        "tests": run_records,
        "s04_integration": s04_integration_record,
        "solid_fixtures": fixture_records,
        "coverage": coverage,
        "performance_measurement": args.measure_performance,
        "performance_iterations": args.performance_iterations,
        "s02_exit_gate": (
            "pass" if gate_pass else "failure"
        ) if args.stage == "S02" else "not_applicable",
        "ami_observation_status": (
            "complete" if gate_pass else "failure"
        ) if args.stage == "S03" else "not_applicable",
        "s04_exit_gate": (
            "pass" if gate_pass else "failure"
        ) if args.stage == "S04" else "not_applicable",
        "s05_exit_gate": (
            "pass" if gate_pass else "failure"
        ) if args.stage == "S05" else "not_applicable",
        "s06_exit_gate": (
            "pass" if gate_pass else "failure"
        ) if args.stage == "S06" else "not_applicable",
        "suite_status": (
            "pass" if suite_pass else "failure"
        ) if args.stage in {"S02", "S04", "S05", "S06"} else (
            "ami_observation_complete" if suite_pass else "failure"
        ),
    }
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"Report: {report_path}")
    if args.stage == "S02":
        print(f"S02 exit gate: {report['s02_exit_gate']}")
        return 0 if report["suite_status"] == "pass" else 1
    if args.stage == "S03":
        print(
            "S03 AMI observation: "
            f"{report['ami_observation_status']}"
        )
        return (
            0
            if report["suite_status"] == "ami_observation_complete"
            else 1
        )
    if args.stage == "S05":
        print(f"S05 exit gate: {report['s05_exit_gate']}")
        return 0 if report["suite_status"] == "pass" else 1
    if args.stage == "S06":
        print(f"S06 exit gate: {report['s06_exit_gate']}")
        return 0 if report["suite_status"] == "pass" else 1
    print(f"S04 exit gate: {report['s04_exit_gate']}")
    return 0 if report["suite_status"] == "pass" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"AMI test runner error: {error}", file=sys.stderr)
        raise SystemExit(2)
