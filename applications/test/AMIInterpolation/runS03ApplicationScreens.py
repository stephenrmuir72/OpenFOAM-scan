#!/usr/bin/env python3
"""Capture isolated B0/C0 solidFoam and mixer-tutorial S03 screens."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid

sys.dont_write_bytecode = True
import runTests
from solidMetrics import (
    SOLID_ACCEPTANCE_LIMITS,
    series_resistance_heat_rate,
    series_resistance_interface_temperature,
    solid_reference_metrics,
)


TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]
B0_PROJECT_DIR = Path("/usr/lib/openfoam/openfoam2606")
MIXER_ROOT = Path(
    "tutorials/incompressible/pimpleFoam/laminar/mixerVesselAMI2D"
)
MIXER_CASES = (
    "mixerVesselAMI2D",
    "mixerVesselAMI2D-topologyChange",
)
TIMEOUT = 300


def find_block_span(text: str, name: str) -> tuple[int, int, int]:
    match = re.search(
        rf"(?m)^[ \t]*{re.escape(name)}[ \t]*(?:\n[ \t]*)?\{{",
        text,
    )
    if not match:
        raise ValueError(f"could not find dictionary block {name}")
    opening = text.find("{", match.start())
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return match.start(), opening, index
    raise ValueError(f"unterminated dictionary block {name}")


def extract_block(text: str, name: str) -> str:
    _, opening, closing = find_block_span(text, name)
    return text[opening + 1:closing]


def configure_ami_cache(
    case_dir: Path,
    cache_size: int,
    rotation_axis: list[float],
    rotation_centre: list[float],
    n_theta_stencil_max: int = 2,
) -> dict:
    mesh_path = case_dir / "system/blockMeshDict"
    mesh = mesh_path.read_text(encoding="utf-8")
    spans = {
        name: find_block_span(mesh, name)
        for name in ("AMI1", "AMI2")
    }
    blocks = {}
    for name, (start, opening, closing) in spans.items():
        block = mesh[start:closing + 1]
        local_opening = opening - start
        for key, value in (
            ("cacheSize", str(cache_size)),
            ("nThetaStencilMax", str(n_theta_stencil_max)),
            ("rotationAxis", runTests.fmt_vector(rotation_axis)),
            ("rotationCentre", runTests.fmt_vector(rotation_centre)),
        ):
            entry = re.compile(
                rf"(?m)^([ \t]*){re.escape(key)}[ \t]+[^;]+;"
            )
            if entry.search(block):
                block, count = entry.subn(
                    rf"\g<1>{key} {value};",
                    block,
                    count=1,
                )
                if count != 1:
                    raise ValueError(f"could not set {name} {key}")
            else:
                block = (
                    block[:local_opening + 1]
                    + f"\n        {key} {value};"
                    + block[local_opening + 1:]
                )
        blocks[name] = block

    first_name = min(spans, key=lambda name: spans[name][0])
    second_name = "AMI2" if first_name == "AMI1" else "AMI1"
    first_start, _, first_end = spans[first_name]
    second_start, _, second_end = spans[second_name]
    between = mesh[first_end + 1:second_start]
    mesh = (
        mesh[:first_start]
        + blocks["AMI2"]
        + between
        + blocks["AMI1"]
        + mesh[second_end + 1:]
    )
    if (
        find_block_span(mesh, "AMI2")[0]
        >= find_block_span(mesh, "AMI1")[0]
    ):
        raise ValueError("moving AMI2 must own the cached interpolation")
    mesh_path.write_text(mesh, encoding="utf-8")
    return {
        "cache_size": cache_size,
        "n_theta_stencil_max": n_theta_stencil_max,
        "cache_owner_patch": "AMI2",
        "rotation_axis": rotation_axis,
        "rotation_centre": rotation_centre,
    }


def configure_decomposition(case_dir: Path, ranks: int) -> tuple[int, int, int]:
    decompositions = {
        2: (2, 1, 1),
        4: (2, 2, 1),
    }
    if ranks not in decompositions:
        raise ValueError(f"unsupported T19 MPI rank count {ranks}")

    path = case_dir / "system/decomposeParDict"
    contents = path.read_text(encoding="utf-8")
    contents, rank_count = re.subn(
        r"(?m)^([ \t]*numberOfSubdomains[ \t]+)\d+([ \t]*;)",
        rf"\g<1>{ranks}\g<2>",
        contents,
        count=1,
    )
    if rank_count != 1:
        raise ValueError("could not set numberOfSubdomains")
    coeffs = " ".join(str(value) for value in decompositions[ranks])
    contents, coeff_count = re.subn(
        r"(?m)^([ \t]*n[ \t]+)\([^)]*\)([ \t]*;)",
        rf"\g<1>({coeffs})\g<2>",
        contents,
        count=1,
    )
    if coeff_count != 1:
        raise ValueError("could not set simpleCoeffs n")
    path.write_text(contents, encoding="utf-8")
    return decompositions[ranks]


def parse_scalar_field(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8")
    internal = text.split("internalField", 1)[1].split("boundaryField", 1)[0]
    uniform = re.search(
        r"\buniform\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        internal,
    )
    if uniform:
        values = [float(uniform.group(1))]
    else:
        match = re.search(
            r"nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
            internal,
            re.S,
        )
        if not match:
            raise ValueError(f"could not parse internal scalar field {path}")
        values = [
            float(value)
            for value in re.findall(
                r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
                match.group(2),
            )
        ]
        if len(values) != int(match.group(1)):
            raise ValueError(f"wrong scalar count in {path}")
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"non-finite scalar in {path}")
    return values


def parse_vector_field(path: Path) -> list[tuple[float, float, float]]:
    text = path.read_text(encoding="utf-8")
    internal = text.split("internalField", 1)[1].split("boundaryField", 1)[0]
    match = re.search(
        r"nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
        internal,
        re.S,
    )
    if not match:
        raise ValueError(f"could not parse internal vector field {path}")
    vectors = []
    for vector_match in re.finditer(r"\(([^()]*)\)", match.group(2)):
        components = [
            float(value)
            for value in vector_match.group(1).split()
        ]
        if len(components) != 3:
            raise ValueError(f"invalid vector in {path}")
        vectors.append(tuple(components))
    if len(vectors) != int(match.group(1)):
        raise ValueError(f"wrong vector count in {path}")
    if not all(math.isfinite(value) for vector in vectors for value in vector):
        raise ValueError(f"non-finite vector in {path}")
    return vectors


def parse_boundary_scalar_field(path: Path, patch_name: str) -> list[float]:
    text = path.read_text(encoding="utf-8")
    boundary = extract_block(text, "boundaryField")
    block = extract_block(boundary, patch_name)
    uniform = re.search(
        r"\bvalue\s+uniform\s+"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        block,
    )
    if uniform:
        return [float(uniform.group(1))]
    match = re.search(
        r"value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        block,
        re.S,
    )
    if not match:
        raise ValueError(f"could not parse {patch_name} values in {path}")
    values = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
            match.group(2),
        )
    ]
    if len(values) != int(match.group(1)):
        raise ValueError(f"wrong {patch_name} value count in {path}")
    return values


def parse_label_list(path: Path) -> list[int]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"\n\s*(\d+)\s*\(\s*([^)]*)\)\s*;?", text, re.S)
    if not match:
        raise ValueError(f"could not parse label list {path}")
    labels = [int(value) for value in re.findall(r"\d+", match.group(2))]
    if len(labels) != int(match.group(1)):
        raise ValueError(f"wrong label count in {path}")
    return labels


def patch_range(case_dir: Path, patch_name: str) -> tuple[int, int]:
    block = extract_block(
        (case_dir / "constant/polyMesh/boundary").read_text(
            encoding="utf-8"
        ),
        patch_name,
    )
    start = re.search(r"\bstartFace\s+(\d+)\s*;", block)
    size = re.search(r"\b(?:nFaces|size)\s+(\d+)\s*;", block)
    if not start or not size:
        raise ValueError(f"missing face range for patch {patch_name}")
    return int(start.group(1)), int(size.group(1))


def latest_time(case_dir: Path) -> str:
    times = []
    for path in case_dir.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0 and math.isfinite(value):
            times.append((value, path.name))
    if not times:
        raise ValueError(f"no written solution time in {case_dir}")
    return max(times)[1]


def read_profile_fields(
    case_dir: Path,
    field_time: str,
) -> tuple[list[float], list[tuple[float, float, float]]]:
    temperature = parse_scalar_field(case_dir / field_time / "T")
    centres = parse_vector_field(case_dir / field_time / "C")
    if len(temperature) == 1 and len(centres) > 1:
        temperature *= len(centres)
    if len(temperature) != len(centres):
        raise ValueError("T and C internal field sizes differ")
    return temperature, centres


def zone_means(
    values: list[float],
    zones: dict[str, list[int]],
) -> dict[str, float]:
    means = {}
    for zone, labels in zones.items():
        if not labels or any(label < 0 or label >= len(values) for label in labels):
            raise ValueError(f"invalid cell labels for zone {zone}")
        means[zone] = sum(values[label] for label in labels)/len(labels)
    return means


def profile_metrics(
    case_dir: Path,
    test_id: str,
    values: dict,
    field_time: str,
) -> dict:
    zones_path = case_dir / "constant/polyMesh/cellZones"
    zones = runTests.parse_cell_zones(zones_path.read_text(encoding="utf-8"))
    if set(zones) != {"left", "right"}:
        raise ValueError("expected cell zones left and right")

    temperature, centres = read_profile_fields(case_dir, field_time)

    metrics = {
        "temperature_min": min(temperature),
        "temperature_max": max(temperature),
        "zone_mean_temperature": zone_means(temperature, zones),
        "cell_count": len(temperature),
    }

    if test_id in {"T14", "T18"}:
        if test_id == "T14":
            metrics["analytical_reference"] = solid_reference_metrics(
                test_id,
                values,
            )
        metrics["uniform_temperature_error"] = max(
            abs(value - values["initialT"]) for value in temperature
        )

    if test_id in {"T15", "T17", "T19"}:
        reference = solid_reference_metrics(test_id, values)
        if test_id == "T17":
            conductivities = (
                reference["conductivity_tensors"]["left"][0][0],
                reference["conductivity_tensors"]["right"][0][0],
            )
            heat_rate = series_resistance_heat_rate(
                values["referenceArea"],
                (values["leftLength"], values["rightLength"]),
                conductivities,
                values["hotT"],
                values["coldT"],
            )
            reference["analytical_heat_rate"] = heat_rate
            reference["interface_temperature"] = (
                series_resistance_interface_temperature(
                    values["hotT"],
                    heat_rate,
                    values["referenceArea"],
                    values["leftLength"],
                    conductivities[0],
                )
            )
            reference["temperature_scale"] = abs(
                values["hotT"] - values["coldT"]
            )
            reference["heat_rate_scale"] = abs(heat_rate)
        else:
            conductivities = (
                values["leftKappa"][0],
                values["rightKappa"][0],
            )
        area = values["referenceArea"]
        heat_rate = reference["analytical_heat_rate"]
        heat_flux = heat_rate/area
        junction_temperature = reference["interface_temperature"]
        expected = []
        for centre in centres:
            x = centre[0]
            if x <= values["leftLength"]:
                expected.append(
                    values["hotT"]
                    - heat_flux*x/conductivities[0]
                )
            else:
                expected.append(
                    junction_temperature
                    - heat_flux
                    *(x - values["leftLength"])
                    / conductivities[1]
                )
        metrics["analytical_reference"] = reference
        metrics["normalized_profile_error"] = (
            max(
                abs(actual - target)
                for actual, target in zip(temperature, expected)
            )
            / max(abs(values["hotT"] - values["coldT"]), 1e-300)
        )
    if test_id == "T16":
        cell_volumes = {
            "left": (
                values["leftLength"]
                * values["referenceArea"]
                /(2*values["leftNy"])
            ),
            "right": (
                values["rightLength"]
                * values["referenceArea"]
                /(2*values["rightNy"])
            ),
        }
        energy_change = 0.0
        for zone in ("left", "right"):
            volumetric_capacity = (
                values[f"{zone}Rho"]*values[f"{zone}Cp"]
            )
            energy_change += sum(
                volumetric_capacity
                *cell_volumes[zone]
                *(temperature[label] - values["initialT"])
                for label in zones[zone]
            )
        metrics["stored_energy_change"] = energy_change
        metrics["volumetric_heat_capacity"] = {
            zone: values[f"{zone}Rho"]*values[f"{zone}Cp"]
            for zone in ("left", "right")
        }

    return metrics


def interface_heat_rate(
    case_dir: Path,
    values: dict,
    temperature: list[float],
    centres: list[tuple[float, float, float]],
    field_time: str,
    normal_conductivities: dict[str, float] | None = None,
) -> dict:
    owner = parse_label_list(case_dir / "constant/polyMesh/owner")
    temperature_path = case_dir / field_time / "T"
    rates = {}
    for patch, zone, sign in (
        ("AMI1", "left", 1),
        ("AMI2", "right", -1),
    ):
        start_face, face_count = patch_range(case_dir, patch)
        patch_values = parse_boundary_scalar_field(temperature_path, patch)
        if len(patch_values) == 1:
            patch_values *= face_count
        if len(patch_values) != face_count:
            raise ValueError(f"wrong face count for {patch}")

        conductivity = (
            normal_conductivities[zone]
            if normal_conductivities is not None
            else values[f"{zone}Kappa"][0]
        )
        face_area = values["referenceArea"]/face_count
        total = 0.0
        for local_face in range(face_count):
            cell = owner[start_face + local_face]
            distance = abs(centres[cell][0] - values["leftLength"])
            if distance <= 0:
                raise ValueError(f"zero owner-to-face distance on {patch}")
            delta = (
                temperature[cell] - patch_values[local_face]
                if sign > 0
                else patch_values[local_face] - temperature[cell]
            )
            total += conductivity*delta/distance*face_area
        rates[patch] = total
    return {
        "owner_to_neighbour": rates,
        "absolute_imbalance": abs(
            rates["AMI1"] - rates["AMI2"]
        ),
        "normalization": max(abs(rates["AMI1"]), abs(rates["AMI2"]), 1e-300),
        "relative_imbalance": (
            abs(rates["AMI1"] - rates["AMI2"])
            / max(abs(rates["AMI1"]), abs(rates["AMI2"]), 1e-300)
        ),
    }


def conformal_interface_heat_rate(
    case_dir: Path,
    values: dict,
    temperature: list[float],
    centres: list[tuple[float, float, float]],
) -> float:
    zones = runTests.parse_cell_zones(
        (case_dir / "constant/polyMesh/cellZones").read_text(
            encoding="utf-8"
        )
    )
    left_cells = [
        cell for cell in zones["left"]
        if centres[cell][0] < values["leftLength"]
    ]
    right_cells = [
        cell for cell in zones["right"]
        if centres[cell][0] > values["leftLength"]
    ]
    left_x = max(centres[cell][0] for cell in left_cells)
    right_x = min(centres[cell][0] for cell in right_cells)
    left_cells = [cell for cell in left_cells if abs(centres[cell][0] - left_x) < 1e-12]
    right_cells = [cell for cell in right_cells if abs(centres[cell][0] - right_x) < 1e-12]
    total = 0.0
    used = set()
    for left_cell in left_cells:
        match = min(
            right_cells,
            key=lambda cell: abs(centres[cell][1] - centres[left_cell][1]),
        )
        if match in used or abs(
            centres[match][1] - centres[left_cell][1]
        ) > 1e-12:
            raise ValueError("could not pair conformal interface cells")
        used.add(match)
        resistance = (
            (values["leftLength"] - left_x)/values["leftKappa"][0]
            + (right_x - values["leftLength"])/values["rightKappa"][0]
        )
        flux = (temperature[left_cell] - temperature[match])/resistance
        total += flux*values["referenceArea"]/len(left_cells)
    if len(used) != len(right_cells):
        raise ValueError("unpaired conformal interface cells")
    return total


def modify_conformal_case(case_dir: Path) -> None:
    mesh_path = case_dir / "system/blockMeshDict"
    mesh = mesh_path.read_text(encoding="utf-8")
    if mesh.count("type cyclicAMI;") != 2:
        raise ValueError("expected two cyclicAMI patches in solid template")
    mesh = mesh.replace("type cyclicAMI;", "type patch;")
    mesh = re.sub(
        r"(?m)^[ \t]*(?:neighbourPatch|transform) [^;]+;\n",
        "",
        mesh,
    )
    if "mergePatchPairs ();" not in mesh:
        raise ValueError("could not find mergePatchPairs in solid template")
    mesh = mesh.replace(
        "mergePatchPairs ();",
        "mergePatchPairs ((AMI1 AMI2));",
    )
    mesh_path.write_text(mesh, encoding="utf-8")

    field_block = re.compile(
        r"(?m)^[ \t]*AMI[12][ \t]*(?:\n[ \t]*)?\{\n"
        r"(?:[ \t]+[^\n]*\n)*?"
        r"[ \t]*\}\n"
    )
    for field in (case_dir / "0/T", case_dir / "0/p"):
        contents = field.read_text(encoding="utf-8")
        contents, count = field_block.subn("", contents)
        if count != 2:
            raise ValueError(f"could not remove cyclic fields from {field}")
        field.write_text(contents, encoding="utf-8")


def configure_solid_case(
    case_dir: Path,
    values: dict,
    mode: str,
    dynamic_ami: bool = False,
    dynamic_topology: bool = True,
) -> None:
    control_path = case_dir / "system/controlDict"
    control = control_path.read_text(encoding="utf-8")
    end_time = values["endTime"]
    delta_t = values["deltaT"]
    if mode == "steady":
        end_time = 100
        schemes_path = case_dir / "system/fvSchemes"
        schemes = schemes_path.read_text(encoding="utf-8")
        schemes = schemes.replace("default Euler;", "default steadyState;", 1)
        schemes_path.write_text(schemes, encoding="utf-8")
    else:
        solver = """solver GAMG;
        smoother GaussSeidel;
        tolerance 1e-12;
        relTol 0;
        maxIter 100;"""
        solution = f"""FoamFile
{{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSolution;
}}
solvers
{{
    h
    {{
        {solver}
    }}
    hFinal
    {{
        {solver}
    }}
}}
PIMPLE
{{
    nOuterCorrectors 1;
    nCorrectors 1;
    nNonOrthogonalCorrectors 0;
    meshFluxCorrection {'true' if dynamic_ami else 'false'};
}}
relaxationFactors
{{
    equations
    {{
        h 1;
    }}
}}
"""
        (case_dir / "system/fvSolution").write_text(
            solution,
            encoding="utf-8",
        )

    control = re.sub(
        r"(?m)^endTime\s+[^;]+;",
        f"endTime {end_time};",
        control,
        count=1,
    )
    control = re.sub(
        r"(?m)^deltaT\s+[^;]+;",
        f"deltaT {delta_t};",
        control,
        count=1,
    )
    control_path.write_text(control, encoding="utf-8")

    if dynamic_ami:
        motion_function = values.get("motionFunction", "rotatingMotion")
        motion_origin = runTests.fmt_vector(
            values.get("motionOrigin", [1, 0.5, 0.05])
        )
        motion_axis = runTests.fmt_vector(
            values.get("motionAxis", [0, 0, 1])
        )
        motion_omega = runTests.fmt_number(
            values.get("motionOmega", 1)
        )
        if motion_function == "rotatingMotion":
            motion_parameter = f"axis {motion_axis};"
        elif motion_function == "oscillatingRotatingMotion":
            motion_parameter = (
                "amplitude "
                + runTests.fmt_vector(values["motionAmplitude"])
                + ";"
            )
        else:
            raise ValueError(
                f"unsupported solid motion function {motion_function}"
            )
        dynamic_mesh_type = (
            "dynamicMotionSolverFvMeshAMI"
            if dynamic_topology
            else "dynamicMotionSolverFvMesh"
        )
        dynamic_mesh = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object dynamicMeshDict;
}
dynamicFvMesh @DYNAMIC_MESH_TYPE@;
motionSolver solidBody;
cellZone right;
solidBodyMotionFunction @MOTION_FUNCTION@;
origin @ORIGIN@;
@MOTION_PARAMETER@
omega @OMEGA@;
"""
        dynamic_mesh = (
            dynamic_mesh
            .replace("@DYNAMIC_MESH_TYPE@", dynamic_mesh_type)
            .replace("@MOTION_FUNCTION@", motion_function)
            .replace("@ORIGIN@", motion_origin)
            .replace("@MOTION_PARAMETER@", motion_parameter)
            .replace("@OMEGA@", motion_omega)
        )
        (case_dir / "constant/dynamicMeshDict").write_text(
            dynamic_mesh,
            encoding="utf-8",
        )


def run_command(
    command: list[str],
    case_dir: Path,
    env: dict[str, str],
    timeout: float,
    log_name: str,
) -> dict:
    log_path = case_dir / log_name
    result = runTests.run_profiled_process(
        command,
        case_dir,
        env,
        timeout,
        case_dir / f"time.{log_name}.txt",
    )
    runTests.write_log(log_path, result)
    return {
        "command": command,
        "return_code": result["return_code"],
        "timed_out": result["timed_out"],
        "duration_seconds": result["duration_seconds"],
        "log": str(log_path),
        "performance": result["performance"],
    }


def solver_warnings(solver_log: str) -> list[str]:
    lines = solver_log.splitlines()
    warnings = []
    index = 0
    while index < len(lines):
        if not lines[index].startswith("--> FOAM Warning :"):
            index += 1
            continue

        block = [lines[index]]
        index += 1
        while (
            index < len(lines)
            and lines[index]
            and lines[index][0].isspace()
            and not lines[index].lstrip().startswith("--> FOAM Warning :")
        ):
            block.append(lines[index])
            index += 1
        warnings.append("\n".join(block))
    return warnings


def collect_profile(
    case_dir: Path,
    test_id: str,
    values: dict,
    time_name: str,
) -> dict:
    metrics = profile_metrics(case_dir, test_id, values, time_name)
    if test_id in {"T14", "T15", "T17", "T19"} and (
        case_dir / "constant/polyMesh/boundary"
    ).is_file():
        boundary = (case_dir / "constant/polyMesh/boundary").read_text(
            encoding="utf-8"
        )
        if re.search(
            r"(?m)^[ \t]*AMI1[ \t]*(?:\n[ \t]*)?\{",
            boundary,
        ):
            temperature, centres = read_profile_fields(case_dir, time_name)
            normal_conductivities = None
            if test_id == "T17":
                tensors = metrics["analytical_reference"][
                    "conductivity_tensors"
                ]
                normal_conductivities = {
                    zone: tensors[zone][0][0]
                    for zone in ("left", "right")
                }
            metrics["interface_heat_rate"] = interface_heat_rate(
                case_dir,
                values,
                temperature,
                centres,
                time_name,
                normal_conductivities,
            )
            if test_id in {"T15", "T17", "T19"}:
                expected = metrics["analytical_reference"][
                    "analytical_heat_rate"
                ]
                metrics["interface_heat_rate"][
                    "analytical_heat_rate_error"
                ] = {
                    patch: abs(rate - expected)
                    / max(abs(expected), 1e-300)
                    for patch, rate in metrics["interface_heat_rate"][
                        "owner_to_neighbour"
                    ].items()
                }
        else:
            temperature, centres = read_profile_fields(case_dir, time_name)
            metrics["conformal_interface_heat_rate"] = (
                conformal_interface_heat_rate(
                    case_dir,
                    values,
                    temperature,
                    centres,
                )
            )
    return metrics


def compare_serial_mpi_fields(
    serial_record: dict,
    parallel_record: dict,
    reference_scale: float,
) -> dict:
    if (
        serial_record["status"] != "runtime_pass"
        or parallel_record["status"] != "runtime_pass"
    ):
        return {
            "status": "not_comparable",
            "reason": "one or both solidFoam runs did not complete",
        }

    try:
        serial_temperature, serial_centres = read_profile_fields(
            Path(serial_record["case_directory"]),
            serial_record["latest_time"],
        )
        parallel_temperature, parallel_centres = read_profile_fields(
            Path(parallel_record["case_directory"]),
            parallel_record["latest_time"],
        )
    except (OSError, ValueError) as error:
        return {"status": "error", "error": str(error)}

    def field_map(temperature, centres):
        result = {}
        for centre, value in zip(centres, temperature):
            key = tuple(round(component, 12) for component in centre)
            if key in result:
                raise ValueError(f"duplicate cell centre {key}")
            result[key] = value
        return result

    try:
        serial = field_map(serial_temperature, serial_centres)
        parallel = field_map(parallel_temperature, parallel_centres)
    except ValueError as error:
        return {"status": "error", "error": str(error)}
    if serial.keys() != parallel.keys():
        return {
            "status": "error",
            "error": "serial and MPI cell-centre sets differ",
            "serial_cell_count": len(serial),
            "parallel_cell_count": len(parallel),
        }
    differences = [
        abs(serial[centre] - parallel[centre])
        for centre in serial
    ]
    return {
        "status": "compared",
        "cells_compared": len(differences),
        "max_abs_temperature_difference": max(differences),
        "normalized_max_abs_temperature_difference": (
            max(differences)/max(abs(reference_scale), 1e-300)
        ),
        "mean_abs_temperature_difference": (
            math.fsum(differences)/len(differences)
        ),
    }


def run_solid_case(
    test_id: str,
    values: dict,
    case_dir: Path,
    env: dict[str, str],
    timeout: float,
    target_id: str,
    conformal: bool = False,
    parallel: bool = False,
    dynamic_ami: bool = False,
    dynamic_topology: bool = True,
    cache_ami: bool = False,
    record_id: str | None = None,
) -> dict:
    template = TEST_DIR / "cases/solidFoam/template"
    runTests.render_solid_case(template, case_dir, values)
    mode = values.get("solverMode", "steady")
    configure_solid_case(
        case_dir,
        values,
        mode,
        dynamic_ami,
        dynamic_topology,
    )
    if conformal:
        modify_conformal_case(case_dir)
    if cache_ami and not dynamic_ami:
        raise ValueError("AMI caching is only enabled for the moving T18 screen")
    if cache_ami and dynamic_topology:
        raise ValueError("AMI caching cannot be combined with face creation")

    parallel_ranks = int(values.get("decompose", 2)) if parallel else 1
    decomposition = (
        configure_decomposition(case_dir, parallel_ranks)
        if parallel
        else None
    )
    cache_configuration = (
        configure_ami_cache(
            case_dir,
            int(values.get("cacheSize", 4)),
            values.get("motionAxis", [0, 0, 1]),
            values.get("motionOrigin", [1, 0.5, 0.05]),
            int(values.get("nThetaStencilMax", 2)),
        )
        if cache_ami
        else None
    )

    foam_dictionary = shutil.which("foamDictionary", path=env.get("PATH"))
    block_mesh = shutil.which("blockMesh", path=env.get("PATH"))
    check_mesh = shutil.which("checkMesh", path=env.get("PATH"))
    solid_foam = shutil.which("solidFoam", path=env.get("PATH"))
    if not all((foam_dictionary, block_mesh, check_mesh, solid_foam)):
        raise ValueError(f"required {target_id} solidFoam tools are missing")

    properties_ok, properties, property_errors = (
        runTests.validate_material_properties(
            foam_dictionary,
            case_dir,
            values,
            env,
            timeout,
            case_dir / "log.materialProperties",
        )
    )
    mixture = runTests.run_process(
        [
            foam_dictionary,
            "-entry",
            "thermoType.mixture",
            "-value",
            str(case_dir / "constant/thermophysicalProperties"),
        ],
        case_dir,
        env,
        timeout,
    )
    runTests.write_log(case_dir / "log.foamDictionary", mixture)
    mixture_name = mixture["output"].strip()

    commands = {}
    mesh_result = run_command(
        [block_mesh, "-case", str(case_dir)],
        case_dir,
        env,
        timeout,
        "log.blockMesh",
    )
    commands["blockMesh"] = mesh_result
    check_result = None
    if mesh_result["return_code"] == 0 and not mesh_result["timed_out"]:
        check_result = run_command(
            [check_mesh, "-case", str(case_dir)],
            case_dir,
            env,
            timeout,
            "log.checkMesh",
        )
        commands["checkMesh"] = check_result

    decompose_result = None
    solver_result = None
    reconstruct_result = None
    debug_switches = (
        ["-debug-switch", "AMICache=1"] if cache_ami else []
    )
    if (
        mesh_result["return_code"] == 0
        and not mesh_result["timed_out"]
        and check_result
        and check_result["return_code"] == 0
        and not check_result["timed_out"]
    ):
        if parallel:
            decompose_par = shutil.which("decomposePar", path=env.get("PATH"))
            mpirun = shutil.which("mpirun", path=env.get("PATH"))
            if not decompose_par or not mpirun:
                raise ValueError("decomposePar and mpirun are required for T19")
            decompose_result = run_command(
                [decompose_par, "-case", str(case_dir), "-force"],
                case_dir,
                env,
                timeout,
                "log.decomposePar",
            )
            commands["decomposePar"] = decompose_result
            if decompose_result["return_code"] == 0 and not decompose_result["timed_out"]:
                solver_command = [
                    mpirun,
                    "-np",
                    str(parallel_ranks),
                    solid_foam,
                    *debug_switches,
                    "-case",
                    str(case_dir),
                    "-parallel",
                ]
                solver_result = run_command(
                    solver_command,
                    case_dir,
                    env,
                    timeout,
                    "log.solidFoam",
                )
                commands["solidFoam"] = solver_result
                if solver_result["return_code"] == 0 and not solver_result["timed_out"]:
                    reconstruct_par = shutil.which(
                        "reconstructPar",
                        path=env.get("PATH"),
                    )
                    if not reconstruct_par:
                        raise ValueError("reconstructPar is required for T19")
                    reconstruct_result = run_command(
                        [reconstruct_par, "-case", str(case_dir)],
                        case_dir,
                        env,
                        timeout,
                        "log.reconstructPar",
                    )
                    commands["reconstructPar"] = reconstruct_result
        else:
            solver_result = run_command(
                [solid_foam, *debug_switches, "-case", str(case_dir)],
                case_dir,
                env,
                timeout,
                "log.solidFoam",
            )
            commands["solidFoam"] = solver_result

    solver_log = ""
    if solver_result and Path(solver_result["log"]).is_file():
        solver_log = Path(solver_result["log"]).read_text(encoding="utf-8")
    runtime_selected = "pureZoneMixture" in solver_log
    dynamic_mesh_type = (
        "dynamicMotionSolverFvMeshAMI"
        if dynamic_topology
        else "dynamicMotionSolverFvMesh"
    )
    dynamic_mesh_selected = bool(
        re.search(
            rf"(?m)^Selecting dynamicFvMesh {dynamic_mesh_type}$",
            solver_log,
        )
    )
    ami_rebuild_events = solver_log.count("AMI: Creating AMI")
    additional_face_records = [
        {
            "patch": match.group(1),
            "count": int(match.group(2)),
        }
        for match in re.finditer(
            r"(?m)^AMI: Patch (\S+) additional faces:\s*(\d+)\s*$",
            solver_log,
        )
    ]
    ami_created_faces = sum(
        record["count"] for record in additional_face_records
    )
    cache_restore_events = solver_log.count(
        "AMI: Cached weights and addresses restored"
    )
    cache_unavailable_events = solver_log.count(
        "AMI: Cached weights and addresses unavailable"
    )
    cache_owner_source_events = len(
        re.findall(
            r"AMI: Creating AMI for source:\s*AMI2 and target:\s*AMI1",
            solver_log,
        )
    )
    cache_interpolation_records = []
    gamg_solver_times = []
    current_time = None
    cache_record_pattern = re.compile(
        r"-- t0:\s*([^\s]+)\s+theta:\s*([^\s]+)"
        r"(?:\s+t1:\s*([^\s]+))?"
        r"\s+i0:\s*(-?\d+)\s+i1:\s*(-?\d+)\s+w:\s*([^\s]+)"
    )
    for line in solver_log.splitlines():
        time_match = re.search(r"\bTime\s*=\s*([^\s]+)", line)
        if time_match:
            current_time = float(time_match.group(1))

        cache_match = cache_record_pattern.search(line)
        if cache_match:
            cache_interpolation_records.append(
                {
                    "time": current_time,
                    "theta0": float(cache_match.group(1)),
                    "theta": float(cache_match.group(2)),
                    "theta1": (
                        float(cache_match.group(3))
                        if cache_match.group(3) is not None
                        else None
                    ),
                    "index0": int(cache_match.group(4)),
                    "index1": int(cache_match.group(5)),
                    "weight": float(cache_match.group(6)),
                }
            )

        if re.match(r"^\s*GAMG:\s+Solving for h\b", line):
            gamg_solver_times.append(current_time)

    gamg_solver_events = len(gamg_solver_times)
    fractional_cache_records = [
        record
        for record in cache_interpolation_records
        if record["index1"] != -1 and 0 < record["weight"] < 1
    ]
    fractional_cache_gamg_times = [
        record["time"]
        for record in fractional_cache_records
        if record["time"] is not None
        and any
        (
            solver_time is not None
            and math.isclose
            (
                record["time"],
                solver_time,
                rel_tol=0,
                abs_tol=1e-12,
            )
            for solver_time in gamg_solver_times
        )
    ]
    require_gamg = bool(values.get("requireGAMG", test_id == "T19"))
    require_fractional_cache = bool(
        values.get("requireFractionalCache", False)
    )
    parallel_completed = (
        not parallel
        or (
            reconstruct_result is not None
            and reconstruct_result["return_code"] == 0
            and not reconstruct_result["timed_out"]
        )
    )
    solver_completed = (
        solver_result is not None
        and solver_result["return_code"] == 0
        and not solver_result["timed_out"]
        and "End" in solver_log
        and runtime_selected
        and properties_ok
        and mixture["return_code"] == 0
        and not mixture["timed_out"]
        and parallel_completed
        and (not require_gamg or gamg_solver_events > 0)
        and (
            not require_fractional_cache
            or bool(fractional_cache_records)
        )
        and (
            not (require_fractional_cache and require_gamg)
            or bool(fractional_cache_gamg_times)
        )
        and (
            not dynamic_ami
            or (
                dynamic_mesh_selected
                and ami_rebuild_events >= 2
                and (
                    not dynamic_topology
                    or ami_created_faces > 0
                )
                and (
                    dynamic_topology
                    or ami_created_faces == 0
                )
            )
        )
        and (
            not cache_ami
            or (
                cache_restore_events > 0
                and cache_owner_source_events > 0
                and len(cache_interpolation_records)
                == cache_restore_events
            )
        )
    )
    latest = None
    field_metrics = None
    metrics_command = None
    if solver_completed:
        try:
            latest = latest_time(case_dir)
            post_process = shutil.which("postProcess", path=env.get("PATH"))
            if post_process:
                metrics_command = run_command(
                    [
                        post_process,
                        "-case",
                        str(case_dir),
                        "-func",
                        "writeCellCentres",
                        "-time",
                        latest,
                    ],
                    case_dir,
                    env,
                    timeout,
                    "log.writeCellCentres",
                )
                commands["writeCellCentres"] = metrics_command
                if (
                    metrics_command["return_code"] == 0
                    and not metrics_command["timed_out"]
                ):
                    field_metrics = collect_profile(
                        case_dir,
                        "T15" if conformal else test_id,
                        values,
                        latest,
                    )
        except (OSError, ValueError) as error:
            field_metrics = {"error": str(error)}

    warnings = solver_warnings(solver_log)
    cell_zone_path = case_dir / "constant/polyMesh/cellZones"
    initial_cell_zones = (
        runTests.parse_cell_zones(
            cell_zone_path.read_text(encoding="utf-8")
        )
        if cell_zone_path.is_file()
        else {}
    )
    cell_zone_counts = {
        zone: len(labels) for zone, labels in initial_cell_zones.items()
    }
    updated_cell_zone_path = (
        case_dir / latest / "polyMesh/cellZones"
        if latest is not None
        else None
    )
    if not solver_completed:
        cell_zone_update_status = "not_run"
        updated_cell_zone_counts = None
    elif (
        updated_cell_zone_path is not None
        and updated_cell_zone_path.is_file()
    ):
        updated_cell_zones = runTests.parse_cell_zones(
            updated_cell_zone_path.read_text(encoding="utf-8")
        )
        updated_cell_zone_counts = {
            zone: len(labels) for zone, labels in updated_cell_zones.items()
        }
        cell_zone_update_status = (
            "unchanged"
            if updated_cell_zones == initial_cell_zones
            else "changed"
        )
    else:
        updated_cell_zone_counts = None
        cell_zone_update_status = (
            "no_topology_output"
            if dynamic_ami and not dynamic_topology
            else "not_written"
        )
    if dynamic_ami and dynamic_topology:
        solver_completed = (
            solver_completed
            and cell_zone_update_status == "unchanged"
        )
    record_id = record_id or (
        f"{test_id}-conformal" if conformal else test_id
    )
    return {
        "id": record_id,
        "variant": values.get("variant", "default"),
        "solver_mode": mode,
        "rank_count": parallel_ranks,
        "decomposition": decomposition,
        "status": "runtime_pass" if solver_completed else "runtime_failure",
        "mixture_dictionary": mixture_name,
        "mixture_dictionary_valid": (
            mixture["return_code"] == 0
            and not mixture["timed_out"]
            and mixture_name == "pureZoneMixture"
        ),
        "runtime_mixture_selected": runtime_selected,
        "dynamic_mesh_selected": (
            dynamic_mesh_selected if dynamic_ami else None
        ),
        "dynamic_mesh_type": (
            dynamic_mesh_type if dynamic_ami else None
        ),
        "dynamic_topology": dynamic_topology if dynamic_ami else None,
        "ami_rebuild_events": ami_rebuild_events,
        "ami_additional_face_records": additional_face_records,
        "ami_created_faces": ami_created_faces,
        "ami_cache_configuration": cache_configuration,
        "cache_restore_events": cache_restore_events,
        "cache_unavailable_events": cache_unavailable_events,
        "cache_owner_source_events": cache_owner_source_events,
        "cache_interpolation_records": cache_interpolation_records,
        "fractional_cache_records": fractional_cache_records,
        "fractional_cache_gamg_times": fractional_cache_gamg_times,
        "cache_fraction_verified": (
            bool(fractional_cache_records)
            if cache_ami
            else None
        ),
        "cache_fraction_gamg_verified": (
            bool(fractional_cache_gamg_times)
            if cache_ami and require_gamg
            else None
        ),
        "cache_reuse_verified": (
            cache_restore_events > 0
            and cache_owner_source_events > 0
            and len(cache_interpolation_records) == cache_restore_events
            if cache_ami
            else None
        ),
        "linear_solver": "GAMG" if gamg_solver_events else None,
        "gamg_solver_events": gamg_solver_events,
        "material_properties_valid": properties_ok,
        "material_properties": properties,
        "material_property_errors": property_errors,
        "cell_zone_counts": cell_zone_counts,
        "cell_zone_update_status": cell_zone_update_status,
        "updated_cell_zone_counts": updated_cell_zone_counts,
        "latest_time": latest,
        "field_metrics": field_metrics,
        "case_directory": str(case_dir),
        "solver_warnings": warnings,
        "commands": commands,
        "solver_log": (
            solver_result["log"] if solver_result else None
        ),
        "solver_attempted": solver_result is not None,
    }


def configure_mixer_copy(source: Path, case_dir: Path) -> list[dict]:
    if not source.is_dir():
        raise FileNotFoundError(f"tutorial case not found: {source}")
    source_files = []
    for path in sorted(source.rglob("*")):
        if path.is_file():
            source_files.append(
                {
                    "path": str(path.relative_to(source)),
                    "size_bytes": path.stat().st_size,
                    "sha256": runTests.sha256_file(path),
                }
            )
    shutil.copytree(source, case_dir)
    control_path = case_dir / "system/controlDict"
    control = control_path.read_text(encoding="utf-8")
    for entry, value in (("endTime", "0.002"), ("writeInterval", "0.001")):
        pattern = rf"(?m)^[ \t]*{entry}\s+[^;]+;"
        if not re.search(pattern, control):
            raise ValueError(f"{entry} not found in {control_path}")
        control = re.sub(
            pattern,
            f"{entry} {value};",
            control,
            count=1,
        )
    control_path.write_text(control, encoding="utf-8")
    return source_files


def run_mixer_case(
    case_name: str,
    source: Path,
    case_dir: Path,
    env: dict[str, str],
    timeout: float,
) -> dict:
    source_files = configure_mixer_copy(source, case_dir)
    run_record = run_command(
        ["sh", "-e", "./Allrun"],
        case_dir,
        env,
        timeout,
        "log.Allrun",
    )
    solver_logs = sorted(case_dir.glob("log.pimpleFoam*"))
    solver_output = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in solver_logs
    )
    redistribute_log = case_dir / "log.redistributePar"
    reconstruct_log = case_dir / "log.reconstructPar"
    topology_case = case_name.endswith("-topologyChange")
    reconstruction_log = redistribute_log if topology_case else reconstruct_log
    reconstruction_completed = (
        reconstruction_log.is_file()
        and "End" in reconstruction_log.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )
    time_steps = len(re.findall(r"(?m)^Time\s*=", solver_output))
    solver_reached_end = "End" in solver_output
    decompose_text = (case_dir / "system/decomposeParDict").read_text(
        encoding="utf-8"
    )
    rank_match = re.search(
        r"(?m)^[ \t]*numberOfSubdomains\s+(\d+)\s*;",
        decompose_text,
    )
    if not rank_match:
        raise ValueError("numberOfSubdomains not found in decomposeParDict")
    return {
        "id": case_name,
        "status": (
            "screen_pass"
            if (
                run_record["return_code"] == 0
                and not run_record["timed_out"]
                and solver_logs
                and time_steps > 0
                and solver_reached_end
                and reconstruction_completed
            )
            else "screen_failure"
        ),
        "end_time": 0.002,
        "write_interval": 0.001,
        "mpi_ranks": int(rank_match.group(1)),
        "input_files": source_files,
        "command": run_record,
        "solver_logs": [str(path) for path in solver_logs],
        "time_steps": time_steps,
        "ami_rebuild_events": solver_output.count("AMI: Creating AMI"),
        "solver_reached_end": solver_reached_end,
        "reconstruction_log": (
            str(reconstruction_log)
            if reconstruction_log.is_file()
            else None
        ),
        "reconstruction_completed": reconstruction_completed,
        "topology_reconstruction_log": (
            str(redistribute_log)
            if topology_case and redistribute_log.is_file()
            else None
        ),
        "topology_reconstruction_completed": (
            reconstruction_completed if topology_case else None
        ),
    }


def target_identity(
    target_id: str,
    project_dir: Path,
    env: dict[str, str],
    external_tool_roots: tuple[Path, ...] = (),
) -> dict:
    expected_project = {
        "B0": B0_PROJECT_DIR,
        "C0": REPO_ROOT,
    }[target_id].resolve(strict=True)
    if project_dir.resolve(strict=True) != expected_project:
        raise ValueError(
            f"{target_id} project directory must be {expected_project}"
        )
    if Path(env.get("WM_PROJECT_DIR", "")).resolve() != expected_project:
        raise ValueError("sourced WM_PROJECT_DIR does not match the target")
    if env.get("WM_PROJECT_VERSION") != "v2606":
        raise ValueError("S03 application screens require v2606")
    if env.get("WM_OPTIONS") != "linux64GccDPInt32Opt":
        raise ValueError("S03 application screens require DP Int32 Opt")

    tool_roots = (expected_project,) + tuple(
        root.resolve() for root in external_tool_roots
    )
    tools = {}
    for tool in (
        "wmake",
        "blockMesh",
        "checkMesh",
        "solidFoam",
        "pimpleFoam",
        "foamDictionary",
        "postProcess",
        "decomposePar",
        "reconstructPar",
        "mpirun",
        "m4",
    ):
        path = shutil.which(tool, path=env.get("PATH"))
        if not path:
            raise ValueError(f"required tool {tool} is unavailable")
        resolved = Path(path).resolve()
        if tool not in {"mpirun", "m4"} and not any(
            resolved.is_relative_to(root) for root in tool_roots
        ):
            raise ValueError(f"{tool} resolves outside {target_id}: {resolved}")
        record = {"path": str(resolved)}
        if tool in {"solidFoam", "pimpleFoam", "blockMesh"}:
            record["sha256"] = runTests.sha256_file(resolved)
        tools[tool] = record

    wmake_version = subprocess.run(
        [tools["wmake"]["path"], "-version"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    commit = None
    dirty = None
    if target_id == "C0":
        commit = subprocess.run(
            ["git", "-C", str(expected_project), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "-C", str(expected_project), "status", "--porcelain"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    return {
        "target_id": target_id,
        "project_dir": str(expected_project),
        "project_version": env["WM_PROJECT_VERSION"],
        "wmake_version": (
            wmake_version.stdout + wmake_version.stderr
        ).strip(),
        "wm_options": env["WM_OPTIONS"],
        "commit": commit,
        "worktree_dirty": dirty,
        "tools": tools,
    }


def load_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-id", required=True, choices=("B0", "C0"))
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=TIMEOUT)
    parser.add_argument(
        "--components",
        nargs="+",
        choices=("solid", "moving-solid", "mixer"),
        default=("solid", "moving-solid", "mixer"),
    )
    return parser.parse_args()


def main() -> int:
    args = load_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError("--timeout must be finite and positive")
    env = os.environ.copy()
    project_dir = args.project_dir.resolve(strict=True)
    identity = target_identity(args.target_id, project_dir, env)

    protected_roots = {
        REPO_ROOT.resolve(),
        B0_PROJECT_DIR.resolve(strict=True),
    }
    run_root = args.run_root.resolve()
    results_dir = args.results_dir.resolve()
    for path in (run_root, results_dir):
        if any(path.is_relative_to(root) for root in protected_roots):
            raise ValueError(f"S03 artifacts cannot be inside an OpenFOAM tree: {path}")
    run_root.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id += "-" + uuid.uuid4().hex[:8]
    run_dir = run_root / f"applications-{run_id}"
    run_dir.mkdir(parents=True, exist_ok=False)
    variants = json.loads(
        (
            TEST_DIR
            / "cases/solidFoam/variants.json"
        ).read_text(encoding="utf-8")
    )

    solid_records = []
    if "solid" in args.components:
        for test_id in ("T14", "T15", "T16", "T17"):
            case_dir = run_dir / "solidFoam" / test_id / "case"
            case_dir.parent.mkdir(parents=True, exist_ok=True)
            solid_records.append(
                run_solid_case(
                    test_id,
                    variants[test_id],
                    case_dir,
                    env,
                    args.timeout,
                    args.target_id,
                )
            )

        conformal_values = variants["T15"].copy()
        conformal_values["rightNy"] = conformal_values["leftNy"]
        conformal_case = run_dir / "solidFoam" / "T15-conformal" / "case"
        conformal_case.parent.mkdir(parents=True, exist_ok=True)
        solid_records.append(
            run_solid_case(
                "T15",
                conformal_values,
                conformal_case,
                env,
                args.timeout,
                args.target_id,
                conformal=True,
            )
        )

        serial_values = variants["T19"].copy()
        serial_case = run_dir / "solidFoam" / "T19" / "serial" / "case"
        serial_case.parent.mkdir(parents=True, exist_ok=True)
        solid_records.append(
            run_solid_case(
                "T19",
                serial_values,
                serial_case,
                env,
                args.timeout,
                args.target_id,
                record_id="T19-serial",
            )
        )

        for ranks, rank_name in ((2, "two-rank"), (4, "four-rank")):
            parallel_values = variants["T19"].copy()
            parallel_values["decompose"] = ranks
            parallel_case = (
                run_dir / "solidFoam" / "T19" / rank_name / "case"
            )
            parallel_case.parent.mkdir(parents=True, exist_ok=True)
            solid_records.append(
                run_solid_case(
                    "T19",
                    parallel_values,
                    parallel_case,
                    env,
                    args.timeout,
                    args.target_id,
                    parallel=True,
                    record_id=f"T19-{rank_name}",
                )
            )

    moving_solid_records = []
    if "moving-solid" in args.components:
        values = variants["T18"].copy()
        values["endTime"] = 0.002
        values["deltaT"] = 0.001
        case_dir = run_dir / "solidFoam" / "T18-moving-AMI" / "case"
        case_dir.parent.mkdir(parents=True, exist_ok=True)
        moving_solid_records.append(
            run_solid_case(
                "T18",
                values,
                case_dir,
                env,
                args.timeout,
                args.target_id,
                dynamic_ami=True,
                record_id="T18-moving-AMI",
            )
        )

        point_motion_values = variants["T18"].copy()
        point_motion_values["endTime"] = 0.002
        point_motion_values["deltaT"] = 0.001
        point_motion_values["motionFunction"] = (
            point_motion_values["cacheMotionFunction"]
        )
        point_motion_values["motionAmplitude"] = (
            point_motion_values["cacheMotionAmplitude"]
        )
        point_motion_values["motionOmega"] = (
            point_motion_values["cacheMotionOmega"]
        )
        control_case = (
            run_dir / "solidFoam" / "T18-moving-point-only-control" / "case"
        )
        control_case.parent.mkdir(parents=True, exist_ok=True)
        moving_solid_records.append(
            run_solid_case(
                "T18",
                point_motion_values.copy(),
                control_case,
                env,
                args.timeout,
                args.target_id,
                dynamic_ami=True,
                dynamic_topology=False,
                record_id="T18-moving-point-only-control",
            )
        )

        cached_case = (
            run_dir / "solidFoam" / "T18-moving-AMI-cache" / "case"
        )
        cached_case.parent.mkdir(parents=True, exist_ok=True)
        moving_solid_records.append(
            run_solid_case(
                "T18",
                point_motion_values,
                cached_case,
                env,
                args.timeout,
                args.target_id,
                dynamic_ami=True,
                dynamic_topology=False,
                cache_ami=True,
                record_id="T18-moving-AMI-cache",
            )
        )

    mixer_records = []
    if "mixer" in args.components:
        for case_name in MIXER_CASES:
            source = project_dir / MIXER_ROOT / case_name
            case_dir = run_dir / "mixer" / case_name
            case_dir.parent.mkdir(parents=True, exist_ok=True)
            mixer_records.append(
                run_mixer_case(
                    case_name,
                    source,
                    case_dir,
                    env,
                    args.timeout,
                )
            )

    serial_mpi_comparisons = []
    if "solid" in args.components:
        serial_t15 = next(
            record for record in solid_records if record["id"] == "T15"
        )
        serial_t19 = next(
            record for record in solid_records
            if record["id"] == "T19-serial"
        )
        parallel_t19_two = next(
            record for record in solid_records
            if record["id"] == "T19-two-rank"
        )
        parallel_t19_four = next(
            record for record in solid_records
            if record["id"] == "T19-four-rank"
        )
        serial_mpi_comparisons.append(
            {
                "id": "T15-serial-vs-T19-two-rank",
                "serial_case": serial_t15["id"],
                "parallel_case": parallel_t19_two["id"],
                **compare_serial_mpi_fields(
                    serial_t15,
                    parallel_t19_two,
                    variants["T15"]["hotT"] - variants["T15"]["coldT"],
                ),
            }
        )
        for parallel_record in (parallel_t19_two, parallel_t19_four):
            serial_mpi_comparisons.append(
                {
                    "id": (
                        f"T19-serial-vs-{parallel_record['id']}"
                    ),
                    "serial_case": serial_t19["id"],
                    "parallel_case": parallel_record["id"],
                    **compare_serial_mpi_fields(
                        serial_t19,
                        parallel_record,
                        (
                            variants["T19"]["hotT"]
                            - variants["T19"]["coldT"]
                        ),
                    ),
                }
            )

    moving_cache_comparison = None
    if "moving-solid" in args.components:
        uncached = next(
            record for record in moving_solid_records
            if record["id"] == "T18-moving-point-only-control"
        )
        cached = next(
            record for record in moving_solid_records
            if record["id"] == "T18-moving-AMI-cache"
        )
        moving_cache_comparison = {
            "id": "T18-moving-cache-vs-no-cache",
            "uncached_case": uncached["id"],
            "cached_case": cached["id"],
            **compare_serial_mpi_fields(
                uncached,
                cached,
                (
                    variants["T18"]["hotT"]
                    - variants["T18"]["coldT"]
                ),
            ),
        }

    records = solid_records + moving_solid_records + mixer_records
    metric_extraction_failures = sum(
        record["status"] == "runtime_pass"
        and (
            not isinstance(record["field_metrics"], dict)
            or "error" in record["field_metrics"]
        )
        for record in solid_records + moving_solid_records
    )
    serial_mpi_comparison_failures = sum(
        comparison["status"] != "compared"
        for comparison in serial_mpi_comparisons
    )
    execution_complete = bool(records) and all(
        not command["timed_out"]
        for record in solid_records + moving_solid_records
        for command in record["commands"].values()
    ) and all(
        not record["command"]["timed_out"]
        for record in mixer_records
    ) and all(
        record["status"] == "runtime_pass"
        for record in solid_records + moving_solid_records
    ) and all(
        record["status"] == "screen_pass"
        for record in mixer_records
    ) and all(
        record["gamg_solver_events"] > 0
        for record in solid_records
        if record["id"].startswith("T19-")
    ) and (
        moving_cache_comparison is None
        or moving_cache_comparison["status"] == "compared"
    ) and metric_extraction_failures == 0 and serial_mpi_comparison_failures == 0

    report_path = (
        results_dir / f"S03-applications-{args.target_id}-{run_id}.json"
    )
    manifest_path = (
        results_dir / f"S03-applications-inputs-{args.target_id}-{run_id}.json"
    )
    if report_path.exists() or manifest_path.exists():
        raise FileExistsError("refusing to overwrite an S03 application report")

    suite_inputs = []
    for path in sorted(TEST_DIR.rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts:
            suite_inputs.append(
                {
                    "path": str(path.relative_to(REPO_ROOT)),
                    "size_bytes": path.stat().st_size,
                    "sha256": runTests.sha256_file(path),
                }
            )
    manifest = {
        "schema_version": 1,
        "stage": "S03",
        "run_id": run_id,
        "target": identity,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command_line": [sys.executable, *sys.argv],
        "timeout_seconds": args.timeout,
        "components": args.components,
        "suite_inputs": suite_inputs,
        "report": str(report_path),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    report = {
        "schema_version": 1,
        "stage": "S03",
        "scope": "solidFoam_runtime_and_mixer_tutorial_screens",
        "run_id": run_id,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "report_path": str(report_path),
        "input_manifest": str(manifest_path),
        "target": identity,
        "solid_acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
        "solid_cases": solid_records,
        "moving_solid_cases": moving_solid_records,
        "solid_serial_mpi_comparisons": serial_mpi_comparisons,
        "moving_solid_cache_comparison": moving_cache_comparison,
        "mixer_tutorials": mixer_records,
        "metric_extraction_failures": metric_extraction_failures,
        "serial_mpi_comparison_failures": serial_mpi_comparison_failures,
        "application_observation_status": (
            "complete" if execution_complete else "incomplete"
        ),
        "solver_failures": sum(
            record["status"] != "runtime_pass"
            for record in solid_records + moving_solid_records
        ),
        "tutorial_failures": sum(
            record["status"] != "screen_pass"
            for record in mixer_records
        ),
    }
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"Report: {report_path}")
    print(
        "S03 application observations: "
        f"{report['application_observation_status']}"
    )
    return 0 if execution_complete else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"S03 application runner error: {error}", file=sys.stderr)
        raise SystemExit(2)
