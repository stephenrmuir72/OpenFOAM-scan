"""Small independent reference calculations for the solid AMI fixtures."""

from __future__ import annotations

from math import isfinite, sqrt
from typing import Sequence


SOLID_ACCEPTANCE_LIMITS = {
    "equilibrium_temperature_change": 1e-10,
    "equilibrium_spurious_heat_rate": 1e-10,
    "steady_profile_error": 1e-3,
    "steady_heat_rate_error": 1e-3,
    "stationary_energy_balance": 1e-8,
    "interface_heat_rate_balance": 1e-8,
    "serial_mpi_difference": 1e-8,
}


def _finite(value: float, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _vector3(values: Sequence[float], name: str) -> tuple[float, float, float]:
    if len(values) != 3:
        raise ValueError(f"{name} must contain three components")
    return tuple(_finite(value, name) for value in values)


def _dot(
    left: Sequence[float],
    right: Sequence[float],
) -> float:
    return sum(a*b for a, b in zip(left, right))


def _cross(
    left: Sequence[float],
    right: Sequence[float],
) -> tuple[float, float, float]:
    return (
        left[1]*right[2] - left[2]*right[1],
        left[2]*right[0] - left[0]*right[2],
        left[0]*right[1] - left[1]*right[0],
    )


def _norm(values: Sequence[float]) -> float:
    return sqrt(_dot(values, values))


def series_resistance_heat_rate(
    area: float,
    lengths: Sequence[float],
    conductivities: Sequence[float],
    hot_temperature: float,
    cold_temperature: float,
) -> float:
    """Return Q for layered 1-D conduction with no contact resistance."""
    area_value = _positive(area, "area")
    if not lengths or len(lengths) != len(conductivities):
        raise ValueError("lengths and conductivities must have equal nonzero size")

    resistance = 0.0
    for index, (length, conductivity) in enumerate(zip(lengths, conductivities)):
        resistance += _positive(length, f"lengths[{index}]") / (
            _positive(conductivity, f"conductivities[{index}]")*area_value
        )

    return (
        _finite(hot_temperature, "hot_temperature")
        - _finite(cold_temperature, "cold_temperature")
    )/resistance


def series_resistance_interface_temperature(
    hot_temperature: float,
    heat_rate: float,
    area: float,
    first_length: float,
    first_conductivity: float,
) -> float:
    """Return the material-junction temperature for a two-layer slab."""
    return (
        _finite(hot_temperature, "hot_temperature")
        - _finite(heat_rate, "heat_rate")
        * _positive(first_length, "first_length")
        / (
            _positive(first_conductivity, "first_conductivity")
            * _positive(area, "area")
        )
    )


def volumetric_heat_capacity(density: float, heat_capacity: float) -> float:
    """Return rho*Cp after checking that both material values are positive."""
    return _positive(density, "density")*_positive(heat_capacity, "heat_capacity")


def rotated_diagonal_tensor(
    eigenvalues: Sequence[float],
    e1: Sequence[float],
    e2: Sequence[float],
    tolerance: float = 1e-12,
) -> tuple[tuple[float, float, float], ...]:
    """Build R*diag(k)*R^T from two orthonormal material axes."""
    if len(eigenvalues) != 3:
        raise ValueError("eigenvalues must contain three components")

    values = tuple(
        _positive(value, f"eigenvalues[{index}]")
        for index, value in enumerate(eigenvalues)
    )
    axis1 = _vector3(e1, "e1")
    axis2 = _vector3(e2, "e2")
    tol = _positive(tolerance, "tolerance")

    if (
        abs(_norm(axis1) - 1) > tol
        or abs(_norm(axis2) - 1) > tol
        or abs(_dot(axis1, axis2)) > tol
    ):
        raise ValueError("e1 and e2 must be orthonormal")

    axis3 = _cross(axis1, axis2)
    basis = (axis1, axis2, axis3)
    return tuple(
        tuple(
            sum(values[axis]*basis[axis][row]*basis[axis][column]
                for axis in range(3))
            for column in range(3)
        )
        for row in range(3)
    )


def tensor_flux(
    conductivity: Sequence[Sequence[float]],
    temperature_gradient: Sequence[float],
) -> tuple[float, float, float]:
    """Return the full-tensor conductive flux, q = -K grad(T)."""
    if len(conductivity) != 3 or any(len(row) != 3 for row in conductivity):
        raise ValueError("conductivity must be a 3 by 3 matrix")
    gradient = _vector3(temperature_gradient, "temperature_gradient")
    matrix = tuple(
        tuple(
            _finite(value, "conductivity")
            for value in row
        )
        for row in conductivity
    )
    return tuple(
        -sum(matrix[row][column]*gradient[column] for column in range(3))
        for row in range(3)
    )


def normalized_error(
    observed: float,
    expected: float,
    reference_scale: float,
) -> float:
    """Return an absolute error normalized by a declared positive scale."""
    scale = _positive(reference_scale, "reference_scale")
    return abs(
        _finite(observed, "observed") - _finite(expected, "expected")
    )/scale


def solid_reference_metrics(test_id: str, values: dict) -> dict:
    """Return analytical reference data for one parameterized solid fixture."""
    area = values["referenceArea"]
    lengths = (values["leftLength"], values["rightLength"])
    conductivities = (values["leftKappa"][0], values["rightKappa"][0])
    delta_temperature = values["hotT"] - values["coldT"]

    if test_id == "T14" or test_id == "T18":
        return {
            "expected_equilibrium_heat_rate": series_resistance_heat_rate(
                area,
                lengths,
                conductivities,
                values["hotT"],
                values["coldT"],
            ),
            "acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
        }

    if test_id == "T15" or test_id == "T19":
        heat_rate = series_resistance_heat_rate(
            area,
            lengths,
            conductivities,
            values["hotT"],
            values["coldT"],
        )
        return {
            "analytical_heat_rate": heat_rate,
            "interface_temperature": series_resistance_interface_temperature(
                values["hotT"],
                heat_rate,
                area,
                lengths[0],
                conductivities[0],
            ),
            "temperature_scale": abs(delta_temperature),
            "heat_rate_scale": abs(heat_rate),
            "acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
        }

    if test_id == "T16":
        return {
            "volumetric_heat_capacity": {
                "left": volumetric_heat_capacity
                (
                    values["leftRho"],
                    values["leftCp"],
                ),
                "right": volumetric_heat_capacity
                (
                    values["rightRho"],
                    values["rightCp"],
                ),
            },
            "acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
        }

    if test_id == "T17":
        gradient = values["referenceGradient"]
        tensors = {
            zone: rotated_diagonal_tensor(
                values[f"{zone}Kappa"],
                values["e1"],
                values["e2"],
            )
            for zone in ("left", "right")
        }
        fluxes = {
            zone: tensor_flux(tensors[zone], gradient)
            for zone in ("left", "right")
        }
        return {
            "conductivity_tensors": tensors,
            "reference_gradient": gradient,
            "reference_fluxes": fluxes,
            "acceptance_limits": SOLID_ACCEPTANCE_LIMITS,
        }

    raise ValueError(f"no solid reference metrics for {test_id}")
