#!/usr/bin/env python3

import json
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
from solidMetrics import (
    normalized_error,
    rotated_diagonal_tensor,
    series_resistance_heat_rate,
    series_resistance_interface_temperature,
    solid_reference_metrics,
    tensor_flux,
    volumetric_heat_capacity,
)


class SolidMetricsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        variants_path = Path(__file__).parent / "cases" / "solidFoam" / "variants.json"
        cls.variants = json.loads(variants_path.read_text(encoding="utf-8"))

    def test_t15_series_resistance_reference(self):
        metrics = solid_reference_metrics("T15", self.variants["T15"])
        self.assertAlmostEqual(metrics["analytical_heat_rate"], 90.9090909091)
        self.assertAlmostEqual(metrics["interface_temperature"], 282.2409090909)

    def test_t16_volumetric_heat_capacity(self):
        metrics = solid_reference_metrics("T16", self.variants["T16"])
        self.assertEqual(
            metrics["volumetric_heat_capacity"],
            {"left": 500000.0, "right": 800000.0},
        )

    def test_t17_rotated_tensor_has_off_diagonal_flux(self):
        metrics = solid_reference_metrics("T17", self.variants["T17"])
        left_flux = metrics["reference_fluxes"]["left"]
        self.assertAlmostEqual(left_flux[0], -2.5)
        self.assertAlmostEqual(left_flux[1], -1.5)
        self.assertNotEqual(left_flux[1], 0)

    def test_zero_temperature_difference_gives_zero_heat_rate(self):
        metrics = solid_reference_metrics("T14", self.variants["T14"])
        self.assertEqual(metrics["expected_equilibrium_heat_rate"], 0)

    def test_normalized_error_uses_declared_scale(self):
        self.assertAlmostEqual(normalized_error(101, 100, 100), 0.01)

    def test_rejects_invalid_reference_inputs(self):
        with self.assertRaises(ValueError):
            series_resistance_heat_rate(0, (1,), (1,), 2, 1)
        with self.assertRaises(ValueError):
            volumetric_heat_capacity(1, float("nan"))
        with self.assertRaises(ValueError):
            rotated_diagonal_tensor((1, 2, 3), (1, 0, 0), (1, 0, 0))
        with self.assertRaises(ValueError):
            series_resistance_interface_temperature(300, 10, 0, 1, 1)
        with self.assertRaises(ValueError):
            tensor_flux(((1, 0), (0, 1)), (1, 0, 0))


if __name__ == "__main__":
    unittest.main()
