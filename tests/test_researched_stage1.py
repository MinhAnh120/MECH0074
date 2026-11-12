"""Check correlations, independent timing, and the load-export contract."""

from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import params as p
from complete_stage1 import comparison, shared_case
from export_bends import export_route
from fluid_assumptions import colebrook, friction_estimate, methane_saturation_pa
from run_stage1 import metrics
from stage1 import simulate, theory


class ResearchedTests(unittest.TestCase):
    def test_methane_fixed_points_and_111K(self):
        # Independent physical fixed points published by Setzmann & Wagner.
        self.assertAlmostEqual(methane_saturation_pa(111.668)/101325, 1, delta=1e-4)
        self.assertAlmostEqual(methane_saturation_pa(90.6941)/11696, 1, delta=1e-4)
        self.assertAlmostEqual(methane_saturation_pa(111), 95873.898949157, places=5)
        for t in (80, 200, float("nan")):
            with self.assertRaises(ValueError):
                methane_saturation_pa(t)

    def test_colebrook_independent_residual_and_monotonicity(self):
        for re, rough in [(1e5, 0), (4.5e6, 15e-6/.79394), (4e6, 50e-6/.79394)]:
            f = colebrook(re, rough)
            residual = f**-.5+2*np.log10(rough/3.7+2.51/(re*np.sqrt(f)))
            self.assertLess(abs(residual), 1e-12)
        self.assertAlmostEqual(colebrook(1e5, 0), .017989773084, places=10)
        bounds = friction_estimate()
        self.assertLess(bounds["lower"], bounds["darcy_factor"])
        self.assertLess(bounds["darcy_factor"], bounds["upper"])
        with self.assertRaises(ValueError):
            colebrook(100, 0)

    def test_shared_configuration_and_phase_crossing(self):
        for mode in ("bulk", "elastic"):
            case = shared_case(mode, friction_estimate()["darcy_factor"])
            self.assertEqual(case.diameter, p.D_IN)
            self.assertEqual(case.pressure_abs, p.P0_ABS)
            result = simulate(case, 100, 1)
            self.assertIsNotNone(result["first_phase_crossing_s"])
            self.assertLess(result["minimum_absolute_pa"], case.phase_limit)

    def test_boundary_travel_measurements(self):
        for mode in ("bulk", "elastic"):
            case = shared_case(mode, 0)
            result = simulate(case, 100, 1)
            rows = comparison(case, result, metrics(case, result), theory(case)["fundamental_hz"])
            for row in rows:
                if row[0].startswith(("reservoir_", "first_arrival")):
                    self.assertIsNotNone(row[3])
                    self.assertLessEqual(abs(row[4]), result["dt"]+1e-10)

    def test_export_contract_with_synthetic_geometry_only(self):
        config = {"length_m": 100, "bend_stations_m": [11, 29, 53, 87],
                  "speed_mode": "elastic", "geometry_source": "SYNTHETIC TEST ONLY",
                  "reaches": 100, "duration_s": .3}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            meta = export_route(config, root, "test-version")
            data = np.loadtxt(root/"bend_histories.csv", delimiter=",", skiprows=1)
            self.assertEqual(data.shape[1], 10)
            np.testing.assert_allclose(np.diff(data[:, 0]), meta["dt_s"], atol=1e-12)
            case = replace(shared_case("elastic", friction_estimate()["darcy_factor"]), length=100)
            direct = simulate(case, 100, .3, stations=config["bend_stations_m"])
            np.testing.assert_allclose(data[:, 1:5], direct["pressure"], atol=1e-5, rtol=0)
            np.testing.assert_allclose(data[:, 5:9], direct["flow"], atol=1e-11, rtol=0)
            self.assertFalse(meta["approved_physical_load"])
            self.assertIn(0, data[:, -1])
            self.assertEqual(json.loads((root/"metadata.json").read_text())["code_version"], "test-version")
            with self.assertRaises(ValueError):
                export_route({**config, "bend_stations_m": None}, root)


if __name__ == "__main__":
    unittest.main()
