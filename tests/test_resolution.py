"""Independent checks for the expanded grid study and its error definition."""

import unittest
import numpy as np
from complete_stage1 import shared_case
from fluid_assumptions import friction_estimate
from resolution import (DEFAULT_REACHES, meets_tolerance, refinement_change,
                        resolution_study, validate_reaches)
from run_stage1 import metrics


class ResolutionTests(unittest.TestCase):
    def test_between_sample_error_is_not_missed(self):
        coarse = {"time": np.array([0., .5, 1.]), "stations": np.array([25.]),
                  "pressure": np.zeros((3, 1))}
        fine = {"time": np.array([0., .25, .5, .75, 1.]), "stations": np.array([25.]),
                "pressure": np.array([[0.], [1.], [0.], [1.], [0.]])}
        measure = {"first_valve_rise_pa": 1., "period_s": 1.}
        change = refinement_change(coarse, fine, measure, measure, 1.)
        self.assertEqual(change["waveform_difference_percent"], 100.)
        self.assertFalse(meets_tolerance(change))

    def test_six_grids_and_explicit_station_sampling(self):
        case = shared_case("elastic", friction_estimate()["darcy_factor"])
        finest, rows, passing = resolution_study(case, DEFAULT_REACHES, metrics)
        self.assertEqual([r["reaches"] for r in rows], list(DEFAULT_REACHES))
        self.assertAlmostEqual(finest["dx"], 100/320)
        self.assertTrue(rows[0]["interpolated_stations"])
        self.assertTrue(all(not r["interpolated_stations"] for r in rows[1:]))
        self.assertIsNone(rows[-1]["difference_from_finest"])
        self.assertIsNone(rows[-1]["meets_criterion_vs_finest"])
        self.assertFalse(rows[0]["meets_criterion_vs_finest"])
        self.assertTrue(rows[-2]["meets_criterion_vs_next"])
        self.assertIn(passing, DEFAULT_REACHES[:-1])

    def test_strict_one_percent_and_missing_measurement(self):
        self.assertTrue(meets_tolerance({"peak": .5, "period": .2, "waveform": .99}))
        self.assertFalse(meets_tolerance({"peak": 1., "period": .2, "waveform": .5}))
        self.assertFalse(meets_tolerance({"peak": .5, "period": None, "waveform": .5}))
        self.assertFalse(meets_tolerance(None))

    def test_refinement_input_validation(self):
        for reaches in ([10, 20], [10, 30, 60], [0, 0, 0], [10., 20., 40.],
                        [True, 2, 4], [40, 20, 10]):
            with self.subTest(reaches=reaches), self.assertRaises(ValueError):
                validate_reaches(reaches)
        self.assertEqual(validate_reaches(DEFAULT_REACHES), DEFAULT_REACHES)


if __name__ == "__main__":
    unittest.main()
