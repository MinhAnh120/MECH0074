"""Independent physical checks for the adapted teaching model."""

from dataclasses import replace
import unittest
import numpy as np
from stage1 import Case, G, estimate_period, simulate, spectrum, theory
from run_stage1 import metrics, refinement_change


def reflected_ramps(case, time, stations):
    """Independent frictionless solution from travelling/reflected finite ramps."""
    t, x = time[:, None], stations[None, :]
    pressure = np.full((time.size, stations.size), case.pressure_abs)
    amplitude = case.density*case.speed*case.velocity
    for m in range(int(np.ceil(time[-1]*case.speed/(2*case.length)))+1):
        outward = ((2*m+1)*case.length-x)/case.speed
        returned = ((2*m+1)*case.length+x)/case.speed
        pressure += amplitude*(-1)**m*(np.clip((t-outward)/case.closure_time, 0, 1)
                                      - np.clip((t-returned)/case.closure_time, 0, 1))
    return pressure


class Stage1Tests(unittest.TestCase):
    def test_analytical_values(self):
        case = Case()
        self.assertAlmostEqual(case.speed, 997.899867315862, places=8)
        self.assertAlmostEqual(case.area*case.velocity, 0.742602669945, places=10)
        self.assertAlmostEqual(theory(case)["surge_pa"], 681066.659443, places=5)
        self.assertAlmostEqual(theory(replace(case, speed_mode="bulk"))["fundamental_hz"], 2.75)

    def test_parameter_validation(self):
        for change in ({"friction": -1}, {"diameter": 0.9}, {"density": float("nan")},
                       {"closure_time": 0}, {"phase_limit": -1}, {"speed_mode": "unknown"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                Case(**change)
        for n in (0, 2.5, True):
            with self.assertRaises(ValueError):
                simulate(Case(), n, 1)
        with self.assertRaises(ValueError):
            simulate(Case(), 100, 1, stations=[101])

    def test_frictional_steady_state(self):
        case = Case(friction=0.020)
        result = simulate(case, 100, 0.5, stations=[0, 10, 50, 100], hold_open=True)
        expected = case.pressure_abs-case.friction*result["stations"]*case.density*case.velocity**2/(2*case.diameter)
        np.testing.assert_allclose(result["pressure"], np.tile(expected, (len(result["time"]), 1)), atol=1e-6, rtol=0)
        np.testing.assert_allclose(result["flow"], case.area*case.velocity, atol=1e-12, rtol=0)

    def test_boundaries(self):
        case = Case(friction=0.020)
        result = simulate(case, 100, 1, stations=[0, 100])
        np.testing.assert_allclose(result["pressure"][:, 0], case.pressure_abs, atol=1e-7, rtol=0)
        expected = case.area*case.velocity*np.maximum(1-result["time"]/case.closure_time, 0)
        np.testing.assert_allclose(result["flow"][:, -1], expected, atol=1e-14, rtol=0)

    def test_worked_frictional_update(self):
        case = Case(friction=0.020)
        dt = 1/case.speed
        result = simulate(case, 100, 2*dt, stations=[99, 100])
        head = (result["pressure"]-case.external_abs)/(case.density*G)
        self.assertAlmostEqual(head[1, -1], 141.481678233121, places=9)
        self.assertAlmostEqual(head[2, 0], 141.484425968198, places=9)
        self.assertAlmostEqual(result["flow"][2, 0], 0.705395080943, places=11)
        self.assertAlmostEqual(result["pressure"][1, -1]-result["pressure"][0, -1], 34125, places=6)

    def test_full_finite_closure_histories(self):
        for mode in ("elastic", "bulk"):
            case = Case(speed_mode=mode)
            result = simulate(case, 100, 1)
            expected = reflected_ramps(case, result["time"], result["stations"])
            np.testing.assert_allclose(result["pressure"], expected, atol=1e-5, rtol=0)
            period = estimate_period(result["time"], result["pressure"][:, -1]-case.pressure_abs, case.closure_time)
            self.assertAlmostEqual(period, theory(case)["period_s"], places=10)

    def test_arrivals_and_bulk_ten_metre_overlap(self):
        case = Case(speed_mode="bulk")
        result = simulate(case, 100, 0.15)
        delta = result["pressure"]-result["pressure"][0]
        for i, x in enumerate(result["stations"]):
            index = np.flatnonzero(abs(delta[:, i]) > 1)[0]
            arrival = (case.length-x)/case.speed
            self.assertLessEqual(abs(result["time"][index]-arrival), result["dt"]+1e-12)
        self.assertAlmostEqual(delta[:, 0].max(), 682500, places=5)

    def test_station_interpolation_is_explicit(self):
        case = Case(friction=0.020)
        result = simulate(case, 100, 0.01, stations=[10.25, 100], hold_open=True)
        self.assertTrue(result["interpolated_stations"])
        expected = case.pressure_abs-case.friction*10.25*case.density*case.velocity**2/(2*case.diameter)
        self.assertAlmostEqual(result["pressure"][0, 0], expected, places=7)

    def test_fft_frequency_and_amplitude(self):
        case = Case()
        time = np.arange(20000)/case.speed
        f = theory(case)["fundamental_hz"]
        signal = 2e6+3e5*np.sin(2*np.pi*f*time)+5e4*np.sin(6*np.pi*f*time)
        frequency, amplitude, peak = spectrum(time, signal)
        self.assertAlmostEqual(peak, f, places=10)
        self.assertAlmostEqual(amplitude[np.argmax(amplitude[1:])+1]/3e5, 1, places=4)
        with self.assertRaises(ValueError):
            spectrum([0, 1, 2, 4], [0, 1, 0, -1])

    def test_phase_diagnostics_do_not_claim_physical_validation(self):
        unknown = simulate(Case(), 100, 0.3)
        self.assertEqual(unknown["phase_status"], "LNG_phase_limit_not_supplied")
        entered = simulate(Case(phase_limit=50000), 100, 0.3)
        self.assertEqual(entered["phase_status"], "phase_limit_crossed")
        negative = simulate(Case(speed_mode="bulk"), 100, 0.3)
        self.assertEqual(negative["phase_status"], "nonpositive_absolute_pressure")

    def test_frictional_grid_refinement(self):
        case = Case(friction=0.020)
        previous = previous_metrics = None
        for n in (100, 200, 400):
            result = simulate(case, n, 1)
            row = metrics(case, result)
            self.assertIsNotNone(row["period_s"])
            if previous is not None:
                change = refinement_change(previous, result, previous_metrics, row,
                                           theory(case)["surge_pa"])
                for value in change.values():
                    self.assertLess(value, 1.0)  # Proposed numerical tolerance, not a UCL threshold.
            previous, previous_metrics = result, row


if __name__ == "__main__":
    unittest.main()
