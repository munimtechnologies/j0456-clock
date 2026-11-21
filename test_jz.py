# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests (standard library only; no network, no astropy)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from jz import changepoint, ephem, geometry, grade, photometry, stats  # noqa: E402

CDELT = 6.548089e-4  # Swift-XRT sky pixel, degrees


class TestStats(unittest.TestCase):
    def test_gamma_a1_both_branches(self):
        for x in (0.3, 2.5, 9.0):
            self.assertAlmostEqual(stats.gammainc_lower(1.0, x), 1.0 - math.exp(-x), places=12)
            self.assertAlmostEqual(stats.gammainc_upper(1.0, x), math.exp(-x), places=12)

    def test_poisson_tail(self):
        want = 1.0 - math.exp(-1.0) * (1.0 + 1.0 + 0.5)
        self.assertAlmostEqual(stats.poisson_tail(3, 1.0), want, places=12)
        self.assertEqual(stats.poisson_tail(0, 1.0), 1.0)

    def test_sigma(self):
        p3 = 0.5 * math.erfc(3.0 / math.sqrt(2.0))
        self.assertAlmostEqual(stats.sigma_from_p(p3), 3.0, places=6)
        self.assertEqual(stats.sigma_from_p(0.7), 0.0)

    def test_upper_limits(self):
        ul0 = -math.log(1.0 - 0.9973)
        self.assertAlmostEqual(stats.source_upper_limit(0, 0.0), ul0, places=6)
        self.assertAlmostEqual(stats.source_upper_limit(0, 2.0), ul0, places=6)
        s = stats.source_upper_limit(1, 0.0, cl=0.9)
        self.assertAlmostEqual(s, 3.8897, places=3)


class TestGeometry(unittest.TestCase):
    def test_reference_pixel_and_orientation(self):
        x, y = geometry.tan_world_to_pixel(74.2075, -20.63, 74.2075, -20.63, 500.5, 500.5,
                                           -CDELT, CDELT)
        self.assertAlmostEqual(x, 500.5, places=9)
        self.assertAlmostEqual(y, 500.5, places=9)
        x, y = geometry.tan_world_to_pixel(74.2075, -20.62, 74.2075, -20.63, 500.5, 500.5,
                                           -CDELT, CDELT)
        self.assertAlmostEqual(y - 500.5, 0.01 / CDELT, delta=0.01)

    def test_eef(self):
        self.assertTrue(0.85 < photometry.EEF < 0.88)

    def test_archive_month(self):
        self.assertEqual(geometry.mjd_to_archive_month(60202.15), "2023_09")
        self.assertEqual(geometry.mjd_to_archive_month(60202.15, 1), "2023_10")


class TestPhotometry(unittest.TestCase):
    def test_area_ratio(self):
        self.assertAlmostEqual(photometry.AREA_RATIO, 1225.0 / 35100.0, places=12)

    def test_count_aperture(self):
        ps = CDELT * 3600.0
        xs = [500.0] * 7 + [500.0 + 100.0 / ps, 500.0 + 50.0 / ps, 500.0 + 300.0 / ps]
        ys = [500.0] * 10
        pis = [100] * 5 + [10, 250] + [100, 100, 100]
        self.assertEqual(photometry.count_aperture(xs, ys, pis, 500.0, 500.0, ps), (5, 1))

    def test_measure(self):
        m = photometry.measure(0, 10, 1000.0)
        self.assertFalse(m["detected"])
        self.assertEqual(m["level"], m["ul3_rate"])
        m = photometry.measure(50, 10, 1000.0)
        self.assertTrue(m["detected"])
        want = (50 - 10 * photometry.AREA_RATIO) / (1000.0 * photometry.EEF)
        self.assertAlmostEqual(m["rate"], want, places=12)
        self.assertFalse(photometry.measure(3, 0, 1000.0)["detected"])


class TestEphem(unittest.TestCase):
    def test_numbers(self):
        self.assertEqual(ephem.RECUR, 190.5)
        self.assertEqual(ephem.DROP_OFFSET, 107.0)
        self.assertEqual(ephem.rise_start(6), 60162.5)
        self.assertEqual(ephem.predicted_drop(6), 60265.5)
        self.assertEqual(ephem.on_window(6), (60197.5, 60257.5))
        self.assertEqual(ephem.off_window(5), (60092.0, 60157.0))
        self.assertEqual(ephem.off_window(6), (60282.5, 60347.5))

    def test_phase_class(self):
        self.assertEqual(ephem.phase_class(60230.0), "ON")
        self.assertEqual(ephem.phase_class(60300.0), "OFF")
        self.assertEqual(ephem.phase_class(60265.5), "DROP_ZONE")
        self.assertEqual(ephem.phase_class(59970.0), "PRE_RISE")
        self.assertEqual(ephem.cycle_of(59970.0), 4)
        self.assertEqual(ephem.phase_class(60963.0), "ON")
        self.assertEqual(ephem.cycle_of(60963.0), 10)
        self.assertEqual(ephem.phase_class(60170.0), "RISE")


def obs(mjd, n_src, n_bkg, exp):
    return {"obsid": "%011d" % int(mjd * 10), "mjd": mjd, "n_src": n_src, "n_bkg": n_bkg,
            "exposure_s": exp}


class TestChangepoint(unittest.TestCase):
    def test_fit_rate_single(self):
        r, _ = changepoint.fit_rate([(0.0, 10, 2.0, 100.0)])
        self.assertAlmostEqual(r, 0.08, places=8)
        r, _ = changepoint.fit_rate([(0.0, 0, 2.0, 100.0)])
        self.assertEqual(r, 0.0)

    def test_step(self):
        rows = [obs(10.0 * i, 20 if i < 5 else 0, 0, 1000.0) for i in range(10)]
        cp = changepoint.changepoint(rows, 0.1, 1.0)
        self.assertEqual(cp["bracket"], [40.0, 50.0])
        self.assertAlmostEqual(cp["r_before"], 0.02, places=4)
        self.assertEqual(cp["r_after"], 0.0)
        self.assertGreater(cp["ts"], 50.0)
        self.assertEqual(cp["n_obs"], 10)

    def test_too_few(self):
        self.assertIsNone(changepoint.changepoint([obs(1.0, 1, 1, 100.0)] * 3, 0.1, 1.0))


def history():
    rows = [obs(m, 50, 30, 3000.0) for m in (59650.0, 59670.0, 59690.0)]
    rows += [obs(m, 30, 10, 3000.0) for m in (59820.0, 59830.0, 59840.0, 59850.0, 59860.0,
                                              59870.0, 59876.0)]
    rows += [obs(m, 0, 10, 3000.0) for m in (59886.0, 59900.0, 59920.0, 59940.0, 59960.0)]
    rows += [obs(m, 15, 10, 3000.0) for m in (60010.0, 60030.0, 60060.0)]
    return rows


OFF5_T = (60120.0, 60130.0, 60140.0, 60150.0)
ON6_T = (60205.0, 60215.0, 60225.0, 60235.0, 60245.0, 60255.0)
REST6_T = (60270.0, 60290.0, 60300.0, 60310.0, 60320.0, 60330.0, 60340.0)


class TestGrade(unittest.TestCase):
    def test_gates_pass(self):
        g = grade.run_gates(history())
        self.assertTrue(g["all_passed"], g["status"])
        self.assertEqual(g["cycle4_changepoint"]["bracket"], [59876.0, 59886.0])

    def test_gate_fails_without_drop(self):
        rows = [r for r in history() if not 59886.0 <= r["mjd"] <= 59960.0]
        rows += [obs(m, 30, 10, 3000.0) for m in (59886.0, 59900.0, 59920.0, 59940.0, 59960.0)]
        g = grade.run_gates(rows)
        self.assertEqual(g["status"]["G3_faint_phase_below_half_cycle4"], "FAIL")
        self.assertFalse(g["all_passed"])

    def test_on_schedule(self):
        rows = history() + [obs(m, 0, 10, 1000.0) for m in OFF5_T]
        rows += [obs(m, 10, 10, 1000.0) for m in ON6_T]
        rows += [obs(m, 0, 10, 1000.0) for m in REST6_T]
        r = grade.grade(rows)
        self.assertEqual(r["outcome"], "FLARE_ON_SCHEDULE")
        self.assertGreater(r["z_on_vs_off"], 5.0)
        self.assertEqual(r["T2_amplitude"]["status"], "FADING_SLOWER")
        # sqrt((45 - 30a) / (150 - 90a)) with a = AREA_RATIO = 0.0349
        self.assertAlmostEqual(r["fade_forecast"]["per_cycle_factor"], 0.5471, delta=0.002)
        self.assertEqual(r["T3_timing"]["status"], "ON_TIME")
        self.assertEqual(r["T3_timing"]["changepoint"]["bracket"], [60255.0, 60270.0])
        self.assertAlmostEqual(r["T3_timing"]["drop_to_drop_recurrence_days"], 191.0, places=9)

    def test_no_flare(self):
        rows = history() + [obs(m, 0, 10, 1000.0) for m in OFF5_T + ON6_T + REST6_T]
        r = grade.grade(rows)
        self.assertEqual(r["outcome"], "NO_FLARE")
        self.assertEqual(r["depth_reading"], "BELOW_FORECAST_LEVEL")
        self.assertEqual(r["T2_amplitude"]["status"], "NOT_APPLICABLE")
        self.assertEqual(r["T3_timing"]["status"], "UNTIMED")

    def test_off_schedule(self):
        rows = history() + [obs(m, 0, 10, 1000.0) for m in OFF5_T + ON6_T]
        rows += [obs(m, 15 if 60290.0 <= m <= 60320.0 else 0, 10, 1000.0) for m in REST6_T]
        r = grade.grade(rows)
        self.assertEqual(r["outcome"], "FLARE_OFF_SCHEDULE")
        self.assertEqual(r["free_phase_scan"]["mjd_first"], 60290.0)
        self.assertEqual(r["free_phase_scan"]["mjd_last"], 60320.0)

    def test_no_data(self):
        self.assertEqual(grade.grade(history())["outcome"], "NO_DATA")

    def test_late_groups(self):
        rows = [obs(m, 0, 5, 1000.0) for m in (60502.0, 60512.0, 60563.0, 60850.0, 60851.0,
                                               61282.0)]
        g = grade.late_groups(rows)
        self.assertEqual(len(g), 4)
        self.assertEqual(g[0]["n_obs"], 2)
        self.assertEqual(g[0]["forecast_phases"], ["OFF", "OFF"])


if __name__ == "__main__":
    unittest.main()
