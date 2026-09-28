#!/usr/bin/env python3
"""Synthetic checks for the intermittent nearest-neighbour correlation."""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "compute_nearest_neighbour.py"
SPEC = importlib.util.spec_from_file_location("nn", SCRIPT)
nn = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = nn
SPEC.loader.exec_module(nn)


class NearestNeighbourCorrelationTests(unittest.TestCase):
    def test_reformation_is_counted(self):
        corr = nn.IntermittentCorrelator(max_lag=2)
        corr.update({1}, 0)
        corr.update(set(), 1)
        corr.update({1}, 2)
        result = corr.finish()
        self.assertEqual(result[0], 1.0)
        self.assertEqual(result[1], 0.0)
        self.assertEqual(result[2], 1.0)

    def test_first_crossing_is_interpolated(self):
        time = nn.np.array([0.0, 1.0, 2.0])
        corr = nn.np.array([1.0, 0.5, 0.25])
        crossing = nn.first_crossing(time, corr)
        expected = 1.0 + (nn.np.exp(-1.0) - 0.5) / (0.25 - 0.5)
        self.assertAlmostEqual(crossing, expected)


if __name__ == "__main__":
    unittest.main()
