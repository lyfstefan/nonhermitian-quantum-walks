"""Protect plotted probabilities and the common geometry of the approved bars."""

import unittest
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from qwalk.figures import build_case_figure, draw_profile, profile_labels


def fixture():
    sites = np.arange(-15, 16)
    profiles = np.column_stack(
        [
            np.exp(-abs(sites - center) / scale)
            for center, scale in [(0, 1), (0, 2), (-2, 1.5), (2, 0.8)]
        ]
    )
    profiles /= profiles.sum(axis=0)
    large = np.arange(-35, 36)
    probability = np.exp(-abs(large[None, :] - 0.02 * np.arange(301)[:, None]))
    probability /= probability.sum(axis=1)[:, None]
    return {
        "eigen_sites": sites,
        "eigen_profiles": profiles,
        "moduli": np.array([1.41, 0.71, 1.47, 0.68]),
        "sites": large,
        "probability": probability,
    }


class Figures(unittest.TestCase):
    def test_actual_bar_heights_and_evolution_values_are_unchanged(self):
        data = fixture()
        before = {k: v.copy() for k, v in data.items()}
        heights = []
        original = Axes3D.bar3d

        def observe(ax, x, y, z, dx, dy, dz, *args, **kwargs):
            heights.append(np.asarray(dz).copy())
            return original(ax, x, y, z, dx, dy, dz, *args, **kwargs)

        with patch.object(Axes3D, "bar3d", observe):
            fig = build_case_figure([("U1", 0.2, data)])
        try:
            use = abs(data["eigen_sites"]) <= 12
            np.testing.assert_array_equal(
                np.column_stack(heights), data["eigen_profiles"][use]
            )
            np.testing.assert_array_equal(
                fig.axes[1].images[0].get_array(),
                data["probability"][:, abs(data["sites"]) <= 30],
            )
            self.assertEqual(fig.axes[1].images[0].get_clim(), (0, 1))
            for key in data:
                np.testing.assert_array_equal(data[key], before[key])
        finally:
            plt.close(fig)

    def test_uniform_camera_height_limits_and_square_heatmaps(self):
        cases = [
            (model, gamma, fixture())
            for model, gamma in [
                ("U1", 0.01),
                ("U2", 0.95),
                ("U1", 0.2),
                ("U2", 0.7),
                ("U1", 0.3),
                ("U2", 0.1),
            ]
        ]
        fig = build_case_figure(cases)
        try:
            self.assertEqual(len(fig.axes), 12)
            for profile, evolution in zip(fig.axes[::2], fig.axes[1::2]):
                self.assertEqual((profile.elev, profile.azim), (22, -55))
                self.assertEqual(profile.get_zlim(), (0, 1.08))
                bbox = evolution.get_position()
                self.assertAlmostEqual(bbox.width * fig.get_figwidth() * 2.54, 2.5)
                self.assertAlmostEqual(bbox.height * fig.get_figheight() * 2.54, 2.5)
        finally:
            plt.close(fig)

    def test_nearby_moduli_in_the_specific_u2_panel_remain_distinct(self):
        self.assertEqual(
            profile_labels("U2", 0.1, [9.8513, 10.0, 0.101509884, 0.100000255]),
            ["9.85", "10.00", "0.1015", "0.1000"],
        )
        self.assertEqual(
            profile_labels("U1", 0.2, [1.4079, 0.7103, 1.4689, 0.6808]),
            ["1.41", "0.71", "1.47", "0.68"],
        )

    def test_invalid_probability_is_rejected(self):
        data = fixture()
        data["eigen_profiles"][0, 0] = -0.1
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        try:
            with self.assertRaisesRegex(ValueError, "nonnegative"):
                draw_profile(ax, "U1", 0.2, data)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
