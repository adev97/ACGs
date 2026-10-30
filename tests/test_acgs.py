import tempfile
import unittest
from pathlib import Path

import numpy as np

from acgs import autocorrelogram
from acgs.cli import good_cluster_ids, read_params, recording_duration


class ACGTests(unittest.TestCase):
    def test_against_brute_force_including_boundaries_and_coincident_spikes(self):
        spikes = np.array([170, 0, 30, 30, 450, 900], dtype=np.int64)
        edges, counts, norm = autocorrelogram(spikes, 30000, 1)
        expected = np.zeros(30, dtype=np.int64)
        for i in range(len(spikes)):
            for j in range(i + 1, len(spikes)):
                delay = abs(int(spikes[i]) - int(spikes[j]))
                if delay <= 450:
                    b = min(delay // 30, 14)
                    expected[15 + b] += 1
                    expected[14 - b] += 1
        np.testing.assert_array_equal(counts, expected)
        np.testing.assert_array_equal(counts, counts[::-1])
        np.testing.assert_array_equal(edges, np.arange(-15, 16))
        centers = (edges[:-1] + edges[1:]) / 2000
        np.testing.assert_allclose(norm, counts / (36 * .001 * (1 - abs(centers))))

    def test_self_pairs_excluded_and_empty_train(self):
        _, counts, norm = autocorrelogram(np.array([0]), 30000, 1)
        self.assertEqual(counts.sum(), 0)
        self.assertEqual(norm.sum(), 0)
        _, counts, norm = autocorrelogram(np.array([], dtype=int), 30000, 1)
        self.assertEqual(counts.sum(), 0)
        self.assertTrue(np.isnan(norm).all())

    def test_independent_spikes_normalize_to_one(self):
        rng = np.random.default_rng(5)
        spikes = np.sort(rng.integers(0, 30000 * 1000, size=50000))
        _, _, norm = autocorrelogram(spikes, 30000, 1000)
        self.assertAlmostEqual(float(norm.mean()), 1, delta=.04)

    def test_invalid_inputs(self):
        for kwargs in ({"sample_rate": 0}, {"duration_s": .001}, {"bin_ms": .333}, {"window_ms": 15.3}):
            params = dict(sample_rate=30000, duration_s=1)
            params.update(kwargs)
            with self.assertRaises(ValueError):
                autocorrelogram(np.array([0]), **params)
        with self.assertRaises(ValueError):
            autocorrelogram(np.array([30000]), 30000, 1)

    def test_automatic_labels_and_duration_from_binary(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            ks = root / "kilosort"
            ks.mkdir()
            (ks / "cluster_KSLabel.tsv").write_text("cluster_id\tKSLabel\n1\tgood\n2\tmua\n")
            (ks / "cluster_group.tsv").write_text("cluster_id\tgroup\n2\tgood\n")
            np.testing.assert_array_equal(good_cluster_ids(ks), [1])
            raw = ks / "test.raw"
            raw.write_bytes(bytes(400))
            params_path = ks / "params.py"
            params_path.write_text("sample_rate = 100\nn_channels_dat = 2\ndtype = 'int16'\ndat_path = ['test.raw']\n")
            params = read_params(params_path)
            duration, _ = recording_duration(ks, params)
            self.assertEqual(duration, 1)


if __name__ == "__main__":
    unittest.main()
