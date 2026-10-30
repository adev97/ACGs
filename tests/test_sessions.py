import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from acgs.cli import main as full_main
from acgs.core import autocorrelogram_seconds
from acgs.sessions import load_session, run_sessions


class SessionTests(unittest.TestCase):
    def test_second_precision_bin_edges_and_validation(self):
        _, counts, _ = autocorrelogram_seconds([0, .000999999, .001], 1)
        self.assertEqual(counts[15], 2)
        self.assertEqual(counts[16], 1)
        for times in ([float("nan")], [-.1], [1]):
            with self.assertRaises(ValueError):
                autocorrelogram_seconds(times, 1)

    def test_independent_session_exports_preserve_ids_duration_and_silent_units(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            ks = root / "ProbeA" / "kilosort"
            ks.mkdir(parents=True)
            (ks / "params.py").write_text("sample_rate = 30000\n")
            (ks / "cluster_KSLabel.tsv").write_text("cluster_id\tKSLabel\n7\tgood\n9\tgood\n12\tmua\n")
            (ks / "cluster_group.tsv").write_text("cluster_id\tgroup\n12\tgood\n")
            np.save(ks / "spike_times.npy", np.array([30], dtype=np.int64))
            np.save(ks / "spike_clusters.npy", np.array([7], dtype=np.int64))
            for number, second in ((1, .999), (2, .001)):
                session = root / "sessions" / f"session_{number:02d}"
                neural = session / "neural" / "ProbeA"
                neural.mkdir(parents=True)
                np.savez(session / "session_info.npz", duration_s=1.0,
                         adc_global_start_s=float(number - 1), adc_global_end_s=float(number),
                         time_reference="session_local_adc_seconds", boundary_convention="[start, end)",
                         session_name=f"test {number}")
                np.save(neural / "adc_spike_times.npy", [second, .2])
                np.save(neural / "spike_clusters.npy", np.array([7, 12]))
            outputs = run_sessions(root / "sessions", ks, [1, 2], root / "out")
            for output in outputs:
                with np.load(output / "acgs.npz") as saved:
                    np.testing.assert_array_equal(saved["cluster_ids"], [7, 9])
                    np.testing.assert_array_equal(saved["spike_counts"], [1, 0])
                    np.testing.assert_array_equal(saved["firing_rates_hz"], [1, 0])
                    self.assertEqual(saved["counts"].sum(), 0)  # No pair across session boundaries.
                    self.assertTrue(np.isnan(saved["normalized"][1]).all())
                    self.assertEqual(saved["duration_s"], 1)
                meta = json.loads((output / "metadata.json").read_text())
                self.assertEqual(meta["time_reference"], "session_local_adc_seconds")
                self.assertTrue(meta["label_file"].endswith("cluster_KSLabel.tsv"))
                self.assertTrue((output / "overview.pdf").exists())
            # Existing full-recording CLI still uses the original sample-index path.
            full_main([str(ks), "--duration-s", "2", "--output", str(root / "full")])
            with np.load(root / "full" / "acgs.npz") as saved:
                self.assertEqual(saved["firing_rates_hz"][0], .5)
            # Reject clock conventions instead of interpreting global times as local.
            info = root / "sessions" / "session_01" / "session_info.npz"
            with np.load(info) as data:
                payload = {k: data[k] for k in data.files}
            payload["time_reference"] = "global_seconds"
            np.savez(info, **payload)
            with self.assertRaises(ValueError):
                load_session(info.parent)


if __name__ == "__main__":
    unittest.main()
