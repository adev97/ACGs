# Kilosort ACGs

Create double-sided, firing-rate-normalized autocorrelograms for a full recording
or for sessions that have already been separated. Edit the settings at the top of
the appropriate runner and run it from your Python editor. No command-line
parameters are required for either runner.

Both runners are currently configured for **1 ms bins and -500 to +500 ms lags**.
Both use the same full-recording **Kilosort `KSLabel == good` cluster IDs**
(131 clusters for the configured recording). No smoothing is applied.

## 1. Install once

Use Python 3.10 or newer. In a terminal, open this repository and install the package:

```powershell
cd "R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Projects\ACGs"
python -m pip install -e .
```

This installs NumPy, Matplotlib, and the local package. Use the same Python
environment in your editor when running the scripts.

## 2. Full-recording ACGs

Open `run_full_recording.py`. Edit the settings at the top:

```python
from pathlib import Path

BIN_MS = 1.0
WINDOW_MS = 500.0  # Maximum lag on EACH side of zero.

RECORDING_FOLDER = Path(
    r"R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Data_by_mouse\m002\small-spread-neuropixel-config\2026-09-25\ad_session_2026_09_25"
)
KILOSORT_DIR = RECORDING_FOLDER / "ProbeA"
OUTPUT_FOLDER = RECORDING_FOLDER / "full_recording_acgs"
DURATION_S = None
```

`RECORDING_FOLDER` must be a `Path` object before using `/` to append folders.
Import `Path` above these settings. The input can be either `ProbeA` or
`ProbeA/kilosort`; the loader automatically finds the Kilosort subfolder.

Run the file using your editor's Run button, or from the repository terminal:

```powershell
python run_full_recording.py
```

The current output setting creates:

```text
<RECORDING_FOLDER>/
  full_recording_acgs/
    full_recording_bin1ms_window500ms_<timestamp>/
      clusters/
      overview.pdf
      overview_01.png ...
      acgs.npz
      acgs.csv
      clusters.csv
      metadata.json
```

Leave `DURATION_S = None` to obtain the full duration from recording metadata
and the raw binary's file size. The raw recording is not loaded. If the binary
is unavailable, set `DURATION_S` to the known complete recording duration in
seconds, including silent periods. The last spike time is not a duration estimate.

The full recording is treated as one continuous observation, including pairs
that span a boundary between the separated sessions.

## 3. Individual-session ACGs

Open `run_session_acgs.py`. Set `RECORDING_FOLDER` to the same recording root
shown above, then edit these settings:

```python
BIN_MS = 1.0
WINDOW_MS = 500.0
SESSION_NUMBERS = [1, 2, 3, 4, 5]
PROBE = "ProbeA"
OUTPUT_FOLDER = RECORDING_FOLDER / "session_acgs"
```

To process only the receptive-field mapping session:

```python
SESSION_NUMBERS = [4]
```

To process any subset, list the desired session numbers, for example `[1, 3, 5]`.
Run the file from your editor, or use:

```powershell
python run_session_acgs.py
```

With all five sessions selected, the current output setting creates:

```text
<RECORDING_FOLDER>/
  session_acgs/
    bin1ms_window500ms_<timestamp>/
      session_01/
      session_02/
      session_03/
      session_04/
      session_05/
```

Each session folder contains the same individual PNG plots, overview PDF/PNGs,
and numerical files listed for the full recording. Session numbering matches
the existing `sessions/session_XX` input folders.

Sessions are processed independently: **no pair crosses a session boundary**.
Normalization uses each cluster's session spike count and the entire session
duration from `session_info.npz`, including silent beginning and end periods.
The same selected cluster IDs are retained in every session. A silent cluster
has zero counts and firing rate, NaN normalized values, and a no-spikes plot note.

## Changing bin size, lag window, and output location

Edit `BIN_MS` and `WINDOW_MS` in the runner you plan to use. For example:

```python
BIN_MS = 1.0
WINDOW_MS = 15.0  # -15 to +15 ms.
```

The window must be an integer multiple of the bin width. For full-recording
sample times, each bin must also represent an integer number of acquisition
samples. The observation duration must exceed the maximum lag.

`OUTPUT_FOLDER` may be an absolute path or a path relative to the runner's
directory. With the settings above, both runners save alongside the recording.
Each run creates a fresh timestamped subfolder so earlier runs are preserved.
The underlying workflow refuses to overwrite a nonempty result folder.

## Input files and cluster identity

Paths below are relative to `RECORDING_FOLDER`:

| Input | Used for |
| --- | --- |
| `ProbeA/kilosort/cluster_KSLabel.tsv` | Select `cluster_id` rows where `KSLabel` is `good`, for both workflows. |
| `ProbeA/kilosort/spike_times.npy` | Full-recording spike times in sample indices. |
| `ProbeA/kilosort/spike_clusters.npy` | Full-recording cluster ID for each spike. |
| `ProbeA/kilosort/params.py` | Acquisition sampling rate, read as literals without executing Python. |
| `ProbeA/concat/binary.json` and its referenced raw binary | Full-recording duration from metadata and binary size. Falls back to `params.py` binary metadata when `binary.json` is absent. |
| `sessions/session_XX/session_info.npz` | Session duration, ADC boundaries, and time-reference convention. |
| `sessions/session_XX/neural/ProbeA/adc_spike_times.npy` | Synchronized session-local spike times in seconds. |
| `sessions/session_XX/neural/ProbeA/spike_clusters.npy` | Session cluster ID for each spike. |

The spike-time and cluster arrays align entry by entry. Selection uses
`cluster_KSLabel.tsv`, not the separate curated `cluster_group.tsv` labels.
Do not substitute template IDs or sequential unit indices for cluster IDs.
For example, `cluster_254_acg.png` belongs to actual cluster 254.

Session metadata must declare `session_local_adc_seconds` and `[start, end)`
boundaries. `spike_times_global.npy` is not used as local session time. ADC
seconds are represented internally as integer nanoseconds rounded to the nearest
1 ns; they are not rounded back to the acquisition sample grid. In exported
files, `sample_rate_hz` still records the acquisition rate (30,000 Hz here).

## Counting and firing-rate normalization

The ACG is symmetric about zero. With 1 ms bins and a 15 ms window there are
30 bins spanning -15 to +15 ms; with a 500 ms window there are 1,000 bins.
Zero is a boundary between the two central bins, not a dedicated zero-lag bin.
Absolute-lag bins are left-inclusive and right-exclusive, with the outer window
endpoint included in the final bin. Positive counts are mirrored to negative lags.

Each spike's pairing with itself is excluded. Distinct spikes at an identical
time are retained, contributing one count to each central bin. Plots show the
normalized bin counts directly, without Gaussian filtering, moving averages,
or interpolation.

For each cluster, let `N` be its spike count in the observation interval, `T` the
interval duration in seconds, `r = N/T` its mean firing rate in Hz, `dt` the bin
width in seconds, and `tau` the bin center in seconds:

```text
normalized_acg(tau) = pair_count(tau) / [r^2 * dt * (T - |tau|)]
```

This divides counts by `N * dt` to obtain conditional firing rate, divides
again by the mean firing rate `r`, and applies the finite-duration correction
`T / (T - |tau|)`. For session ACGs, both `N` and `T` are session-specific.

The dashed line at **1** is the expected baseline for independent stationary
spikes. It is not a statistical significance threshold. Values above/below 1
indicate more/fewer pairs than that baseline. Rate changes within an interval,
sampling discretization, and the inclusive outer endpoint can affect the baseline.

## Result files

| File | Contents |
| --- | --- |
| `clusters/cluster_<id>_acg.png` | Individual normalized ACG, actual cluster ID, spike count, and mean firing rate. |
| `overview.pdf`, `overview_XX.png` | Overview pages, up to 20 clusters per page, with independent y scales. |
| `acgs.npz` | Cluster IDs, bin edges/centers, raw counts, normalized values, spike counts, rates, sampling rate, and duration. |
| `acgs.csv` | One row per cluster and bin, with raw and normalized counts. |
| `clusters.csv` | Cluster IDs, spike counts, and mean firing rates. |
| `metadata.json` | Parameters, selected IDs, source paths, time reference, normalization, and duration provenance. |

Load numerical results by replacing the path below with a completed run folder:

```python
from pathlib import Path
import numpy as np

result_folder = Path(r"R:\path\to\a\completed\run\session_04")
data = np.load(result_folder / "acgs.npz")
print(data["cluster_ids"])
print(data["normalized"].shape)  # (number of clusters, number of bins)
```

## Optional full-recording command line

The script settings above are the normal workflow. The full-recording package
also accepts command-line arguments; they are independent of the runner settings:

```powershell
python -m acgs "R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Data_by_mouse\m002\small-spread-neuropixel-config\2026-09-25\ad_session_2026_09_25\ProbeA" --bin-ms 1 --window-ms 500 --output outputs/full_recording_500ms
```

The installed `kilosort-acgs` command accepts the same arguments. If omitted,
the command-line defaults are 1 ms bins and a 15 ms window. `--duration-s`
overrides the full-recording duration when needed. Choose a new output folder
for another command-line run.

## Validation

From the repository, run:

```powershell
python -m unittest discover -s tests -v
```

Checks cover pair counting, bin boundaries, normalization, cluster-label selection,
session timing, silent units, independent session processing, and the full-recording
workflow.
