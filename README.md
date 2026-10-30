# Kilosort ACGs

A small Python package for full-recording autocorrelograms. Defaults: **1 ms bins,
−15 to +15 ms lags, Kilosort `KSLabel == good`, firing-rate normalization**.
Support for the already separated sessions can be added using the same computation
function; this version runs the full recording only.

## Edit settings in a script

Open `run_full_recording.py` in the repository root and edit the settings at the top:

```python
BIN_MS = 1.0
WINDOW_MS = 15.0  # -15 to +15 ms
```

Run that file from your Python editor, or run `python run_full_recording.py`.
No command-line parameters are needed. The input recording path is also set at
its top. Each run saves to a new folder under `outputs/`, named with the bin size,
window, and timestamp, so you can rerun it without overwriting earlier results.
Output paths are relative to the script, regardless of the current working folder.
The same Kilosort good-cluster selection and firing-rate normalization are used.
Install the package dependencies once using the instructions below.

## Install and run

From this repository, using Python 3.10 or newer:

```powershell
python -m pip install -e .
python -m acgs "R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Data_by_mouse\m002\small-spread-neuropixel-config\2026-09-25\ad_session_2026_09_25\ProbeA" --output outputs/full_recording
```

The installed command `kilosort-acgs` accepts the same arguments. Both the `ProbeA`
folder and its `kilosort` subfolder are accepted. Optional flags:
`--bin-ms 1 --window-ms 15 --duration-s <seconds>`. An existing nonempty output
folder is protected against overwriting; choose a new output folder for another run.

## Cluster selection and time base

Read spike assignments from `spike_clusters.npy` and spike sample indices from
`spike_times.npy`. Select cluster IDs from **`cluster_KSLabel.tsv`** where
`KSLabel` is `good`. The separate Phy `group` labels are not used. For the specified
recording, this selects **131 clusters**, whereas the curated `group` labels contain
14 `good` entries. Do not substitute template IDs for cluster IDs.

Read the sample rate from `params.py` as literal data, without executing it.
For the specified recording it is 30,000 Hz. Calculate the observation duration
from `concat/binary.json` and the binary file size (subtracting any binary offset),
without reading the raw recording. Fall back to `params.py` binary metadata when
`concat/binary.json` is absent. If the raw binary is unavailable or has multiple
segments, supply the true observation duration using `--duration-s`; the final
spike time is not used as a duration estimate. All good-cluster spikes must fall
within that duration. This full-recording workflow treats time as one continuous
observation; it does not exclude pauses or session boundaries.

## Counting and normalization

There are **30 bins**, with edges −15, −14, …, 0, …, +14, +15 ms.
Zero lies between the two central bins; centers are −14.5, …, −0.5, +0.5, …, +14.5 ms.
For each distinct pair, calculate the absolute sample difference and count it on
both sides. Absolute-lag bins are `[0,1)`, `[1,2)`, …, `[14,15]` ms; the outer
endpoint is included. Mirror the positive counts to obtain exactly symmetric
negative counts. Each spike's pairing with itself is excluded. Distinct spikes
at the same sample are retained, with one count in each central bin. Bins must
be an integer number of samples, and the window an integer number of bins.

Let `N` be the cluster spike count, `T` the full recording duration in seconds,
`r = N/T` its mean firing rate in Hz, `dt` the bin width in seconds, and `tau` the
bin center in seconds. The saved normalized ACG is dimensionless:

```text
normalized_acg(tau) = pair_count(tau) / [r² × dt × (T − |tau|)]
```

Equivalently, divide counts by `N × dt` to obtain the conditional firing rate,
then divide by mean firing rate `r`, with a finite-duration overlap correction
`T/(T−|tau|)`. For an independent stationary train, the baseline is approximately
**1**; the plots mark it with a dashed line. Sampling discretization and the
inclusive outer endpoint can produce small deviations from a continuous-time
baseline. Normalization uses the full recording mean rate; rate changes over the
recording can affect the ACG. Empty clusters produce NaN normalized values.

## Outputs

- `clusters/cluster_<id>_acg.png`: individual normalized ACG plots.
- `overview.pdf` and `overview_XX.png`: 20 clusters per page, with independent y scales.
- `acgs.npz`: cluster IDs, bin edges/centers, raw counts, normalized ACGs,
  spike counts, firing rates, sample rate, and duration.
- `acgs.csv`: long-form raw and normalized counts for every cluster and bin.
- `clusters.csv`: cluster spike counts and mean firing rates.
- `metadata.json`: source paths, selected IDs, parameters, duration provenance,
  normalization formula, and bin convention.

```python
import numpy as np
data = np.load("outputs/full_recording/acgs.npz")
print(data["cluster_ids"], data["normalized"].shape)
```

Run validation with `python -m unittest discover -s tests -v`.
