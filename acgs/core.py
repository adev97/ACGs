"""Pair counting and firing-rate normalization, independent of file format."""

import numpy as np


def autocorrelogram(spike_samples, sample_rate, duration_s, bin_ms=1.0, window_ms=15.0):
    """Return (edges_ms, counts, normalized).

    Count distinct spike pairs in both directions, excluding self pairs.
    Bins on each side use absolute lags [0, bin), [bin, 2*bin), etc.;
    the final bin includes the window endpoint. Positive counts are mirrored
    exactly. Distinct spikes at the same sample contribute one count to each
    central bin. The origin is a bin boundary, not a dedicated zero-lag bin.

    Normalize by rate**2 * bin_width_s * (duration_s - abs(bin_center_s)).
    This corrects the finite observation overlap at each lag. An independent
    stationary spike train has an expected value approximately equal to 1.
    """
    values = np.asarray(spike_samples).reshape(-1)
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError("Spike times must be integer sample indices.")
    if not all(np.isfinite(v) and v > 0 for v in (sample_rate, duration_s, bin_ms, window_ms)):
        raise ValueError("Sample rate, duration, bin width, and window must be positive and finite.")
    if duration_s <= window_ms / 1000:
        raise ValueError("Recording duration must exceed the lag window.")
    n_side = int(round(window_ms / bin_ms))
    if n_side < 1 or not np.isclose(n_side * bin_ms, window_ms):
        raise ValueError("Window must be an integer multiple of bin width.")
    step = bin_ms * sample_rate / 1000
    limit = window_ms * sample_rate / 1000
    if not np.isclose(step, round(step), rtol=0, atol=1e-7):
        raise ValueError("Bin width must correspond to an integer number of samples.")
    step, limit = int(round(step)), int(round(limit))
    if step < 1:
        raise ValueError("Bin width must be at least one sample.")
    spikes = np.sort(values.astype(np.int64, copy=False))
    if spikes.size and (spikes[0] < 0 or spikes[-1] >= duration_s * sample_rate):
        raise ValueError("Spike samples fall outside the recording duration.")
    positive = np.zeros(n_side, dtype=np.int64)
    # Only visit neighbors that could fall within the window; no N x N matrix.
    for offset in range(1, len(spikes)):
        delays = spikes[offset:] - spikes[:-offset]
        delays = delays[delays <= limit]
        if not len(delays):
            break
        bins = np.minimum(delays // step, n_side - 1)
        positive += np.bincount(bins, minlength=n_side)
    counts = np.concatenate((positive[::-1], positive))
    edges = np.arange(-n_side, n_side + 1, dtype=float) * bin_ms
    centers_s = (edges[:-1] + edges[1:]) / 2000
    rate_hz = len(spikes) / duration_s
    expected = rate_hz**2 * (bin_ms / 1000) * (duration_s - np.abs(centers_s))
    normalized = np.full(counts.shape, np.nan)
    np.divide(counts, expected, out=normalized, where=expected > 0)
    return edges, counts, normalized
