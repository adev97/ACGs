"""Command line workflow for full-recording Kilosort outputs."""

import argparse
import ast
import csv
import json
from pathlib import Path

import numpy as np

from .core import autocorrelogram


def read_params(path):
    """Read literals without executing the Kilosort params.py file."""
    result = {}
    for node in ast.parse(path.read_text(encoding="utf-8-sig")).body:
        if isinstance(node, ast.Assign):
            try:
                value = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    result[target.id] = value
    return result


def recording_duration(ks, params):
    """Use the actual binary length, never the last detected spike."""
    binary_metadata = ks.parent / "concat" / "binary.json"
    if binary_metadata.exists():
        kwargs = json.loads(binary_metadata.read_text(encoding="utf-8"))["kwargs"]
        if not np.isclose(float(kwargs["sampling_frequency"]), params["sample_rate"]):
            raise ValueError("Binary metadata and Kilosort sample rates disagree.")
        channels = int(kwargs["num_channels"])
        dtype = np.dtype(kwargs["dtype"])
        offset = int(kwargs.get("file_offset", 0))
        paths = [Path(p) for p in kwargs["file_paths"]]
        paths = [p if p.is_absolute() else binary_metadata.parent / p for p in paths]
        if len(paths) != 1:
            raise ValueError("Multiple binary segments: provide --duration-s explicitly.")
        raw = paths[0]
    else:
        candidates = params.get("dat_path", [])
        if isinstance(candidates, str):
            candidates = [candidates]
        if len(candidates) != 1:
            raise ValueError("Cannot resolve recording binary: provide --duration-s.")
        raw = Path(candidates[0])
        if not raw.is_absolute():
            raw = ks / raw
        channels = int(params["n_channels_dat"])
        dtype = np.dtype(params["dtype"])
        offset = int(params.get("offset", 0))
    frame_bytes = channels * dtype.itemsize
    size = raw.stat().st_size - offset
    if size <= 0 or size % frame_bytes:
        raise ValueError("Recording binary size is not an integer number of frames.")
    return size // frame_bytes / params["sample_rate"], str(raw)


def good_cluster_ids(ks):
    # Kilosort's automatic labels, deliberately independent of Phy group edits.
    with (ks / "cluster_KSLabel.tsv").open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    ids = [int(row["cluster_id"]) for row in rows if row["KSLabel"].strip().lower() == "good"]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Expected a nonempty set of unique KSLabel=good cluster IDs.")
    return np.array(sorted(ids), dtype=np.int64)


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kilosort_dir", type=Path, help="Kilosort folder or its ProbeA parent")
    parser.add_argument("--output", type=Path, default=Path("outputs/full_recording"))
    parser.add_argument("--bin-ms", type=float, default=1.0)
    parser.add_argument("--window-ms", type=float, default=15.0, help="Maximum absolute lag")
    parser.add_argument("--duration-s", type=float, help="Full observation duration; normally read from binary size")
    args = parser.parse_args(argv)
    ks = args.kilosort_dir.resolve()
    if not (ks / "spike_times.npy").exists() and (ks / "kilosort").is_dir():
        ks = ks / "kilosort"
    params = read_params(ks / "params.py")
    sample_rate = float(params["sample_rate"])
    if args.duration_s is None:
        duration_s, duration_source = recording_duration(ks, params)
    else:
        duration_s, duration_source = args.duration_s, "explicit --duration-s"
    cluster_ids = good_cluster_ids(ks)
    times = np.load(ks / "spike_times.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    assignments = np.load(ks / "spike_clusters.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    if times.shape != assignments.shape:
        raise ValueError("Spike times and cluster assignments have different lengths.")
    if not np.issubdtype(times.dtype, np.integer) or not np.issubdtype(assignments.dtype, np.integer):
        raise ValueError("Spike times and cluster IDs must be integer arrays.")
    keep = np.isin(assignments, cluster_ids)
    selected_ids = np.asarray(assignments[keep])
    selected_times = np.asarray(times[keep])
    del times, assignments, keep
    order = np.argsort(selected_ids, kind="stable")
    selected_times = selected_times[order]
    selected_ids = selected_ids[order]
    starts = np.searchsorted(selected_ids, cluster_ids, side="left")
    ends = np.searchsorted(selected_ids, cluster_ids, side="right")
    del order, selected_ids
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise FileExistsError(f"Output folder is not empty: {args.output}. Choose a new --output.")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    print(f"{len(cluster_ids)} KSLabel=good clusters; {duration_s:.6f} s; {sample_rate:g} Hz", flush=True)
    all_counts, all_normalized, summaries = [], [], []
    for index, (cluster, start, end) in enumerate(zip(cluster_ids, starts, ends)):
        spikes = selected_times[start:end]
        edges, counts, normalized = autocorrelogram(spikes, sample_rate, duration_s, args.bin_ms, args.window_ms)
        all_counts.append(counts)
        all_normalized.append(normalized)
        summaries.append((int(cluster), len(spikes), len(spikes) / duration_s))
        print(f"[{index + 1}/{len(cluster_ids)}] cluster {cluster}: {len(spikes)} spikes", flush=True)
    all_counts, all_normalized = np.array(all_counts), np.array(all_normalized)
    centers = (edges[:-1] + edges[1:]) / 2
    np.savez_compressed(args.output / "acgs.npz", cluster_ids=cluster_ids, bin_edges_ms=edges,
                        lag_centers_ms=centers, counts=all_counts, normalized=all_normalized,
                        spike_counts=np.array([s[1] for s in summaries]),
                        firing_rates_hz=np.array([s[2] for s in summaries]),
                        sample_rate_hz=sample_rate, duration_s=duration_s)
    write_csv(args.output / "clusters.csv", ["cluster_id", "n_spikes", "firing_rate_hz"], summaries)
    write_csv(args.output / "acgs.csv", ["cluster_id", "lag_left_ms", "lag_right_ms", "lag_center_ms", "pair_count", "normalized_acg"],
              ((int(c), edges[j], edges[j+1], centers[j], int(all_counts[i,j]), all_normalized[i,j])
               for i,c in enumerate(cluster_ids) for j in range(len(centers))))
    metadata = dict(kilosort_dir=str(ks), label_file=str(ks / "cluster_KSLabel.tsv"), label="good",
                    cluster_count=len(cluster_ids), cluster_ids=cluster_ids.tolist(), sample_rate_hz=sample_rate,
                    duration_s=duration_s, duration_source=duration_source, bin_ms=args.bin_ms,
                    window_ms=args.window_ms, self_pairs_excluded=True,
                    normalization="counts / (rate_hz**2 * bin_width_s * (duration_s - abs(lag_center_s)))",
                    bin_convention="Mirror absolute-lag bins, left-inclusive and right-exclusive; final bin includes window endpoint; zero is a boundary.")
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    def draw(ax, i):
        cluster, n_spikes, rate = summaries[i]
        ax.stairs(all_normalized[i], edges, fill=True, color="#286f9b", linewidth=0.7)
        ax.axhline(1, color="#bb6b39", linestyle="--", linewidth=0.8)
        ax.axvline(0, color="0.5", linewidth=0.5)
        ax.set(xlim=(-args.window_ms, args.window_ms), ylim=(0, None),
               title=f"Cluster {cluster} | {rate:.2f} Hz | {n_spikes:,} spikes",
               xlabel="Lag (ms)", ylabel="Rate-normalized ACG")
        ax.spines[["top", "right"]].set_visible(False)

    plot_dir = args.output / "clusters"
    plot_dir.mkdir()
    for i, cluster in enumerate(cluster_ids):
        fig, ax = plt.subplots(figsize=(6, 3.6), layout="constrained")
        draw(ax, i)
        fig.savefig(plot_dir / f"cluster_{cluster}_acg.png", dpi=160)
        plt.close(fig)
    with PdfPages(args.output / "overview.pdf") as pdf:
        for page, first in enumerate(range(0, len(cluster_ids), 20), start=1):
            fig, axes = plt.subplots(5, 4, figsize=(16, 17), layout="constrained")
            fig.suptitle(f"Full-recording ACGs | KSLabel=good | {args.bin_ms:g} ms bins | page {page}")
            for i, ax in enumerate(axes.flat, start=first):
                if i < min(first + 20, len(cluster_ids)):
                    draw(ax, i)
                    ax.title.set_fontsize(9)
                else:
                    ax.set_visible(False)
            pdf.savefig(fig)
            fig.savefig(args.output / f"overview_{page:02d}.png", dpi=110)
            plt.close(fig)
    print(f"Saved plots and numerical results to {args.output.resolve()}", flush=True)
