"""ACGs for already separated sessions with local ADC-second spike times."""

from pathlib import Path

import numpy as np

from .cli import generate_acgs, good_cluster_ids, read_params


def load_session(session_dir, probe="ProbeA"):
    """Read session duration and local times; reject unknown clock conventions."""
    session_dir = Path(session_dir)
    info_path = session_dir / "session_info.npz"
    with np.load(info_path, allow_pickle=False) as info:
        duration_s = float(info["duration_s"])
        time_reference = str(info["time_reference"].item())
        if time_reference != "session_local_adc_seconds":
            raise ValueError(f"Unsupported session time reference: {time_reference}")
        if str(info["boundary_convention"].item()) != "[start, end)":
            raise ValueError("Expected session boundaries [start, end).")
        start_s = float(info["adc_global_start_s"])
        end_s = float(info["adc_global_end_s"])
        if not np.isfinite(duration_s) or duration_s <= 0 or not np.isclose(
            end_s - start_s, duration_s, rtol=0, atol=1e-8
        ):
            raise ValueError("Session duration disagrees with its recorded boundaries.")
        session_name = str(info["session_name"].item())
    neural = session_dir / "neural" / probe
    times = np.load(neural / "adc_spike_times.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    assignments = np.load(neural / "spike_clusters.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    if times.shape != assignments.shape:
        raise ValueError(f"Spike array lengths disagree in {session_dir.name}.")
    metadata = dict(
        session_folder=str(session_dir.resolve()), session_name=session_name,
        session_info_file=str(info_path.resolve()), duration_source=str(info_path.resolve()),
        spike_times_file=str((neural / "adc_spike_times.npy").resolve()),
        spike_clusters_file=str((neural / "spike_clusters.npy").resolve()),
        time_reference=time_reference, boundary_convention="[start, end)",
        adc_global_start_s=start_s, adc_global_end_s=end_s,
        internal_time_resolution_s=1e-9,
        normalization_scope="session spike count and full session duration",
    )
    return times, assignments, duration_s, metadata


def run_sessions(sessions_root, full_kilosort_dir, session_numbers, output_root,
                 bin_ms=1.0, window_ms=15.0, probe="ProbeA"):
    """Keep the full-recording good IDs stable across all requested sessions."""
    sessions_root, full_kilosort_dir = Path(sessions_root), Path(full_kilosort_dir)
    if not (full_kilosort_dir / "spike_times.npy").exists() and (full_kilosort_dir / "kilosort").is_dir():
        full_kilosort_dir = full_kilosort_dir / "kilosort"
    numbers = list(session_numbers)
    if not numbers or any(not isinstance(n, int) or isinstance(n, bool) or n < 1 for n in numbers):
        raise ValueError("SESSION_NUMBERS must be a nonempty list of positive integers.")
    if len(numbers) != len(set(numbers)):
        raise ValueError("SESSION_NUMBERS must not contain duplicates.")
    sample_rate = float(read_params(full_kilosort_dir / "params.py")["sample_rate"])
    cluster_ids = good_cluster_ids(full_kilosort_dir)
    # Check every input before starting output generation.
    for number in numbers:
        session_dir = sessions_root / f"session_{number:02d}"
        load_session(session_dir, probe)
        output = Path(output_root) / session_dir.name
        if output.exists() and any(output.iterdir()):
            raise FileExistsError(f"Output folder is not empty: {output}")
    results = []
    for number in numbers:
        session_dir = sessions_root / f"session_{number:02d}"
        times, assignments, duration_s, metadata = load_session(session_dir, probe)
        metadata.update(kilosort_dir=str(full_kilosort_dir.resolve()),
                        label_file=str((full_kilosort_dir / "cluster_KSLabel.tsv").resolve()))
        print(f"\n{session_dir.name}: {metadata['session_name']}", flush=True)
        result = generate_acgs(
            times, assignments, cluster_ids, sample_rate, duration_s,
            Path(output_root) / session_dir.name, bin_ms, window_ms,
            metadata=metadata, title=f"{session_dir.name} ACGs", times_in_seconds=True,
        )
        results.append(result)
    return results
