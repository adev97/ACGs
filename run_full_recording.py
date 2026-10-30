"""Edit the settings below, then run this file to generate full-recording ACGs."""

# ---------- Edit these settings ----------
from pathlib import Path


BIN_MS = 1.0                 # Width of each bin in milliseconds.
WINDOW_MS = 500.0             # Each side of zero: 15 means -15 to +15 ms.

RECORDING_FOLDER = Path(r"R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Data_by_mouse\m002\small-spread-neuropixel-config\2026-09-25\ad_session_2026_09_25")
KILOSORT_DIR = RECORDING_FOLDER / "ProbeA"
OUTPUT_FOLDER = RECORDING_FOLDER / "full_recording_acgs"    # Relative to this script, or an absolute path.
DURATION_S = None            # None reads the full duration from recording metadata.
# ----------------------------------------

from datetime import datetime
from pathlib import Path

from acgs.cli import main as run_acgs


def main():
    # A separate folder for each run preserves earlier results.
    output_root = Path(OUTPUT_FOLDER)
    if not output_root.is_absolute():
        output_root = Path(__file__).resolve().parent / output_root
    run_name = (
        f"full_recording_bin{BIN_MS:g}ms_window{WINDOW_MS:g}ms_"
        f"{datetime.now():%Y%m%d_%H%M%S_%f}"
    )
    output_dir = output_root / run_name
    arguments = [
        str(KILOSORT_DIR),
        "--bin-ms", str(BIN_MS),
        "--window-ms", str(WINDOW_MS),
        "--output", str(output_dir),
    ]
    if DURATION_S is not None:
        arguments.extend(["--duration-s", str(DURATION_S)])
    print(f"Bin size: {BIN_MS:g} ms | Window: -{WINDOW_MS:g} to +{WINDOW_MS:g} ms")
    print(f"Output: {output_dir}")
    run_acgs(arguments)


if __name__ == "__main__":
    main()
