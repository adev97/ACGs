"""Edit these settings, then run this file to make one set of ACGs per session."""

# ---------- Edit these settings ----------
from pathlib import Path


BIN_MS = 1.0
WINDOW_MS = 500.0            # Each side of zero: 500 means -500 to +500 ms.
SESSION_NUMBERS = [1, 2, 3, 4, 5]  # Use [4] to process only session_04.
PROBE = "ProbeA"

RECORDING_FOLDER = Path(r"R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Data_by_mouse\m002\small-spread-neuropixel-config\2026-09-25\ad_session_2026_09_25")
OUTPUT_FOLDER = RECORDING_FOLDER / "session_acgs"  # Relative to this script, or an absolute path.
# Uses the same full-recording cluster_KSLabel.tsv good IDs as the existing ACGs.
# ----------------------------------------

from datetime import datetime
from pathlib import Path

from acgs.sessions import run_sessions


def main():
    recording = Path(RECORDING_FOLDER)
    output = Path(OUTPUT_FOLDER)
    if not output.is_absolute():
        output = Path(__file__).resolve().parent / output
    output = output / (
        f"bin{BIN_MS:g}ms_window{WINDOW_MS:g}ms_"
        f"{datetime.now():%Y%m%d_%H%M%S_%f}"
    )
    print(f"Bin size: {BIN_MS:g} ms | Window: -{WINDOW_MS:g} to +{WINDOW_MS:g} ms")
    print(f"Output: {output}")
    run_sessions(recording / "sessions", recording / PROBE / "kilosort",
                 SESSION_NUMBERS, output, BIN_MS, WINDOW_MS, PROBE)


if __name__ == "__main__":
    main()
