"""USAGE:
    Import and call ensure_headless_display() before creating a render window.

REQUIREMENTS:
    On Linux without DISPLAY, xvfb-run and Xvfb must be installed.
"""

import os
from pathlib import Path
import shutil
import sys


def ensure_headless_display() -> None:
    if os.environ.get("DISPLAY"):
        return

    xvfb_run = shutil.which("xvfb-run")
    if xvfb_run is None:
        raise SystemExit(
            "No graphical display is available and xvfb-run was not found. "
            "Install Xvfb to render headlessly."
        )

    os.execvpe(
        xvfb_run,
        [
            xvfb_run,
            "-a",
            sys.executable,
            str(Path(sys.argv[0]).resolve()),
            *sys.argv[1:],
        ],
        os.environ,
    )
