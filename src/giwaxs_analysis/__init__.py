"""giwaxs_analysis — tools for processing GIWAXS data.

This package collects routines for loading detector frames, applying a
pyFAI geometry calibration, performing 1D / 2D integrations, plotting
in-situ heatmaps, and fitting Bragg peaks (including Scherrer analysis).

The submodules are intentionally thin — each exposes the small set of
functions that the notebooks in `notebooks/` actually call. When in
doubt, look at the example notebooks for end-to-end usage.

Submodules
----------
calibration
    Load and validate pyFAI PONI files and detector masks.
io
    Discover and load detector frames (typically via fabio).
integration
    1D radial / azimuthal integrations and 2D q-space reshaping.
plotting
    2D detector images, q-space maps, line cuts, in-situ heatmaps.
peak_fitting
    Peak fitting (Gaussian / pseudo-Voigt) and Scherrer crystallite-size
    analysis.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("giwaxs-analysis")
except PackageNotFoundError:  # pragma: no cover — package not installed
    __version__ = "0.0.0+unknown"

__all__ = [
    "__version__",
    "calibration",
    "io",
    "integration",
    "plotting",
    "peak_fitting",
]

from . import calibration, integration, io, peak_fitting, plotting  # noqa: E402,F401
