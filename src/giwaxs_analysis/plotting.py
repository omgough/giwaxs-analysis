"""Plotting helpers for 2D detector frames, q-space maps, line cuts and heatmaps.

These wrap matplotlib so the notebooks can stay short. Every function
returns the ``(fig, ax)`` it drew on so callers can keep tweaking.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from .calibration import Calibration


def plot_detector(
    frame: "np.ndarray",
    calib: "Calibration | None" = None,
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "viridis",
) -> "tuple[Figure, Axes]":
    """Show a raw detector frame, optionally with mask overlaid.

    Set ``log=True`` (default) for log-intensity, which is almost always
    what you want for diffraction data.
    """
    raise NotImplementedError


def plot_qmap(
    frame: "np.ndarray",
    calib: "Calibration",
    *,
    ax: "Axes | None" = None,
    log: bool = True,
) -> "tuple[Figure, Axes]":
    """Re-grid a frame into (qz, qxy) and plot it."""
    raise NotImplementedError


def plot_1d(
    q: "np.ndarray",
    intensity: "np.ndarray",
    *,
    ax: "Axes | None" = None,
    label: str | None = None,
    log_y: bool = False,
) -> "tuple[Figure, Axes]":
    """Plot a 1D radial / sector profile."""
    raise NotImplementedError


def plot_heatmap(
    q: "np.ndarray",
    time: "np.ndarray",
    intensity_2d: "np.ndarray",
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "magma",
) -> "tuple[Figure, Axes]":
    """In-situ heatmap: q on x, time (or scan index) on y, intensity as colour."""
    raise NotImplementedError


def line_cut(
    cake_q: "np.ndarray",
    cake_chi: "np.ndarray",
    cake_intensity: "np.ndarray",
    *,
    direction: str = "qz",
    width: float = 0.05,
) -> "tuple[np.ndarray, np.ndarray]":
    """Take an in-plane or out-of-plane line cut through a cake.

    ``direction`` is one of ``"qz"`` (out-of-plane) or ``"qxy"``
    (in-plane). ``width`` is the chi-window half-width in degrees.
    """
    raise NotImplementedError
