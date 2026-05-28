"""Plotting helpers for 2D detector frames, q-space maps, line cuts and heatmaps.

These wrap matplotlib so the notebooks can stay short. Every function
returns the ``(fig, ax)`` it drew on so callers can keep tweaking.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from .calibration import Calibration
    from .peak_fitting import DoublePeakFit, PeakFit


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


def plot_peak_fit(
    q: "np.ndarray",
    intensity: "np.ndarray",
    fit: "PeakFit",
    *,
    q_range: tuple[float, float] | None = None,
    title: str | None = None,
    show_residuals: bool = True,
) -> "tuple[Figure, list[Axes]]":
    """Two-panel plot of a peak fit: data + model on top, residuals on bottom.

    Parameters
    ----------
    q, intensity
        The 1D profile that was fit.
    fit
        The :class:`giwaxs_analysis.peak_fitting.PeakFit` result.
    q_range
        Optional ``(q_lo, q_hi)`` to draw as vertical dashed lines,
        showing the window the fit was performed in.
    title
        Plot title. Defaults to a string built from the fit label and model.
    show_residuals
        If False, only the top panel is drawn.

    Returns
    -------
    (fig, axes)
        ``axes`` is a list of one or two Axes depending on ``show_residuals``.
    """
    import matplotlib.pyplot as plt

    from .peak_fitting import linear_gaussian, linear_lorentzian

    q = np.asarray(q, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    model_fn = linear_gaussian if fit.model == "gaussian" else linear_lorentzian
    fitted = model_fn(
        q,
        fit.background_slope,
        fit.background_intercept,
        fit.amplitude,
        fit.q_center,
        fit.sigma,
    )

    if show_residuals:
        fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1]})
        ax_main, ax_res = axes
    else:
        fig, ax_main = plt.subplots(figsize=(10, 4))
        ax_res = None
        axes = [ax_main]

    label_str = fit.label if fit.label is not None else ""
    ax_main.plot(q, intensity, "b-", label=f"data {label_str}".strip())
    ax_main.plot(q, fitted, "r-", label=f"{fit.model} fit")
    if q_range is not None:
        ax_main.axvline(q_range[0], color="k", linestyle="--", linewidth=0.8)
        ax_main.axvline(q_range[1], color="k", linestyle="--", linewidth=0.8)
    ax_main.set_ylabel("Intensity")
    ax_main.legend()
    if title is None:
        title = f"Peak fit — {fit.label or ''} ({fit.model})".strip()
    ax_main.set_title(title)

    if ax_res is not None:
        residuals = intensity - fitted
        ax_res.plot(q, residuals, "g-")
        ax_res.axhline(0, color="k", linestyle="--", linewidth=0.8)
        ax_res.set_xlabel(r"$q$ (Å$^{-1}$)")
        ax_res.set_ylabel("Residuals")
    else:
        ax_main.set_xlabel(r"$q$ (Å$^{-1}$)")

    fig.tight_layout()
    return fig, axes


def plot_double_peak_fit(
    q: "np.ndarray",
    intensity: "np.ndarray",
    fit: "DoublePeakFit",
    *,
    q_range: tuple[float, float] | None = None,
    title: str | None = None,
) -> "tuple[Figure, Axes]":
    """Plot a two-Lorentzian fit, with each component drawn separately.

    Each Lorentzian component is shown stacked on the linear background,
    matching the style of the original notebook.
    """
    import matplotlib.pyplot as plt

    from .peak_fitting import linear_double_lorentzian

    q = np.asarray(q, dtype=float)
    m, c, A1, mu1, g1, A2, mu2, g2 = fit.raw_params

    linear = m * q + c
    lor1 = A1 * g1 ** 2 / ((q - mu1) ** 2 + g1 ** 2)
    lor2 = A2 * g2 ** 2 / ((q - mu2) ** 2 + g2 ** 2)
    combined = linear_double_lorentzian(q, *fit.raw_params)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(q, intensity, "b-", label=f"data {fit.label or ''}".strip())
    ax.plot(q, combined, "r-", label="double-Lorentzian fit")
    ax.plot(q, linear, "k--", linewidth=0.8, label="linear background")
    ax.plot(q, lor1 + linear, "g--", linewidth=0.8, label="Lorentzian 1")
    ax.plot(q, lor2 + linear, "m--", linewidth=0.8, label="Lorentzian 2")
    if q_range is not None:
        ax.axvline(q_range[0], color="k", linestyle=":", linewidth=0.8)
        ax.axvline(q_range[1], color="k", linestyle=":", linewidth=0.8)

    ax.set_xlabel(r"$q$ (Å$^{-1}$)")
    ax.set_ylabel("Intensity")
    if title is None:
        title = f"Double-Lorentzian fit — {fit.label or ''}".strip()
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig, ax
