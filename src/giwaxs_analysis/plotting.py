"""Plotting helpers for 2D detector frames, q-space maps, line cuts and heatmaps.

These wrap matplotlib so the notebooks can stay short. Every plotting
function returns the ``(fig, ax)`` (or ``(fig, list_of_ax)``) it drew on
so callers can keep tweaking.

Implemented
-----------
* ``plot_detector`` — single 2D detector frame, optionally with mask applied.
* ``plot_detector_grid`` — grid of frames sharing a colour scale.
* ``plot_cake`` — (q, χ) cake from :func:`giwaxs_analysis.integration.cake`.
* ``plot_peak_fit`` / ``plot_double_peak_fit`` — fit result + residuals.

Stubs (not implemented)
-----------------------
* ``plot_qmap`` — GIWAXS qz / qxy reshape (needs incidence angle).
* ``plot_1d`` / ``plot_heatmap`` — currently inlined in the notebooks.
* ``line_cut`` — extract in-/out-of-plane cuts from a cake.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from .calibration import Calibration
    from .peak_fitting import DoublePeakFit, PeakFit


# ---------------------------------------------------------------------------
# Detector frames
# ---------------------------------------------------------------------------


def plot_detector(
    frame: np.ndarray,
    calib: "Calibration | None" = None,
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    title: str | None = None,
) -> "tuple[Figure, Axes]":
    """Show a raw detector frame.

    Parameters
    ----------
    frame
        2D detector image, shape ``(H, W)``.
    calib
        Optional :class:`Calibration`. If given, masked pixels are
        rendered as NaN (transparent) instead of skewing the colour scale.
    ax
        Existing axes to draw on. If None, a new figure is created.
    log
        Log intensity scale (almost always what you want for diffraction).
    cmap
        Matplotlib colormap name.
    vmin, vmax
        Intensity range. None = auto from the data.
    title
        Plot title.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    frame = np.asarray(frame, dtype=float)
    if calib is not None:
        frame = np.where(calib.mask, np.nan, frame)

    if vmin is None:
        vmin = max(np.nanmin(frame), 1.0) if log else np.nanmin(frame)
    if vmax is None:
        vmax = np.nanmax(frame)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    else:
        fig = ax.figure

    norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(frame, cmap=cmap, norm=norm)
    ax.set_xlabel("pixel x")
    ax.set_ylabel("pixel y")
    if title:
        ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Intensity (counts)")

    return fig, ax


def plot_detector_grid(
    frames: "np.ndarray | list[np.ndarray]",
    calib: "Calibration | None" = None,
    *,
    labels: list[str] | None = None,
    ncols: int = 3,
    log: bool = True,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
) -> "tuple[Figure, list[Axes]]":
    """Plot multiple detector frames in a grid with a shared colour scale.

    Useful for eyeballing a whole folder of frames at once.

    Parameters
    ----------
    frames
        Either a 3D array ``(N, H, W)`` or a list of 2D arrays.
    calib
        Optional :class:`Calibration` — masked pixels become NaN.
    labels
        Per-panel titles. Defaults to ``"Frame 0", "Frame 1", ...``.
    ncols
        Number of columns in the grid.
    log, cmap, vmin, vmax
        Same as :func:`plot_detector`. ``vmin``/``vmax`` default to the
        full-stack min/max so the colour scale is comparable across panels.

    Returns
    -------
    (fig, axes) : Figure and list of Axes.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    stack = np.asarray([np.asarray(f, dtype=float) for f in frames])
    if calib is not None:
        stack = np.where(calib.mask, np.nan, stack)

    if vmin is None:
        vmin = max(np.nanmin(stack), 1.0) if log else np.nanmin(stack)
    if vmax is None:
        vmax = np.nanmax(stack)
    norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)

    n = len(stack)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False
    )
    flat_axes = axes.flatten()

    im = None
    for i, ax in enumerate(flat_axes):
        if i < n:
            im = ax.imshow(stack[i], cmap=cmap, norm=norm)
            label = labels[i] if labels else f"Frame {i}"
            ax.set_title(label)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.axis("off")

    if im is not None:
        fig.colorbar(im, ax=flat_axes[:n], shrink=0.6, label="Intensity (counts)")

    fig.tight_layout()
    return fig, list(flat_axes)


# ---------------------------------------------------------------------------
# Cake (q, χ)
# ---------------------------------------------------------------------------


def plot_cake(
    q: np.ndarray,
    chi: np.ndarray,
    I_2d: np.ndarray,
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    title: str | None = None,
) -> "tuple[Figure, Axes]":
    """Plot a (q, χ) cake from :func:`giwaxs_analysis.integration.cake`.

    If geometry is right, scattering rings appear as vertical lines.
    Curved or tilted lines = something's off (sample tilt, wrong PONI,
    incorrect detector distance...).

    Parameters
    ----------
    q, chi, I_2d
        Output of :func:`giwaxs_analysis.integration.cake`. ``I_2d`` has
        shape ``(len(chi), len(q))``.
    log, cmap, vmin, vmax, ax, title
        Same as :func:`plot_detector`.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    I_plot = np.where(I_2d > 0, I_2d, np.nan) if log else np.asarray(I_2d, dtype=float)

    if vmin is None:
        vmin = max(np.nanmin(I_plot), 1.0) if log else np.nanmin(I_plot)
    if vmax is None:
        vmax = np.nanmax(I_plot)

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5))
    else:
        fig = ax.figure

    norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(
        I_plot,
        origin="lower",
        aspect="auto",
        extent=(q.min(), q.max(), chi.min(), chi.max()),
        cmap=cmap,
        norm=norm,
    )
    ax.set_xlabel(r"$q$ (Å$^{-1}$)")
    ax.set_ylabel(r"$\chi$ (°)")
    if title:
        ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Intensity (counts)")

    return fig, ax


# ---------------------------------------------------------------------------
# Stubs — not yet implemented
# ---------------------------------------------------------------------------

def plot_qmap(
    qxy: np.ndarray,
    qz: np.ndarray,
    I_2d: np.ndarray,
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "RdYlBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    qxy_range: tuple[float, float] | None = None,
    qz_range: tuple[float, float] | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (8, 4),
    masked_color: str | tuple | None = None,
) -> "tuple[Figure, Axes]":
    """Plot a GIWAXS q-map from :func:`giwaxs_analysis.integration.giwaxs_reshape`.

    ``qxy`` on x (in-plane), ``qz`` on y (out-of-plane).

    Parameters
    ----------
    qxy, qz, I_2d
        Output of :func:`giwaxs_analysis.integration.giwaxs_reshape`.
    log
        Log intensity scale (almost always what you want for diffraction).
    cmap
        Matplotlib colormap name. Default ``"RdYlBu_r"`` gives a
        dark-blue → cream → red look. Try ``"jet"``, ``"turbo"``,
        ``"nipy_spectral"`` for alternatives.
    vmin, vmax
        Intensity range for the colour scale. ``None`` = auto from data.
    qxy_range, qz_range
        Optional ``(min, max)`` tuples (Å⁻¹) to crop the displayed area.
    title
        Plot title.
    figsize
        Figure size when ``ax`` is None. Default is wider than tall to
        match typical GIWAXS aspect (qxy span > qz span).
    masked_color
        Colour for masked / NaN pixels. ``None`` (default) uses the
        cmap's lowest colour so the mask blends into the background.
        Pass a matplotlib colour string (e.g. ``"black"``) to override.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    I_plot = np.where(I_2d > 0, I_2d, np.nan) if log else np.asarray(I_2d, dtype=float)

    if vmin is None:
        vmin = max(np.nanmin(I_plot), 1.0) if log else np.nanmin(I_plot)
    if vmax is None:
        vmax = np.nanmax(I_plot)

    # Copy the cmap so set_bad doesn't mutate matplotlib's globals.
    cmap_obj = plt.colormaps[cmap].copy()
    if masked_color is None:
        masked_color = cmap_obj(0.0)   # match the low-intensity background
    cmap_obj.set_bad(masked_color)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, layout="constrained")
    else:
        fig = ax.figure

    norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(
        I_plot,
        origin="lower",
        extent=(qxy.min(), qxy.max(), qz.min(), qz.max()),
        cmap=cmap_obj,
        norm=norm,
        interpolation="nearest",
    )

    ax.set_xlabel(r"$q_{xy}$ (Å$^{-1}$)")
    ax.set_ylabel(r"$q_z$ (Å$^{-1}$)")
    ax.set_aspect("equal")

    # Strip the spines and tick marks but keep the labels for a clean look.
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, which="both")

    if qxy_range is not None:
        ax.set_xlim(*qxy_range)
    if qz_range is not None:
        ax.set_ylim(*qz_range)
    if title:
        ax.set_title(title, fontsize=11)

    cbar = fig.colorbar(im, ax=ax, label="Intensity (counts)", pad=0.02, fraction=0.046)
    cbar.outline.set_visible(False)

    return fig, ax

def plot_1d(
    q: np.ndarray,
    intensity: np.ndarray,
    *,
    ax: "Axes | None" = None,
    label: str | None = None,
    log_y: bool = False,
) -> "tuple[Figure, Axes]":
    """Plot a 1D radial / sector profile. *Not implemented — inline in notebooks.*"""
    raise NotImplementedError


def plot_heatmap(
    q: np.ndarray,
    time: np.ndarray,
    intensity_2d: np.ndarray,
    *,
    ax: "Axes | None" = None,
    log: bool = True,
    cmap: str = "magma",
) -> "tuple[Figure, Axes]":
    """In-situ heatmap. *Not implemented — inline in notebooks.*"""
    raise NotImplementedError


def line_cut(
    cake_q: np.ndarray,
    cake_chi: np.ndarray,
    cake_intensity: np.ndarray,
    *,
    direction: str = "qz",
    width: float = 0.05,
) -> "tuple[np.ndarray, np.ndarray]":
    """Take a line cut through a cake. *Not implemented yet.*"""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Peak fit plotting (already implemented — used by 05_scherrer_analysis)
# ---------------------------------------------------------------------------


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