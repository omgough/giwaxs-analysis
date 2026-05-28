"""Peak fitting and Scherrer crystallite-size analysis.

This module fits an isolated Bragg peak on a linear background. Two
single-peak shapes are supported (Gaussian, Lorentzian), plus a
double-Lorentzian for overlapping peaks. Fitting is via
:func:`scipy.optimize.curve_fit`, which gives us standard errors on the
fitted parameters via the covariance matrix.

The Scherrer equation is implemented in q-space:

.. math::

    D_{hkl} = \\frac{2 \\pi K}{\\Delta q_{hkl}}

where :math:`\\Delta q_{hkl}` is the peak FWHM in :math:`\\text{Å}^{-1}`
and :math:`K` is the Scherrer shape factor (1 by default; common
alternatives are 0.9 for spherical and 0.94 for cubic crystallites).

Note on FWHM conventions
------------------------
The FWHM-to-shape-parameter conversion depends on which lineshape was
fit:

* Gaussian ``A·exp(-((x-μ)/σ)² / 2)``  →  FWHM = 2√(2·ln 2)·σ ≈ 2.3548·σ
* Lorentzian ``A / (π·σ·(1 + ((x-μ)/σ)²))``  →  FWHM = 2·σ

This module uses the correct formula per fitted lineshape automatically;
you don't need to do the conversion yourself.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

# scipy is imported lazily inside the fit functions so that the model
# functions, dataclasses, and the scherrer() helper remain importable
# without it.

# Conversion factor: Gaussian FWHM = GAUSSIAN_FWHM_FACTOR · σ
GAUSSIAN_FWHM_FACTOR = 2.0 * np.sqrt(2.0 * np.log(2.0))  # 2.3548...


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------


def linear_gaussian(x, m, c, A, mu, sigma):
    """Linear background plus a Gaussian peak."""
    return m * x + c + A * np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))


def linear_lorentzian(x, m, c, A, mu, sigma):
    """Linear background plus a Lorentzian peak.

    Uses the convention ``A / (π·σ·(1 + ((x-μ)/σ)²))`` to match the
    original notebooks; with this convention FWHM = 2·σ.
    """
    return m * x + c + A / (np.pi * sigma * (1.0 + ((x - mu) / sigma) ** 2))


def linear_double_lorentzian(x, m, c, A1, mu1, gamma1, A2, mu2, gamma2):
    """Linear background plus two Lorentzian peaks.

    Uses the gamma-parameterised form ``A·γ² / ((x-μ)² + γ²)``; here
    FWHM = 2·γ for each component.
    """
    lor1 = A1 * gamma1 ** 2 / ((x - mu1) ** 2 + gamma1 ** 2)
    lor2 = A2 * gamma2 ** 2 / ((x - mu2) ** 2 + gamma2 ** 2)
    return m * x + c + lor1 + lor2


# ---------------------------------------------------------------------------
# Result objects
# ---------------------------------------------------------------------------


@dataclass
class PeakFit:
    """Result of fitting a single Bragg peak on a linear background.

    Attributes
    ----------
    q_center, fwhm, amplitude
        Peak position (Å⁻¹), full-width at half-maximum (Å⁻¹), and amplitude.
    q_center_err, fwhm_err
        1-sigma standard errors from the fit covariance matrix.
    background_slope, background_intercept
        Linear background parameters ``m`` and ``c``.
    sigma
        The raw shape parameter the model was fit with (interpretation
        depends on ``model``).
    model
        Which lineshape was used: ``"gaussian"`` or ``"lorentzian"``.
    label
        Optional identifier — useful when fitting many curves in a loop.
    success
        Whether the fit converged.
    message
        Human-readable status / error message.
    """

    q_center: float
    fwhm: float
    amplitude: float
    q_center_err: float = np.nan
    fwhm_err: float = np.nan
    background_slope: float = np.nan
    background_intercept: float = np.nan
    sigma: float = np.nan
    model: Literal["gaussian", "lorentzian"] = "gaussian"
    label: str | None = None
    success: bool = True
    message: str = ""

    def as_row(self) -> dict:
        """Return a flat dict suitable for ``pandas.DataFrame``."""
        return {
            "label": self.label,
            "model": self.model,
            "q_center": self.q_center,
            "q_center_err": self.q_center_err,
            "fwhm": self.fwhm,
            "fwhm_err": self.fwhm_err,
            "amplitude": self.amplitude,
            "sigma": self.sigma,
            "background_slope": self.background_slope,
            "background_intercept": self.background_intercept,
            "success": self.success,
            "message": self.message,
        }


@dataclass
class DoublePeakFit:
    """Result of fitting two overlapping Lorentzian peaks on a linear background."""

    q_center_1: float
    fwhm_1: float
    amplitude_1: float
    q_center_2: float
    fwhm_2: float
    amplitude_2: float
    q_center_1_err: float = np.nan
    fwhm_1_err: float = np.nan
    q_center_2_err: float = np.nan
    fwhm_2_err: float = np.nan
    background_slope: float = np.nan
    background_intercept: float = np.nan
    raw_params: np.ndarray = field(default_factory=lambda: np.zeros(8))
    label: str | None = None
    success: bool = True
    message: str = ""


# ---------------------------------------------------------------------------
# Single-peak fitting
# ---------------------------------------------------------------------------


def fit_peak(
    q: np.ndarray,
    intensity: np.ndarray,
    *,
    q_range: tuple[float, float],
    model: Literal["gaussian", "lorentzian"] = "gaussian",
    label: str | None = None,
    initial_sigma: float = 0.1,
) -> PeakFit:
    """Fit a single Bragg peak inside ``q_range``.

    Parameters
    ----------
    q, intensity
        1D arrays of equal length. ``q`` typically in Å⁻¹.
    q_range
        ``(q_lo, q_hi)`` window to fit inside. Keep it tight so only
        one peak sits in the window.
    model
        ``"gaussian"`` or ``"lorentzian"``. Pseudo-Voigt is not yet
        supported — use one of the two for now and let me know if you
        need pV.
    label
        Optional identifier carried through to the result.
    initial_sigma
        Starting guess for the peak width parameter. The default of
        0.1 Å⁻¹ works for most GIWAXS peaks.

    Returns
    -------
    PeakFit
        Fit result with parameters, FWHM and 1-σ errors. If the fit
        fails, ``success=False`` and the FWHM/centre fields are NaN.
    """
    q = np.asarray(q, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if q.shape != intensity.shape:
        raise ValueError(
            f"q and intensity must have the same shape, got {q.shape} vs {intensity.shape}"
        )

    q_lo, q_hi = q_range
    mask = (q >= q_lo) & (q <= q_hi)
    if mask.sum() < 5:
        return PeakFit(
            q_center=np.nan, fwhm=np.nan, amplitude=np.nan,
            model=model, label=label, success=False,
            message=f"Fewer than 5 points inside q_range {q_range}.",
        )

    x_fit = q[mask]
    y_fit = intensity[mask]

    p0 = [
        1.0,                                    # slope m
        float(np.mean(y_fit)),                  # intercept c
        float(np.max(y_fit) - np.min(y_fit)),   # amplitude A
        float(x_fit[np.argmax(y_fit)]),         # centre mu
        float(initial_sigma),                   # sigma (>0)
    ]
    # Force sigma positive, and mu inside the fit window.
    lo = [-np.inf, -np.inf, -np.inf, q_lo, 1e-6]
    hi = [np.inf, np.inf, np.inf, q_hi, np.inf]

    func = linear_gaussian if model == "gaussian" else linear_lorentzian

    from scipy.optimize import OptimizeWarning, curve_fit  # lazy import

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            params, cov = curve_fit(func, x_fit, y_fit, p0=p0, bounds=(lo, hi))
    except (RuntimeError, ValueError) as exc:
        return PeakFit(
            q_center=np.nan, fwhm=np.nan, amplitude=np.nan,
            model=model, label=label, success=False, message=str(exc),
        )

    m, c, A, mu, sigma = params
    errs = np.sqrt(np.diag(cov)) if cov is not None else np.full(5, np.nan)
    mu_err, sigma_err = errs[3], errs[4]

    if model == "gaussian":
        fwhm = GAUSSIAN_FWHM_FACTOR * sigma
        fwhm_err = GAUSSIAN_FWHM_FACTOR * sigma_err
    else:  # lorentzian
        fwhm = 2.0 * sigma
        fwhm_err = 2.0 * sigma_err

    return PeakFit(
        q_center=float(mu),
        fwhm=float(fwhm),
        amplitude=float(A),
        q_center_err=float(mu_err),
        fwhm_err=float(fwhm_err),
        background_slope=float(m),
        background_intercept=float(c),
        sigma=float(sigma),
        model=model,
        label=label,
        success=True,
        message="ok",
    )


def fit_peaks_batch(
    q: np.ndarray,
    intensity_2d: np.ndarray,
    *,
    q_range: tuple[float, float],
    model: Literal["gaussian", "lorentzian"] = "gaussian",
    labels: list[str] | None = None,
    initial_sigma: float = 0.1,
) -> list[PeakFit]:
    """Apply :func:`fit_peak` to every column of a 2D intensity array.

    Parameters
    ----------
    q
        Shared q-axis (shape ``(N,)``).
    intensity_2d
        Intensities, shape ``(N, M)`` — one column per scan / time-point.
        This matches what you get from
        :func:`giwaxs_analysis.io.load_integrations_csv` and from
        :func:`giwaxs_analysis.integration.batch_integrate`.
    q_range, model, initial_sigma
        Passed through to :func:`fit_peak`.
    labels
        Optional list of length ``M`` to identify each column. Defaults
        to the integer index as a string.

    Returns
    -------
    list[PeakFit]
        One fit result per column. Failed fits are included with
        ``success=False`` rather than silently dropped.
    """
    q = np.asarray(q, dtype=float)
    intensity_2d = np.asarray(intensity_2d, dtype=float)

    if intensity_2d.ndim != 2:
        raise ValueError(f"intensity_2d must be 2D, got shape {intensity_2d.shape}")
    if intensity_2d.shape[0] != q.shape[0]:
        raise ValueError(
            f"intensity_2d rows ({intensity_2d.shape[0]}) must match q length ({q.shape[0]})"
        )

    n_cols = intensity_2d.shape[1]
    if labels is None:
        labels = [str(i) for i in range(n_cols)]
    elif len(labels) != n_cols:
        raise ValueError(f"Expected {n_cols} labels, got {len(labels)}")

    fits = []
    for i in range(n_cols):
        fits.append(
            fit_peak(
                q,
                intensity_2d[:, i],
                q_range=q_range,
                model=model,
                label=labels[i],
                initial_sigma=initial_sigma,
            )
        )
    return fits


# ---------------------------------------------------------------------------
# Double-Lorentzian (overlapping peaks)
# ---------------------------------------------------------------------------


def fit_double_lorentzian(
    q: np.ndarray,
    intensity: np.ndarray,
    *,
    q_range: tuple[float, float],
    initial_centres: tuple[float, float],
    initial_widths: tuple[float, float] = (0.07, 0.02),
    label: str | None = None,
) -> DoublePeakFit:
    """Fit two overlapping Lorentzian peaks on a linear background.

    For when one peak isn't enough — e.g. closely-spaced (h00) and
    (h00') reflections, or a shoulder.

    Parameters
    ----------
    q, intensity
        1D arrays of equal length.
    q_range
        ``(q_lo, q_hi)`` window to fit inside, encompassing both peaks.
    initial_centres
        ``(μ1, μ2)`` initial guesses for the peak centres — these
        matter; eyeball them from a plot first.
    initial_widths
        Initial guesses for the two HWHM values (``γ``). Defaults are
        reasonable for typical GIWAXS reflections.
    """
    q = np.asarray(q, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    q_lo, q_hi = q_range
    mask = (q >= q_lo) & (q <= q_hi)
    x_fit = q[mask]
    y_fit = intensity[mask]

    mu1, mu2 = initial_centres
    g1, g2 = initial_widths
    p0 = [
        1.0,                          # slope m
        float(np.mean(y_fit)),        # intercept c
        float(np.max(y_fit) / 2.0),   # A1
        float(mu1), float(g1),
        float(np.max(y_fit) / 2.0),   # A2
        float(mu2), float(g2),
    ]

    from scipy.optimize import OptimizeWarning, curve_fit  # lazy import

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            params, cov = curve_fit(linear_double_lorentzian, x_fit, y_fit, p0=p0)
    except (RuntimeError, ValueError) as exc:
        nan = float("nan")
        return DoublePeakFit(
            q_center_1=nan, fwhm_1=nan, amplitude_1=nan,
            q_center_2=nan, fwhm_2=nan, amplitude_2=nan,
            label=label, success=False, message=str(exc),
        )

    m, c, A1, mu1, g1, A2, mu2, g2 = params
    errs = np.sqrt(np.diag(cov)) if cov is not None else np.full(8, np.nan)

    return DoublePeakFit(
        q_center_1=float(mu1),
        fwhm_1=float(2.0 * abs(g1)),
        amplitude_1=float(A1),
        q_center_2=float(mu2),
        fwhm_2=float(2.0 * abs(g2)),
        amplitude_2=float(A2),
        q_center_1_err=float(errs[3]),
        fwhm_1_err=float(2.0 * errs[4]),
        q_center_2_err=float(errs[6]),
        fwhm_2_err=float(2.0 * errs[7]),
        background_slope=float(m),
        background_intercept=float(c),
        raw_params=np.asarray(params),
        label=label,
        success=True,
        message="ok",
    )


# ---------------------------------------------------------------------------
# Scherrer equation
# ---------------------------------------------------------------------------


def scherrer(
    fwhm_q: float | np.ndarray,
    *,
    K: float = 1.0,
    unit: Literal["A", "nm"] = "nm",
) -> float | np.ndarray:
    """Scherrer crystallite size from a peak FWHM in q.

    Implements

    .. math::

        D = \\frac{2 \\pi K}{\\Delta q}

    Parameters
    ----------
    fwhm_q
        Peak full-width at half-maximum in Å⁻¹. Scalar or array.
    K
        Scherrer shape factor. ``1.0`` matches the convention used in
        the original notebooks; ``0.9`` (spherical) and ``0.94`` (cubic)
        are also common.
    unit
        Output unit: ``"nm"`` (default) or ``"A"`` (Ångström).

    Returns
    -------
    float or ndarray
        Crystallite coherence length in the requested unit.
    """
    fwhm_q = np.asarray(fwhm_q, dtype=float)
    d_angstrom = (2.0 * np.pi * K) / fwhm_q
    if unit == "A":
        result = d_angstrom
    elif unit == "nm":
        result = d_angstrom / 10.0
    else:
        raise ValueError(f"unit must be 'A' or 'nm', got {unit!r}")
    # Preserve scalar-in / scalar-out behaviour.
    if result.ndim == 0:
        return float(result)
    return result


# ---------------------------------------------------------------------------
# Backwards-compatibility alias
# ---------------------------------------------------------------------------

# The earlier scaffold advertised this name; keep it pointing at the new function.
fit_peaks_over_time = fit_peaks_batch
