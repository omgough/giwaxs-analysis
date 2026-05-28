"""Peak fitting and Scherrer crystallite-size analysis.

The three ``peak_fitting_scherrer`` notebooks all do variants of:
    1. cut a window around a Bragg peak,
    2. fit a pseudo-Voigt or Gaussian to the profile,
    3. report position, FWHM, and Scherrer crystallite size.

These helpers expose that as a few named functions instead of a wall
of inline scipy calls.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import numpy as np


@dataclass
class PeakFit:
    """Result of fitting a single Bragg peak."""

    q_center: float           # peak position, Å⁻¹
    fwhm: float               # full width at half max, Å⁻¹
    amplitude: float
    background: float
    model: str                # "gaussian" | "pseudo_voigt"
    success: bool
    residual: "np.ndarray"


def fit_peak(
    q: "np.ndarray",
    intensity: "np.ndarray",
    *,
    q_range: tuple[float, float],
    model: Literal["gaussian", "pseudo_voigt"] = "pseudo_voigt",
) -> PeakFit:
    """Fit a Bragg peak inside ``q_range``.

    Parameters
    ----------
    q, intensity
        1D profile, typically from :func:`giwaxs_analysis.integration.radial_integrate`.
    q_range
        (q_lo, q_hi) window to fit inside. Keep it tight — one peak per fit.
    model
        Lineshape. Pseudo-Voigt is the standard default for diffraction.
    """
    raise NotImplementedError


def scherrer(
    fwhm_q: float,
    *,
    wavelength_A: float,
    shape_factor: float = 0.9,
) -> float:
    """Scherrer crystallite size from a peak's FWHM in q.

    Parameters
    ----------
    fwhm_q
        Peak full-width at half-maximum in Å⁻¹.
    wavelength_A
        X-ray wavelength in Å (read it off the PONI file).
    shape_factor
        Scherrer K, default 0.9 (spherical crystallites).

    Returns
    -------
    float
        Crystallite size in Å.
    """
    raise NotImplementedError


def fit_peaks_over_time(
    q: "np.ndarray",
    intensity_2d: "np.ndarray",
    *,
    q_range: tuple[float, float],
    model: Literal["gaussian", "pseudo_voigt"] = "pseudo_voigt",
) -> list[PeakFit]:
    """Apply :func:`fit_peak` row-by-row to an in-situ ``(N, npt)`` stack.

    Returns one :class:`PeakFit` per time-point.
    """
    raise NotImplementedError
