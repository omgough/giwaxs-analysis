"""1D and 2D integrations of GIWAXS detector frames.

The three workflows in your existing notebooks all reduce to a small
set of calls into :mod:`pyFAI`. The functions here name those clearly
so a notebook reads as ``radial_integrate(...)`` instead of an opaque
``ai.integrate1d(...)``.

Conventions
-----------
* ``q`` is in inverse Angstroms (Å⁻¹) by default — matching what
  ``pyFAI`` produces with ``unit="q_A^-1"``.
* Intensity is returned as a 1D or 2D numpy array; q (and chi where
  relevant) are returned as separate axes so they're easy to drop
  straight into a DataFrame.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import numpy as np

    from .calibration import Calibration


def radial_integrate(
    frame: "np.ndarray",
    calib: "Calibration",
    *,
    npt: int = 1000,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
) -> tuple["np.ndarray", "np.ndarray"]:
    """Full-azimuthal radial integration → (q, I)."""
    raise NotImplementedError


def sector_integrate(
    frame: "np.ndarray",
    calib: "Calibration",
    *,
    chi_range: tuple[float, float],
    npt: int = 1000,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
) -> tuple["np.ndarray", "np.ndarray"]:
    """Sector (azimuthal-slice) integration → (q, I).

    ``chi_range`` is in degrees. Use this for in-plane vs out-of-plane
    GIWAXS comparisons.
    """
    raise NotImplementedError


def cake(
    frame: "np.ndarray",
    calib: "Calibration",
    *,
    npt_rad: int = 1000,
    npt_azim: int = 360,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
) -> tuple["np.ndarray", "np.ndarray", "np.ndarray"]:
    """Reshape detector frame into a (q, χ) cake → (q, chi, I_2d)."""
    raise NotImplementedError


def batch_integrate(
    stack: "np.ndarray",
    calib: "Calibration",
    *,
    mode: Literal["radial", "sector"] = "radial",
    **kwargs,
) -> tuple["np.ndarray", "np.ndarray"]:
    """Integrate a (N, H, W) stack → (q, I_2d) of shape (N, npt).

    For in-situ time series — the heatmap notebooks should call this
    instead of looping over frames manually.
    """
    raise NotImplementedError
