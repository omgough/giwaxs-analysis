"""1D and 2D integrations of GIWAXS detector frames.

The three workflows in your existing notebooks all reduce to a small
set of calls into :mod:`pyFAI`. The functions here name those clearly
so a notebook reads as ``radial_integrate(...)`` instead of an opaque
``ai.integrate1d(...)``.

Conventions
-----------
* ``q`` is in inverse Angstroms (Å⁻¹) by default — matching what
  ``pyFAI`` produces with ``unit="q_A^-1"``.
* ``chi`` (azimuthal angle) is in degrees, between –180° and 180°.
* Intensity is returned as a 1D or 2D numpy array; q (and chi where
  relevant) are returned as separate axes so they're easy to drop
  straight into a DataFrame.
* The mask in ``calib.mask`` is applied automatically — you don't need
  to pass it through every call.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from .calibration import Calibration


def radial_integrate(
    frame: np.ndarray,
    calib: "Calibration",
    *,
    npt: int = 1000,
    q_range: tuple[float, float] | None = None,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
    polarization_factor: float | None = 0.99,
    correct_solid_angle: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Full-azimuthal radial integration → ``(q, I)``.

    Integrates the detector frame over the full azimuthal angle χ at
    each q bin, returning a 1D intensity profile I(q).

    Parameters
    ----------
    frame
        2D detector image, shape ``(H, W)``.
    calib
        Calibration bundle from :func:`giwaxs_analysis.calibration.load_calibration`.
    npt
        Number of radial bins.
    q_range
        Optional ``(q_min, q_max)`` to limit the radial range. In units of ``unit``.
    unit
        pyFAI radial unit. Default is q in Å⁻¹.
    polarization_factor
        Polarization correction factor (between –1 and 1). 0.99 is
        typical for horizontally-polarized synchrotron sources. Pass
        ``None`` to disable the correction (e.g. for unpolarized lab sources).
    correct_solid_angle
        Apply the per-pixel solid-angle correction.

    Returns
    -------
    (q, I) : tuple of np.ndarray
        Two 1D arrays of length ``npt``.
    """
    q, I = calib.integrator.integrate1d(
        frame,
        npt=npt,
        mask=calib.mask,
        radial_range=q_range,
        unit=unit,
        polarization_factor=polarization_factor,
        correctSolidAngle=correct_solid_angle,
        method="bbox",
    )
    return q, I


def sector_integrate(
    frame: np.ndarray,
    calib: "Calibration",
    *,
    chi_range: tuple[float, float],
    npt: int = 1000,
    q_range: tuple[float, float] | None = None,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
    polarization_factor: float | None = 0.99,
    correct_solid_angle: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Sector (azimuthal-slice) integration → ``(q, I)``.

    Integrates over a wedge of χ instead of the full circle. Use this
    for in-plane vs out-of-plane GIWAXS comparisons.

    Parameters
    ----------
    chi_range
        Azimuthal range ``(chi_min, chi_max)`` in **degrees**, between
        –180° and 180°. Conventional GIWAXS wedges (depend on detector
        orientation): in-plane around χ ≈ 0° or 180°, out-of-plane
        around χ ≈ ±90°.

    See :func:`radial_integrate` for the remaining parameters.

    Returns
    -------
    (q, I) : tuple of np.ndarray
        Two 1D arrays of length ``npt``.
    """
    q, I = calib.integrator.integrate1d(
        frame,
        npt=npt,
        mask=calib.mask,
        azimuth_range=chi_range,
        radial_range=q_range,
        unit=unit,
        polarization_factor=polarization_factor,
        correctSolidAngle=correct_solid_angle,
        method="bbox",
    )
    return q, I


def cake(
    frame: np.ndarray,
    calib: "Calibration",
    *,
    npt_rad: int = 1000,
    npt_azim: int = 360,
    q_range: tuple[float, float] | None = None,
    chi_range: tuple[float, float] | None = None,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
    polarization_factor: float | None = 0.99,
    correct_solid_angle: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reshape detector frame into a (q, χ) cake → ``(q, chi, I_2d)``.

    The cake view is the strongest visual check on calibration: rings
    in the raw image become straight vertical lines if the geometry
    is right; curved or tilted lines mean something's off.

    Parameters
    ----------
    npt_rad, npt_azim
        Number of radial and azimuthal bins.
    chi_range
        Optional ``(chi_min, chi_max)`` to limit the azimuthal range.
        Defaults to pyFAI's full range.

    See :func:`radial_integrate` for the remaining parameters.

    Returns
    -------
    (q, chi, I_2d) : tuple of np.ndarray
        ``q`` shape ``(npt_rad,)``, ``chi`` shape ``(npt_azim,)``,
        ``I_2d`` shape ``(npt_azim, npt_rad)``.
    """
    I_2d, q, chi = calib.integrator.integrate2d(
        frame,
        npt_rad=npt_rad,
        npt_azim=npt_azim,
        mask=calib.mask,
        radial_range=q_range,
        azimuth_range=chi_range,
        unit=unit,
        polarization_factor=polarization_factor,
        correctSolidAngle=correct_solid_angle,
        method="bbox",
    )
    return q, chi, I_2d


def batch_integrate(
    stack: np.ndarray,
    calib: "Calibration",
    *,
    mode: Literal["radial", "sector"] = "radial",
    **kwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate an ``(N, H, W)`` stack → ``(q, I_2d)`` of shape ``(N, npt)``.

    For in-situ time series. Loops over the stack and calls the
    per-frame integrator. ``kwargs`` are forwarded.

    Parameters
    ----------
    stack
        3D array of detector frames, shape ``(N, H, W)``.
    calib
        Calibration bundle.
    mode
        ``"radial"`` for full-azimuthal integration (:func:`radial_integrate`)
        or ``"sector"`` for a χ-wedge (:func:`sector_integrate`, requires
        ``chi_range`` in ``kwargs``).
    **kwargs
        Forwarded to the per-frame integrator: ``npt``, ``q_range``,
        ``unit``, ``polarization_factor``, ``correct_solid_angle``, and
        (for ``mode="sector"``) ``chi_range``.

    Returns
    -------
    (q, I_2d) : tuple of np.ndarray
        ``q`` shape ``(npt,)``, ``I_2d`` shape ``(N, npt)``.
    """
    if stack.ndim != 3:
        raise ValueError(f"stack must be 3D (N, H, W); got shape {stack.shape}.")

    if mode == "radial":
        integrate = radial_integrate
    elif mode == "sector":
        integrate = sector_integrate
    else:
        raise ValueError(f"Unknown mode {mode!r}; expected 'radial' or 'sector'.")

    q_ref: np.ndarray | None = None
    intensities: list[np.ndarray] = []
    for frame in stack:
        q, I = integrate(frame, calib, **kwargs)
        if q_ref is None:
            q_ref = q
        intensities.append(I)

    return q_ref, np.stack(intensities, axis=0)