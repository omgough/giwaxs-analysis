"""1D and 2D integrations of GIWAXS detector frames.

The three workflows in your existing notebooks all reduce to a small
set of calls into :mod:`pyFAI`. The functions here name those clearly
so a notebook reads as ``azimuthal_integrate(...)`` instead of an
opaque ``ai.integrate1d(...)``.

Naming follows GIWAXS convention:

* **azimuthal integration** — integrate over χ at each q → I(q).
* **radial integration** — integrate over a narrow q-range → I(χ).
  (To be added when needed.)
* **cake** — the full 2D (q, χ) unwrapping.

Conventions
-----------
* ``q`` is in inverse Angstroms (Å⁻¹) by default — matching what
  ``pyFAI`` produces with ``unit="q_A^-1"``.
* ``chi`` (azimuthal angle) is in degrees, between –180° and 180°.
* The mask in ``calib.mask`` is applied automatically.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from .calibration import Calibration


def azimuthal_integrate(
    frame: np.ndarray,
    calib: "Calibration",
    *,
    chi_range: tuple[float, float] | None = None,
    q_range: tuple[float, float] | None = None,
    npt: int = 1000,
    unit: Literal["q_A^-1", "q_nm^-1", "2th_deg"] = "q_A^-1",
    polarization_factor: float | None = 0.99,
    correct_solid_angle: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate over χ at each q → ``(q, I)``.

    Pass ``chi_range`` to restrict the azimuthal wedge (e.g. for
    in-plane vs out-of-plane comparisons); leave it as ``None`` for a
    full-azimuthal integration.

    Parameters
    ----------
    frame
        2D detector image, shape ``(H, W)``.
    calib
        Calibration bundle from :func:`giwaxs_analysis.calibration.load_calibration`.
    chi_range
        Azimuthal range ``(chi_min, chi_max)`` in **degrees**, between
        –180° and 180°. ``None`` = full circle.
    q_range
        Optional ``(q_min, q_max)`` to limit the radial range. In units of ``unit``.
    npt
        Number of radial bins.
    unit
        pyFAI radial unit. Default is q in Å⁻¹.
    polarization_factor
        Polarization correction factor (–1 to 1). 0.99 is typical for
        horizontally-polarized synchrotron sources. Pass ``None`` to disable.
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

    See :func:`azimuthal_integrate` for the remaining parameters.

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
    **kwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """Azimuthally integrate an ``(N, H, W)`` stack → ``(q, I_2d)``.

    For in-situ time series. Loops over the stack and calls
    :func:`azimuthal_integrate` per frame. ``kwargs`` are forwarded.

    Parameters
    ----------
    stack
        3D array of detector frames, shape ``(N, H, W)``.
    calib
        Calibration bundle.
    **kwargs
        Forwarded to :func:`azimuthal_integrate`: ``chi_range``,
        ``q_range``, ``npt``, ``unit``, ``polarization_factor``,
        ``correct_solid_angle``.

    Returns
    -------
    (q, I_2d) : tuple of np.ndarray
        ``q`` shape ``(npt,)``, ``I_2d`` shape ``(N, npt)``.
    """
    if stack.ndim != 3:
        raise ValueError(f"stack must be 3D (N, H, W); got shape {stack.shape}.")

    q_ref: np.ndarray | None = None
    intensities: list[np.ndarray] = []
    for frame in stack:
        q, I = azimuthal_integrate(frame, calib, **kwargs)
        if q_ref is None:
            q_ref = q
        intensities.append(I)

    return q_ref, np.stack(intensities, axis=0)