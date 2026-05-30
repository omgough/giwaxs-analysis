"""Loading and validating pyFAI calibration artefacts (PONI + mask).

The pyFAI calibration produces two files we need everywhere downstream:

* a ``.poni`` file describing the detector geometry (sample-detector
  distance, beam centre, rotation, wavelength, detector type), and
* a mask image (``.edf`` / ``.npy`` / ``.tif``) marking dead pixels and
  inter-module gaps.

The helpers here just centralise loading + sanity-checking these so the
notebooks don't each re-implement the boilerplate.

See ``docs/calibration.md`` for the GUI procedure that produces them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pyFAI.azimuthalIntegrator import AzimuthalIntegrator


@dataclass(frozen=True)
class Calibration:
    """Bundles a pyFAI integrator with its mask and source-file paths.

    Carry one of these around instead of passing PONI/mask paths through
    every function. Construct with :func:`load_calibration`.
    """

    integrator: "AzimuthalIntegrator"
    mask: np.ndarray
    poni_path: Path
    mask_path: Path


def load_calibration(poni_path: str | Path, mask_path: str | Path) -> Calibration:
    """Load a pyFAI PONI + mask pair into a :class:`Calibration`.

    Parameters
    ----------
    poni_path
        Path to the ``.poni`` file written by ``pyFAI-calib2``.
    mask_path
        Path to a mask image (``.edf``, ``.npy``, ``.tif``). Non-zero
        pixels are treated as masked (pyFAI convention).

    Returns
    -------
    Calibration
        Geometry + mask ready to use for integration.

    Raises
    ------
    FileNotFoundError
        If either file does not exist.
    ValueError
        If the mask shape doesn't match the detector defined in the PONI.
    """
    import pyFAI
    import fabio

    poni_path = Path(poni_path)
    mask_path = Path(mask_path)

    if not poni_path.is_file():
        raise FileNotFoundError(f"PONI file not found: {poni_path}")
    if not mask_path.is_file():
        raise FileNotFoundError(f"Mask file not found: {mask_path}")

    integrator = pyFAI.load(str(poni_path))

    # fabio handles .edf / .tif transparently; np.load for .npy
    if mask_path.suffix.lower() == ".npy":
        mask = np.load(mask_path)
    else:
        mask = fabio.open(str(mask_path)).data

    # pyFAI expects a boolean / 0-1 mask where True == "ignore this pixel".
    mask = mask.astype(bool)

    # Sanity-check shape against the detector. pyFAI's detector knows its
    # native shape; catching a mismatch here saves a confusing error later.
    detector_shape = integrator.detector.shape
    if detector_shape is not None and mask.shape != detector_shape:
        raise ValueError(
            f"Mask shape {mask.shape} does not match detector shape "
            f"{detector_shape}. Did you load the wrong mask for this PONI?"
        )

    return Calibration(
        integrator=integrator,
        mask=mask,
        poni_path=poni_path,
        mask_path=mask_path,
    )


def summarise(calib: Calibration) -> str:
    """Return a one-line human summary of the geometry.

    Useful as the first cell of a notebook to confirm you loaded the
    right calibration. Example output::

        Pilatus2M  dist=152.3 mm  beam=(89.4, 92.1) mm  λ=0.9763 Å  E=12.70 keV  masked=4.2%
    """
    ai = calib.integrator

    # pyFAI stores distances in metres and wavelength in metres; convert
    # to mm and Å for human-readable output.
    dist_mm = ai.dist * 1e3
    poni1_mm = ai.poni1 * 1e3
    poni2_mm = ai.poni2 * 1e3
    wavelength_A = ai.wavelength * 1e10

    # Photon energy via E = hc / λ.  hc in keV·Å ≈ 12.3984.
    energy_keV = 12.3984 / wavelength_A

    detector_name = type(ai.detector).__name__
    masked_pct = 100 * calib.mask.sum() / calib.mask.size

    return (
        f"{detector_name}  dist={dist_mm:.1f} mm  "
        f"beam=({poni1_mm:.1f}, {poni2_mm:.1f}) mm  "
        f"λ={wavelength_A:.4f} Å  E={energy_keV:.3f} keV  "
        f"masked={masked_pct:.1f}%"
    )