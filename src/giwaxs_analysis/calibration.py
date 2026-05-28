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

if TYPE_CHECKING:
    import numpy as np
    from pyFAI.azimuthalIntegrator import AzimuthalIntegrator


@dataclass(frozen=True)
class Calibration:
    """Bundles a pyFAI integrator with its mask and source-file paths.

    Carry one of these around instead of passing PONI/mask paths through
    every function. Construct with :func:`load_calibration`.
    """

    integrator: "AzimuthalIntegrator"
    mask: "np.ndarray"
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
        pixels are treated as masked.

    Returns
    -------
    Calibration
        Geometry + mask ready to use for integration.

    Raises
    ------
    FileNotFoundError
        If either file does not exist.
    """
    raise NotImplementedError("Port the PONI/mask loading from the existing notebooks.")


def summarise(calib: Calibration) -> str:
    """Return a one-line human summary of the geometry (distance, beam centre, λ).

    Useful as the first cell of a notebook to confirm you loaded the
    right calibration.
    """
    raise NotImplementedError
