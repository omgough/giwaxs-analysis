"""File discovery and detector-frame loading.

Most beamlines emit a directory of detector frames per scan, with a
naming pattern like ``samplename_00001.edf``. The functions here wrap
:mod:`fabio` so the notebooks can ask for "all frames in this folder"
or "the i-th frame of this scan" without re-implementing globbing.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    import numpy as np


def list_frames(directory: str | Path, pattern: str = "*.edf") -> list[Path]:
    """Return the sorted list of detector-frame paths in ``directory``.

    Parameters
    ----------
    directory
        Folder to scan.
    pattern
        Glob pattern, defaulting to ``*.edf``. Use e.g. ``*.tif`` for
        Pilatus / Eiger output.
    """
    raise NotImplementedError


def load_frame(path: str | Path) -> "np.ndarray":
    """Load a single detector frame as a 2D numpy array via fabio."""
    raise NotImplementedError


def load_stack(paths: Iterable[str | Path]) -> "np.ndarray":
    """Load an iterable of frame paths into a single (N, H, W) array.

    Useful for in-situ time series where you want a stack to feed
    straight into :func:`giwaxs_analysis.integration.batch_integrate`.
    """
    raise NotImplementedError
