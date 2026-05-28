"""File discovery, frame loading, and tabular-CSV loading.

Two flavours of I/O live here:

* **Detector frames** — wrappers around :mod:`fabio` for loading a
  single frame or a directory of frames into a stack.
* **Integrations CSV** — the common pattern in your notebooks of a
  CSV where column 0 is the q-axis and remaining columns are
  integrated intensities (one column per sample / scan / time-point).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Iterable

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass  # fabio is imported lazily inside the functions that use it


# ---------------------------------------------------------------------------
# Detector frame loading (stubs — fill in when porting integration notebooks)
# ---------------------------------------------------------------------------


def list_frames(directory: str | Path, pattern: str = "*.edf") -> list[Path]:
    """Return the sorted list of detector-frame paths in ``directory``."""
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"Not a directory: {directory}")
    return sorted(directory.glob(pattern))


def load_frame(path: str | Path) -> np.ndarray:
    """Load a single detector frame as a 2D numpy array via fabio."""
    import fabio  # lazy import — fabio is heavy and not needed for the CSV workflow

    return fabio.open(str(path)).data


def load_stack(paths: Iterable[str | Path]) -> np.ndarray:
    """Load an iterable of frame paths into a single ``(N, H, W)`` array."""
    frames = [load_frame(p) for p in paths]
    if not frames:
        raise ValueError("No frame paths provided.")
    return np.stack(frames, axis=0)


# ---------------------------------------------------------------------------
# Integrations CSV
# ---------------------------------------------------------------------------


@dataclass
class IntegrationsCSV:
    """Container for a loaded integrations CSV.

    Attributes
    ----------
    q
        1D q-axis (Å⁻¹), shape ``(N,)``.
    intensities
        2D intensities, shape ``(N, M)``, one column per dataset.
    labels
        Length-``M`` list of column labels (taken from the CSV header).
    source
        Source path, for traceability.
    """

    q: np.ndarray
    intensities: np.ndarray
    labels: list[str]
    source: Path


def load_integrations_csv(
    path: str | Path,
    *,
    skip_data_rows: int = 2,
    q_column: int = 0,
) -> IntegrationsCSV:
    """Load a CSV of pre-computed 1D integrations.

    Expected layout
    ---------------
    * Row 0 is the column header.
    * The next ``skip_data_rows`` rows are metadata / placeholders and
      are skipped. This matches the format the existing notebooks
      produced (which used ``df.iloc[2:, ...]`` to skip two rows).
    * Column ``q_column`` (default 0) is the q-axis in Å⁻¹.
    * All other columns are intensities, one per scan / sample.

    If your CSVs have no metadata rows, pass ``skip_data_rows=0``.

    Parameters
    ----------
    path
        Path to the CSV file.
    skip_data_rows
        Number of rows to skip *after* the header.
    q_column
        Index of the column containing the q-axis.

    Returns
    -------
    IntegrationsCSV
        ``q``, ``intensities``, ``labels``, ``source`` bundled together.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"CSV not found: {path}")

    df = pd.read_csv(path)
    if skip_data_rows:
        df = df.iloc[skip_data_rows:].reset_index(drop=True)

    q = df.iloc[:, q_column].to_numpy(dtype=float)

    intensity_cols = [c for i, c in enumerate(df.columns) if i != q_column]
    intensities = df[intensity_cols].to_numpy(dtype=float)
    labels = [str(c) for c in intensity_cols]

    return IntegrationsCSV(q=q, intensities=intensities, labels=labels, source=path)
