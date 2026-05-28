# giwaxs-analysis

Python tools for processing **Grazing-Incidence Wide-Angle X-ray Scattering (GIWAXS)** data — from raw detector frames through to 1D integrations, in-situ heatmaps, peak fitting and Scherrer crystallite-size analysis. Built on top of [pyFAI](https://pyfai.readthedocs.io/) and [fabio](https://fabio.readthedocs.io/).

This package grew out of a set of notebooks written during my PhD on perovskite/organic thin-film characterisation. The hope is to make the same workflows usable for other groups working with GIWAXS data from synchrotron and lab sources.

## What's in here

| Module | What it does |
| --- | --- |
| `giwaxs_analysis.calibration` | Load and validate pyFAI PONI files + detector masks |
| `giwaxs_analysis.io` | Discover and load detector frames (EDF, TIFF, ...) |
| `giwaxs_analysis.integration` | 1D radial, sector and cake integrations; batch-integrate in-situ stacks |
| `giwaxs_analysis.plotting` | Detector images, q-space maps, line cuts, in-situ heatmaps |
| `giwaxs_analysis.peak_fitting` | Peak fitting (Gaussian / pseudo-Voigt) and Scherrer analysis |

End-to-end examples live in [`notebooks/`](notebooks/).

## Installation

The cleanest path on Linux/macOS is via conda, because pyFAI's wheels can be fiddly:

```bash
conda env create -f environment.yml
conda activate giwaxs
pip install -e .
```

If you'd rather stay in pure-pip land:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Optional extras:

```bash
pip install -e ".[plotnine]"   # for the plotnine-based plotting helpers
pip install -e ".[fityk]"      # for the fityk peak-fitting backend
pip install -e ".[dev]"        # for tests / linting / jupyter
```

## Prerequisite: calibrate your geometry

Before any of the analysis tools will produce meaningful q-axes, you need a `.poni` file (geometry) and a mask image. Generate these once per beamtime with `pyFAI-calib2`. The full GUI procedure is written up in [`docs/calibration.md`](docs/calibration.md).

## Quickstart

```python
from giwaxs_analysis import calibration, io, integration, plotting

# 1. load the geometry produced by pyFAI-calib2
calib = calibration.load_calibration("calib.poni", "mask.edf")

# 2. load a detector frame
frame = io.load_frame("data/sample_00001.edf")

# 3. integrate
q, I = integration.radial_integrate(frame, calib)

# 4. plot
plotting.plot_1d(q, I, log_y=True)
```

For an in-situ time series:

```python
paths = io.list_frames("data/insitu_scan/")
stack = io.load_stack(paths)
q, I_2d = integration.batch_integrate(stack, calib, mode="radial")
plotting.plot_heatmap(q, time=range(len(stack)), intensity_2d=I_2d)
```

## Example notebooks

See [`notebooks/`](notebooks/) for runnable examples covering each workflow. Each notebook is self-contained and uses the helpers from this package — they're meant as both demos and copy-paste starting points.

## Status

Early alpha. The API surface is small and likely to change. Issues and pull requests welcome — particularly from anyone working with different detector geometries.

## Citation

If you use this in published work, please cite it via the metadata in [`CITATION.cff`](CITATION.cff). GitHub will surface a "Cite this repository" button on the repo page.

## License

[MIT](LICENSE)
