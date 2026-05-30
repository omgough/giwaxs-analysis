# giwaxs-analysis

Python tools for processing **Grazing-Incidence Wide-Angle X-ray Scattering (GIWAXS)** data - from raw detector frames through to 1D / wedge integrations, general GIWAXS data visualisation, _in situ_ heatmaps / waterfall plots, peak fitting and Scherrer crystallite-size analysis. Built on top of [pyFAI](https://pyfai.readthedocs.io/) and [fabio](https://fabio.readthedocs.io/).

This package is adapted from a set of notebooks written during my DPhil, focusing on _in situ_ GIWAXS during thin film deposition using the [MINERVA chamber](https://pubs.aip.org/aip/rsi/article/88/10/103901/834213/MINERVA-A-facility-to-study-Microstructure-and). The hope is to make the same workflows usable for other groups working with GIWAXS data from synchrotron, particularly _in situ_ GIXD data, be it from MINERVA or otherwise. The examples given here are _in situ_ GIWAXS during thin film deposition via vacuum thermal evaporation (VTE), but can also be used for Grazing-Incidence Small-Angle X-ray Scattering (GISAXS) and any kind of _in situ_ series (_e.g._ annealing, different deposition methods, _etc._)


## Repository layout

```
giwaxs-analysis/
├── src/giwaxs_analysis/        # the package itself
│   ├── __init__.py
│   ├── calibration.py          # load / validate PONI files and masks
│   ├── io.py                   # discover and load detector frames
│   ├── integration.py          # 1D radial / sector / cake integrations
│   ├── plotting.py             # detector images, line cuts, heatmaps
│   └── peak_fitting.py         # Gaussian / pseudo-Voigt fits, Scherrer
├── notebooks/                  # runnable end-to-end examples
│   └── 06_peak_fitting_scherrer.ipynb
├── data/                       # small sample data for the notebooks
│   ├── insitu/                 # mini in-situ series (thickness scan)
│   └── README.md
├── docs/
│   └── calibration.md          # how to generate a .poni with pyFAI-calib2
├── tests/
│   └── test_imports.py
├── pyproject.toml              # package metadata + dependencies
├── environment.yml             # conda environment (recommended install)
├── CITATION.cff                # how to cite this work
├── LICENSE                     # MIT
└── README.md
```


## What's in here

| Module | What it does |
| --- | --- |
| `giwaxs_analysis.calibration` | Load and validate pyFAI PONI files + detector masks |
| `giwaxs_analysis.io` | Discover and load detector frames (EDF, TIFF, ...) |
| `giwaxs_analysis.integration` | 1D azimuthal (full + sector) and 2D cake integrations; batch-integrate _in situ_ stacks |
| `giwaxs_analysis.plotting` | Detector images, q-space maps, line cuts, _in situ_ heatmaps |
| `giwaxs_analysis.peak_fitting` | Peak fitting (Gaussian / Lorentzian) and Scherrer analysis |

End-to-end examples live in [`notebooks/`](notebooks/).


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

For an _in situ_ time series:

```python
paths = io.list_frames("data/insitu_scan/")
stack = io.load_stack(paths)
q, I_2d = integration.batch_integrate(stack, calib, mode="radial")
plotting.plot_heatmap(q, time=range(len(stack)), intensity_2d=I_2d)
```

## Example notebooks

See [`notebooks/`](notebooks/) for runnable examples covering each workflow. Each notebook is self-contained and uses the helpers from this package — they're meant as both demos and copy-paste starting points.


## Citation

If you use this package in published work, please cite **both**:

1. The associated thesis: O. Gough, *In situ microstructural characterisation of organic solar cells*, DPhil thesis, University of Oxford, 2025.
2. The pyFAI library this package is built on — see [pyFAI's citation page](https://pyfai.readthedocs.io/en/stable/publications.html) (the JAC 2015 paper is the standard reference for azimuthal integration).

Machine-readable metadata for both lives in [`CITATION.cff`](CITATION.cff); GitHub will surface a "Cite this repository" button on the repo page.

## License

[MIT](LICENSE)
