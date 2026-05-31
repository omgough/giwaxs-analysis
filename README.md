# giwaxs-analysis

Python tools for processing **Grazing-Incidence Wide-Angle X-ray Scattering (GIWAXS)** data — from raw detector frames through to 1D / wedge integrations, general GIWAXS data visualisation, *in situ* heatmaps / waterfall plots, peak fitting and Scherrer crystallite-size analysis. Built on top of [pyFAI](https://pyfai.readthedocs.io/) and [fabio](https://fabio.readthedocs.io/).

This package is adapted from a set of notebooks written during my DPhil, focusing on *in situ* GIWAXS during thin film deposition using the [MINERVA chamber](https://pubs.aip.org/aip/rsi/article/88/10/103901/834213/MINERVA-A-facility-to-study-Microstructure-and). The hope is to make the same workflows usable for other groups working with GIWAXS data from synchrotron, particularly *in situ* GIXD data, be it from MINERVA or otherwise. The examples given here are *in situ* GIWAXS during thin film deposition via vacuum thermal evaporation (VTE), but can also be used for Grazing-Incidence Small-Angle X-ray Scattering (GISAXS) and any kind of *in situ* series (*e.g.* annealing, different deposition methods, *etc.*).

## Contents

- [New to GIWAXS? Start here](#new-to-giwaxs-start-here)
- [Repository layout](#repository-layout)
- [What's in here](#whats-in-here)
- [Installation](#installation)
- [What you'll need before you start](#what-youll-need-before-you-start)
- [Prerequisite: calibrate your geometry](#prerequisite-calibrate-your-geometry)
- [Quickstart](#quickstart)
- [Example notebooks](#example-notebooks)
- [Troubleshooting / FAQ](#troubleshooting--faq)
- [Contributing & getting help](#contributing--getting-help)
- [Citation](#citation)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## New to GIWAXS? Start here

<details>
<summary>30-second primer on what this package expects from you (click to expand)</summary>

GIWAXS measurements produce 2D detector images of the X-ray scattering pattern from a thin film. To turn those into something physically meaningful, you need three things:

- **Detector frames** — the raw 2D images from the beamline (EDF, TIFF, HDF5, CBF, etc. — anything [fabio](https://fabio.readthedocs.io/) can read).
- **A geometry file (`.poni`)** — a small text file produced by `pyFAI-calib2` that describes the experimental geometry: sample-to-detector distance, beam centre, detector tilt, wavelength. Without this the package can't map detector pixels to scattering vector *q*.
- **A detector mask** — a 2D image marking dead pixels, beamstops, gaps and any region you want excluded from integration.

Once you have those, this package handles the rest: integration to 1D *I(q)* curves or 2D cake / q-space maps, batch processing of *in situ* time series, plotting, peak fitting and _in situ_ trends.

If you've never made a PONI file before, see [`docs/calibration.md`](docs/calibration.md) for a step-by-step walkthrough with `pyFAI-calib2`.

</details>

## Repository layout

Feel free to fork this repository and edit the modules as you choose! Module names self explanatory; the most likely module you may want to change is `plotting.py` to plot the data in the style you prefer. 

```
giwaxs-analysis/
├── src/giwaxs_analysis/        # the package itself
│   ├── __init__.py
│   ├── calibration.py          # load / validate PONI files and masks
│   ├── io.py                   # discover and load detector frames
│   ├── integration.py          # 1D radial / sector / cake integrations
│   ├── plotting.py             # detector images, line cuts, heatmaps
│   └── peak_fitting.py         # Gaussian / Lorentzian fits, Scherrer
├── notebooks/                  # runnable end-to-end examples
│   ├── 01_check_calibration.ipynb
│   ├── 02_plot_raw_frames.ipynb
│   ├── 03_reciprocal_space_map.ipynb
│   ├── 04_azimuthal_integration.ipynb
│   └── 05_scherrer_analysis.ipynb
├── data/                       # small sample data for the notebooks
│   ├── insitu/                 # mini in situ series (thickness scan)
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

| Module                         | What it does                                                                            |
| ------------------------------ | --------------------------------------------------------------------------------------- |
| `giwaxs_analysis.calibration`  | Load and validate pyFAI PONI files + detector masks                                     |
| `giwaxs_analysis.io`           | Discover and load detector frames (EDF, TIFF, ...)                                      |
| `giwaxs_analysis.integration`  | 1D azimuthal (full + sector) and 2D cake integrations; batch-integrate *in situ* stacks |
| `giwaxs_analysis.plotting`     | Detector images, q-space maps, line cuts, *in situ* heatmaps                            |
| `giwaxs_analysis.peak_fitting` | Peak fitting (Gaussian / Lorentzian) and Scherrer analysis                              |

End-to-end examples live in [`notebooks/`](notebooks).

## Installation

Requires **Python ≥ 3.10**. Conda is recommended because `pyFAI` and `fabio` depend on system libraries (HDF5, libtiff) that are far easier to install via `conda-forge` than via PyPI wheels — especially on Linux and Apple Silicon.

### Option 1 — conda (recommended)

```bash
git clone https://github.com/omgough/giwaxs-analysis.git
cd giwaxs-analysis
conda env create -f environment.yml
conda activate giwaxs
```

`environment.yml` already runs `pip install -e .`, so the `giwaxs_analysis` package is importable from any directory once the env is active.

### Option 2 — pip / venv

```bash
git clone https://github.com/omgough/giwaxs-analysis.git
cd giwaxs-analysis
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

If `pyFAI` or `fabio` fail to build on pip, fall back to the conda route above.

### Optional extras

```bash
pip install -e ".[plotnine]"   # extra plotting backend used in some notebooks
pip install -e ".[fityk]"      # optional peak-fitting backend
pip install -e ".[dev]"        # pytest, ruff, jupyterlab, nbstripout
pip install -e ".[all]"        # everything above
```

### Verify the install

```bash
python -c "import giwaxs_analysis; print(giwaxs_analysis.__name__, 'ok')"
pytest -q
```

## What you'll need before you start

To process your own data, gather these three things first:

1. **Detector frames** — any format [fabio](https://fabio.readthedocs.io/) can read (EDF, TIFF, HDF5, CBF, etc.). For an *in situ* series, all frames in one folder with a consistent naming pattern (e.g. `scan_00001.edf`, `scan_00002.edf`, …) is easiest.
2. **A `.poni` file** — your experimental geometry, produced by `pyFAI-calib2` from a calibration standard (LaB₆, CeO₂, AgBeh). One per beamtime / geometry change. See [`docs/calibration.md`](docs/calibration.md).
3. **A mask file** — a 2D image (EDF / TIFF / NumPy `.npy`) marking pixels to exclude: beamstop, dead modules, detector gaps, sample shadow. Usually drawn in `pyFAI-calib2` or `pyFAI-drawmask`.

Drop them somewhere sensible (e.g. a `calibration/` subfolder alongside your data) and you're ready for the [Quickstart](#quickstart).

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

For an *in situ* time series:

```python
paths = io.list_frames("data/insitu_scan/")
stack = io.load_stack(paths)
q, I_2d = integration.batch_integrate(stack, calib, mode="radial")
plotting.plot_heatmap(q, time=range(len(stack)), intensity_2d=I_2d)
```

## Example notebooks

See [`notebooks/`](notebooks) for runnable examples covering each workflow. Each notebook is self-contained and uses the helpers from this package — they're meant as both demos and copy-paste starting points.

The bundled sample data in [`data/insitu/`](data) is enough to run every notebook end-to-end without needing your own beamtime data first. Good for kicking the tyres.


## Contributing & getting help

Issues, questions and pull requests are all welcome.

- **Asking a usage question?** Open a [discussion](https://github.com/omgough/giwaxs-analysis/discussions) (or an issue tagged `question`). If you can share a small frame + PONI + mask, it makes diagnosis much easier.
- **Submitting code?** Fork, branch, run `pip install -e ".[dev]"`, make sure `pytest` and `ruff check .` both pass, then open a PR. Keep changes focused; one PR per concern.

For the calibration step itself (PONI generation), the [pyFAI documentation](https://pyfai.readthedocs.io/) and `pyFAI-calib2` are the canonical references — questions about calibration physics are best asked there.

## Citation

If you use this package in published work, please cite **both**:

1. The associated thesis: O. Gough, *In situ microstructural characterisation of organic solar cells*, DPhil thesis, University of Oxford, 2025.
2. The pyFAI library this package is built on — see [pyFAI's citation page](https://pyfai.readthedocs.io/en/stable/publications.html) (the JAC 2015 paper is the standard reference for azimuthal integration).

Machine-readable metadata for both lives in [`CITATION.cff`](CITATION.cff); GitHub will surface a "Cite this repository" button on the repo page.

## License

[MIT](LICENSE)

## Acknowledgements

This work was carried out while O. Gough was a [Wolfson-Marriott Postgraduate Scholar](https://www.wolfson.ox.ac.uk/news/wolfson-dphil-and-team-awarded-block-allocation-group-status-for-synchrotron-experiment/) at Wolfson College, Oxford, and a [UKRI studentship](https://gtr.ukri.org/projects?ref=studentship-2606412) recipient in the [Advanced Functional Materials and Devices (AFMD) group](https://www.physics.ox.ac.uk/research/group/advanced-functional-materials-and-devices-afmd-group) under the supervision of Prof. Moritz Riede.

The demo data included in this repository was collected at the Advanced Light Source, Lawrence Berkeley National Laboratory whilst on a Doctoral Fellowship in the [Su Materials Lab](https://sumaterialslab.lbl.gov/) under the supervision of Dr. Gregory Su. This research used resources of the Advanced Light Source, which is a DOE Office of Science User Facility under contract no. DE-AC02-05CH11231. O. Gough was supported in part by an ALS Doctoral Fellowship in Residence.
