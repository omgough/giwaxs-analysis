# data/

Small sample data for running the example notebooks.

## Contents

- **`insitu/`** — 5 raw detector frames taken at different film thicknesses, used as a mini _in situ_ series.
- **`calibrant.edf`** — a raw image of silver behenate (AgBh) collected on a Pilatus 2M detector. AgBh is a standard calibrant with well-known d-spacings, used to fit the detector geometry.
- **`mask.edf`** — the corresponding detector mask. This is just a 2D array of 1s and 0s marking which pixels to ignore (dead pixels, beamstop shadow, module gaps, etc.).
- **`calibrant.poni`** — the pyFAI geometry file (PONI = Point Of Normal Incidence) produced from the AgBh image. It's a plain text file describing the sample-detector distance, beam centre, detector tilt and wavelength.
- **`integrations.csv`** - integrations of the example _in situ_ data.

## A note on file formats and calibrant types

Raw detector images usually come as either `.tif` or `.edf`, depending on the facility; Diamond Light Source typically writes TIFFs, ESRF typically writes EDF. Check what your beamline produces. Both are handled by `giwaxs_analysis.io`.

The example data here shows a silver behenate calibration image, but this may differ for your experiment. LaB6 is also common for GIWAXS.

## Generating your own

The `.poni` and mask are tied to a specific beamtime — you'll need to make your own for your data. See [`docs/calibration.md`](../docs/calibration.md) for the full step-by-step using `pyFAI-calib2`. The `pyFAI-calib2` documentation can be found [here](https://pyfai.readthedocs.io/en/latest/man/pyFAI-calib2.html), and it is highly recommended to use their GUI for generating the `.poni` file.

## Acknowledgements

This demo data used resources of the Advanced Light Source, which is a
DOE Office of Science User Facility under contract no. DE-AC02-05CH11231.
O. Gough was supported in part by an ALS Doctoral Fellowship in Residence.

The demo data included in this repository was collected at the Advanced
Light Source, Lawrence Berkeley National Laboratory.