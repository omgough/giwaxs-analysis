# Calibrating the geometry with pyFAI-calib2

Before any of the integration / plotting tools in this package will give you a meaningful q-axis, you need to calibrate the experimental geometry. This produces two files used everywhere downstream:

- a `.poni` file describing the detector geometry (sample-detector distance, beam centre, rotations, wavelength, detector type), and
- a **mask image** marking dead pixels and the gaps between detector modules.

Both are loaded together by `giwaxs_analysis.calibration.load_calibration`.

## Prerequisites

- pyFAI installed and on PATH (it ships the `pyFAI-calib2` GUI):

  ```bash
  conda install -c conda-forge pyFAI
  # or
  pip install pyFAI
  ```

- A **calibrant image** — a detector frame from a known reference sample (e.g. LaB₆, CeO₂, silver behenate). Most beamlines collect one at the start of every beamtime.
- Beam energy (or wavelength), detector type, and approximate sample-detector distance from the logbook.

## Step-by-step

Open the GUI:

```bash
pyFAI-calib2
```

### 1. Experiment settings

Fill in:

- **Beam energy** (or wavelength) — read from the beamline logbook.
- **Calibrant** — pick from the dropdown (LaB6_SRM660a, CeO2, AgBh, ...).
- **Detector** — pick the model you used, or define one if it's missing.
- **Calibration image** — point to the detector frame of the calibrant.

Click **Next**.

### 2. Mask

Mask dead pixels and inter-module gaps so they don't get integrated:

- Use the pen tool to paint over bad pixels by hand, **or**
- Threshold below (and above) a given intensity to catch dead/hot pixels in bulk.
- For multi-module detectors (Pilatus, Eiger), mask the gaps between modules.

**Save the mask file** via the save button in the mask sub-tab — you'll need it for every downstream script. Click **Next**.

### 3. Peak picking

Starting from the innermost ring and working outward:

- Pick **~10 peaks** across each ring.
- Make sure they sit exactly on the ring — erase any stray points.
- Add more points to a ring if the fit looks loose.

**Save the picked rings.** This is important if you ever need to come back and re-calibrate without redoing the picking. Click **Next**.

### 4. Geometry fitting

- The distance in the geometry sub-tab should match the sample-detector distance from the logbook.
- All three rotation angles should be very close to zero (the detector is mounted nearly flat).
- Set **Rotation 1, 2, 3 → 0°** and click **Fit**. The fit shouldn't change much — this is asserting a flat detector.

### 5. Cake and integration

Switch to the **Cake and Integration** tab and check intensity vs 2θ. The calibrant peaks should line up sharply at their expected positions. If they do:

- **Save as PONI file...** (point of normal incidence). This is your `.poni` file.

You now have everything you need:

- `your_calibration.poni`
- `your_mask.edf` (or `.npy`)
- (optional) `your_calibration.npt` — the picked-rings file, useful for re-fitting later

## Using them downstream

```python
from giwaxs_analysis import calibration

calib = calibration.load_calibration("your_calibration.poni", "your_mask.edf")
print(calibration.summarise(calib))
```

Pass the `calib` object into every integration / plotting call from there on.

## Troubleshooting

- **"My peaks don't sit at the expected q"** — usually a wrong beam energy or wrong calibrant; double-check the experiment settings tab.
- **"The fit moves the beam centre far off"** — too few picked points on the inner rings, or stray picks on noise; redo peak picking.
- **"Rotations come out non-zero"** — your detector probably isn't actually flat; either accept it, or set them to zero and re-fit (with the caveat that you're then assuming a flat detector).
