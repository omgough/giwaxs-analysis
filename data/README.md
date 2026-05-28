# data/

This folder is for **small, redistributable sample data** that makes the example notebooks runnable out of the box. Suggested contents:

- One representative detector frame (a few MB, e.g. an `.edf` or `.tif`)
- The matching `.poni` file
- The matching mask
- (Optional) a tiny in-situ stack (~10 frames) for the heatmap example

If your sample data is too big to commit (> ~10 MB total), drop a `download_sample.py` script here instead that pulls from a Zenodo deposit or your group's data archive. Document the expected file layout in this README.

Large / raw / scratch data should **not** live here — it's `.gitignore`d. Keep that on shared storage or your beamline's data system.
