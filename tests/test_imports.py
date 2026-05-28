"""Smoke tests: the package imports and exposes its public API.

These don't exercise any actual analysis — they just confirm the
package structure is healthy. As you port functions into the modules,
add real unit tests next to these (use the sample dataset in `data/`).
"""

from __future__ import annotations


def test_package_imports():
    import giwaxs_analysis

    assert giwaxs_analysis.__version__ != ""


def test_submodules_present():
    import giwaxs_analysis

    expected = {"calibration", "io", "integration", "plotting", "peak_fitting"}
    assert expected.issubset(set(dir(giwaxs_analysis)))


def test_public_callables_exist():
    """The functions advertised in the README should be importable, even if NotImplemented."""
    from giwaxs_analysis import calibration, integration, io, peak_fitting, plotting

    assert callable(calibration.load_calibration)
    assert callable(io.load_frame)
    assert callable(io.load_stack)
    assert callable(io.list_frames)
    assert callable(integration.radial_integrate)
    assert callable(integration.sector_integrate)
    assert callable(integration.cake)
    assert callable(integration.batch_integrate)
    assert callable(plotting.plot_detector)
    assert callable(plotting.plot_1d)
    assert callable(plotting.plot_heatmap)
    assert callable(peak_fitting.fit_peak)
    assert callable(peak_fitting.scherrer)
