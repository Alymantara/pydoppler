import json

import numpy as np
import pytest

import pydoppler
from pydoppler.runtime import _status, covariance_modes, phase_mask, read_input


def test_covariance_inverse_and_whitening_are_unit_invariant():
    covariance = np.array([[4., 1.], [1., 2.]])
    residual = np.array([2., -1.])
    metrics = []
    for scale in (1., 1e-15, 1e3):
        values, vectors = covariance_modes(covariance * scale**2)
        whitened = vectors.T @ (residual * scale) / np.sqrt(values)
        metrics.append(whitened @ whitened)
    assert metrics == pytest.approx([residual @ np.linalg.solve(covariance, residual)] * 3)


def test_covariance_rejects_negative_variance_and_handles_null_modes():
    with pytest.raises(ValueError, match="positive semidefinite"):
        covariance_modes(np.array([[1., 2.], [2., 1.]]))
    values, vectors = covariance_modes(np.ones((2, 2)))
    assert len(values) == 1
    assert vectors.shape == (2, 1)


def test_phase_mask_wraps_and_preserves_open_boundaries():
    mask = phase_mask(np.array([.0, .04, .05, .5, .95, .99, 1.5, -.5]), .95, 1.05)
    assert mask.tolist() == [False, False, True, True, True, False, True, True]
    assert phase_mask(np.arange(4) / 4, 0, 0).all()


@pytest.mark.parametrize("field,value,match", [
    ("fit_statistic", "NaN", "nonfinite"),
    ("converged", "0", "did not converge"),
    ("optimizer_status", "11", "did not converge"),
    ("fit_statistic", "2.0", "did not reach"),
])
def test_status_rejects_invalid_and_unconverged_results(tmp_path, field, value, match):
    entries = dict(converged="1", optimizer_status="0", fit_statistic=".5",
                   fit_target="1", iterations="20")
    entries[field] = value
    path = tmp_path / "dop.status"
    path.write_text("\n".join(f"{k} {v}" for k, v in entries.items()))
    with pytest.raises(RuntimeError, match=match):
        _status(path)


def test_status_distinguishes_optimizer_convergence_from_missed_fit_target(tmp_path):
    path = tmp_path / "dop.status"
    path.write_text("converged 0\noptimizer_status 1\nfit_statistic 2\n"
                    "fit_target 1\niterations 30\n")
    with pytest.raises(RuntimeError, match="Optimizer converged, but the alpha schedule"):
        _status(path)


def _synthetic_input(dop, scale=1.):
    # Independent analytic S-wave in the documented PUBLIC coordinate system.
    phases = np.arange(12) / 12 + .017
    velocities = np.linspace(-600, 600, 49)
    center = -150 * np.cos(2 * np.pi * phases) + 220 * np.sin(2 * np.pi * phases)
    flux = scale * np.exp(-.5 * ((velocities[:, None] - center) / 60)**2)
    widths = np.full(len(phases), .002)
    with (dop.workdir / "dopin").open("w") as handle:
        handle.write(f"{len(phases)} {len(velocities)} 6562.8\n0 0 0 synthetic\n")
        np.savetxt(handle, phases[None, :], fmt="%.16e")
        handle.write("1\n")
        np.savetxt(handle, widths[None, :], fmt="%.16e")
        np.savetxt(handle, (velocities * 1e5)[None, :], fmt="%.16e")
        np.savetxt(handle, flux, fmt="%.16e")
    return flux


def test_rapid_rerun_executes_changed_parameters_and_records_failure(tmp_path):
    dop = pydoppler.spruit(workdir=tmp_path)
    dop.verbose = False
    dop.map_size = 20
    dop.ns = 3
    dop.nim = 2000
    dop.ac = 1e-3
    dop.pb0 = dop.pb1 = 0
    dop.clim = 1.6
    _synthetic_input(dop)
    dop.Syncdop()
    first = json.loads((tmp_path / "run.json").read_text())
    assert first["success"]
    assert first["solver"]["backend"] == "numpy"
    assert len(first["solver"]["source_sha256"]) == 64
    assert {"python", "numpy", "scipy"} <= set(first["solver"])
    dop.nim = 1
    dop.ac = 1e-12
    with pytest.raises(RuntimeError, match="converge"):
        dop.Syncdop()
    second = json.loads((tmp_path / "run.json").read_text())
    assert not second["success"]
    assert second["parameters"] != first["parameters"]
    assert dop.fit_diagnostics is None
    with pytest.raises(RuntimeError):
        dop.Dopmap(plot=False)


def test_run_timeout_stops_the_inversion_and_records_failure(tmp_path):
    dop = pydoppler.spruit(workdir=tmp_path)
    dop.verbose = False
    dop.map_size = 20
    dop.clim = 1.6
    _synthetic_input(dop)
    dop.run_timeout = 1e-9
    with pytest.raises(RuntimeError, match="run_timeout"):
        dop.Syncdop()
    manifest = json.loads((tmp_path / "run.json").read_text())
    assert not manifest["success"]
    assert not (tmp_path / "dop.out").exists()
    dop.run_timeout = 600.0
    assert dop.Syncdop()["converged"] == 1


def test_thread_count_is_validated_and_does_not_change_the_map(tmp_path):
    dop = pydoppler.spruit(workdir=tmp_path)
    dop.verbose = False
    dop.map_size = 20
    dop.clim = 1.6
    _synthetic_input(dop)
    dop.threads = 0
    with pytest.raises(ValueError, match="threads"):
        dop.Syncdop()
    maps = []
    for threads in (1, 3):
        dop.threads = threads
        dop.Syncdop()
        maps.append((tmp_path / "dop.out").read_bytes())
    assert maps[0] == maps[1]


def test_read_input_retains_tiny_physical_flux(tmp_path):
    dop = pydoppler.spruit(workdir=tmp_path)
    expected = _synthetic_input(dop, scale=1e-15)
    _, _, _, actual, _ = read_input(tmp_path / "dopin")
    np.testing.assert_allclose(actual, expected, rtol=1e-15, atol=0)


def test_residuals_require_real_uncertainties(tmp_path):
    dop = pydoppler.spruit(workdir=tmp_path)
    data = np.array([[2., 3.], [4., 5.]])
    model = data - 1
    residual, standardized = dop.Residuals(dm=data, dmr=model)
    np.testing.assert_array_equal(residual, np.ones_like(data))
    assert standardized is None
    _, standardized = dop.Residuals(dm=data, dmr=model, sigma=np.full_like(data, .5))
    np.testing.assert_array_equal(standardized, np.full_like(data, 2.))
    with pytest.raises(ValueError, match="positive"):
        dop.Residuals(dm=data, dmr=model, sigma=np.zeros_like(data))
