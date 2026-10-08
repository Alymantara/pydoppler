"""End-to-end checks of covariance, provenance snapshots, and allocation limits."""

import json
import shutil

import numpy as np
import pytest

import pydoppler
from pydoppler import runtime
from pydoppler.pydoppler import _parse_dopout


def _prepared(dopdir):
    pydoppler.copy_fortran_code(dopdir)
    dop = pydoppler.spruit(workdir=dopdir)
    dop.verbose = False
    dop.iw = 1
    dop.map_size = 20
    dop.ns = 3
    dop.nim = 4000
    dop.ac = 1e-3
    dop.pb0 = dop.pb1 = 0.
    dop.lam0, dop.delw = 6550., 20.
    wave = np.linspace(6510., 6590., 101)
    ratio = wave / dop.lam0
    velocity = 299792.458 * (ratio**2 - 1.) / (ratio**2 + 1.)
    phases = np.arange(8) / 8 + .013
    center = -150 * np.cos(2 * np.pi * phases) + 220 * np.sin(2 * np.pi * phases)
    dop.wave = [wave]
    dop.flux = 10. + np.exp(-.5 * ((velocity[None, :] - center[:, None]) / 110.)**2)
    dop.pha = dop.input_phase = phases
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6510, 6525, 6575, 6590],
              flux_err=np.full_like(dop.flux, .1), estimate_errors=False)
    return dop


@pytest.mark.skipif(not shutil.which("make") or not shutil.which("gfortran"),
                    reason="Fortran toolchain required")
def test_full_covariance_run_matches_independent_residual_statistic_and_snapshots(tmp_path):
    dop = _prepared(tmp_path)
    diagnostics = dop.Syncdop()
    assert diagnostics["weight_mode"] == 2
    residuals, whitened = dop.Residuals(whiten=True)
    rank = sum(dop.covariance_ranks)
    measured = np.sqrt(np.sum(whitened**2) / rank)
    assert measured == pytest.approx(diagnostics["fit_statistic"], rel=1e-10)
    _, standardized = dop.Residuals()
    manifest = json.loads((tmp_path / "run.json").read_text())
    assert manifest["success"]
    assert manifest["effective_iw"] == 2
    assert len(manifest["precision_sha256"]) == 64

    # Changing preparation and settings must not change an already fitted run's
    # uncertainty model or phase mask. Deep copying also guards in-place edits.
    dop.normalised_flux_covariance[0].variance *= 4.
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6510, 6525, 6575, 6590],
              flux_err=np.full_like(dop.flux, .2), estimate_errors=False)
    dop.normalised_flux_covariance[0].variance *= 9.
    dop.covariance_rcond = .9
    dop.valid_phase_mask[:] = False
    after_residuals, after_whitened = dop.Residuals(whiten=True)
    _, after_standardized = dop.Residuals()
    np.testing.assert_allclose(after_residuals, residuals, rtol=0, atol=0)
    np.testing.assert_allclose(after_whitened, whitened, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(after_standardized, standardized, rtol=0, atol=0)


def test_memory_limit_is_checked_before_dense_covariance_preparation(tmp_path, monkeypatch):
    dop = _prepared(tmp_path)
    dop.max_memory_bytes = 1
    monkeypatch.setattr(runtime.shutil, "which", lambda name: "/available/" + name)

    def forbidden_dense():
        raise AssertionError("Dense covariance was allocated before the memory guard.")

    monkeypatch.setattr(dop.normalised_flux_covariance[0], "to_dense", forbidden_dense)
    with pytest.raises(MemoryError, match="before preparation"):
        dop.Syncdop()
    assert not (tmp_path / "dop.precision").exists()
    assert not json.loads((tmp_path / "run.json").read_text())["success"]


def test_weighted_data_normalization_requires_explicit_conditional_assumption(tmp_path, monkeypatch):
    dop = _prepared(tmp_path)
    dop.norm = 1
    monkeypatch.setattr(runtime.shutil, "which", lambda name: "/available/" + name)
    with pytest.raises(ValueError, match="conditional_normalization"):
        dop.Syncdop()
    assert not (tmp_path / "dop.precision").exists()
    dop.conditional_normalization = True
    parameters = runtime._parameters(dop, effective_iw=2)
    assert "opt-in flat integrated flux normalization" in parameters


def test_preprocessing_provenance_matches_actual_arrays_and_rejects_changed_input(tmp_path):
    from pydoppler.preprocessing import array_digest

    dop = _prepared(tmp_path)
    metadata = json.loads((tmp_path / "dopin.json").read_text())
    assert metadata["continuum_frame"] == "observed_wavelength"
    assert metadata["continuum_bands_angstrom"] == [6510., 6525., 6575., 6590.]
    assert metadata["continuum_polynomial_degree"] == 1
    assert metadata["source_grid"] == "loaded_common"
    assert metadata["sources"][0]["flux_sha256"] == array_digest(dop.flux[0])
    assert metadata["sources"][0]["wavelength_sha256"] == array_digest(dop.wave[0])
    assert metadata["sources"][0]["noise_source"] == "provided_standard_deviations"
    current_hash = runtime._sha(tmp_path / "dopin")
    matched = runtime.preprocessing_provenance(tmp_path, current_hash)
    assert matched["status"] == "matched"
    assert matched["metadata"] == metadata
    with (tmp_path / "dopin").open("a") as handle:
        handle.write("\n")  # Any external edit invalidates the sidecar's exact-file link.
    unavailable = runtime.preprocessing_provenance(tmp_path, runtime._sha(tmp_path / "dopin"))
    assert unavailable["status"] == "unavailable"
    assert "changed" in unavailable["reason"]


@pytest.mark.skipif(not shutil.which("make") or not shutil.which("gfortran"),
                    reason="Fortran toolchain required")
def test_diagonal_manual_input_does_not_reuse_stale_preparation_covariance(tmp_path):
    dop = _prepared(tmp_path)
    phases, widths, velocity, flux, sigma = runtime.read_input(tmp_path / "dopin", weighted=True)
    sigma *= 2.  # An external/manual input update leaves wrapper covariance stale.
    header = (tmp_path / "dopin").read_text().splitlines()[:2]
    with (tmp_path / "dopin").open("w") as handle:
        handle.write("\n".join(header) + "\n")
        np.savetxt(handle, phases[None, :], fmt="%.17e")
        handle.write("1\n")
        np.savetxt(handle, widths[None, :], fmt="%.17e")
        np.savetxt(handle, velocity[None, :], fmt="%.17e")
        np.savetxt(handle, flux, fmt="%.17e")
        np.savetxt(handle, sigma, fmt="%.17e")
    dop.covariance_mode = "diagonal"
    diagnostics = dop.Syncdop()
    assert diagnostics["weight_mode"] == 1
    assert dop._fit_covariance is None
    residuals, standardized = dop.Residuals()
    np.testing.assert_allclose(standardized, residuals / sigma, rtol=1e-14, atol=0)
    with pytest.raises(ValueError, match="full-covariance fit"):
        dop.Residuals(whiten=True)


@pytest.mark.skipif(not shutil.which("make") or not shutil.which("gfortran"),
                    reason="Fortran toolchain required")
def test_four_quadrant_spots_gamma_sign_and_physical_flux_units(tmp_path):
    """Inject analytic S-waves, without using the backend's projection matrix."""
    dop = _prepared(tmp_path)
    dop.map_size = 24
    dop.delw = 19.9  # Avoid placing the window edge exactly on a native sample.
    rest_wave = np.linspace(6510., 6590., 101)
    light_speed = 299792.458
    ratio = rest_wave / dop.lam0
    native_velocity = light_speed * (ratio**2 - 1.) / (ratio**2 + 1.)
    phases = np.arange(12) / 12 + .013
    baseline_map = None
    baseline_statistic = None
    binary_mtime = None
    cases = [(150., 220., 90., 1.), (-150., 220., -120., 1.),
             (-150., -220., 65., 1.), (150., -220., -80., 1.),
             (150., 220., 90., 1e-15)]
    for vx, vy, gamma, flux_scale in cases:
        # Public phase convention, independent of both Python and Fortran
        # internal map-axis layout. Positive gamma redshifts wavelengths.
        center = -vx * np.cos(2 * np.pi * phases) + vy * np.sin(2 * np.pi * phases)
        redshift = np.sqrt((1. + gamma / light_speed) / (1. - gamma / light_speed))
        dop.wave = [rest_wave * redshift]
        dop.flux = (10. + np.exp(-.5 * (
            (native_velocity[None, :] - center[:, None]) / 110.
        )**2)) * flux_scale
        dop.pha = dop.input_phase = phases
        dop.gama = gamma
        dop.Dopin(plot=False, poly_degree=1,
                  continuum_band=np.array([6510, 6525, 6575, 6590]) * redshift,
                  flux_err=np.full_like(dop.flux, .1 * flux_scale),
                  estimate_errors=False)
        centroid = ((dop.normalised_flux * dop.vell).sum(axis=1)
                    / dop.normalised_flux.sum(axis=1))
        # This is much tighter than a spectral pixel, and catches applying
        # gamma twice or using the opposite sign even before inversion.
        np.testing.assert_allclose(centroid, center, rtol=0, atol=.1)
        diagnostics = dop.Syncdop()
        assert diagnostics["converged"] == 1
        assert diagnostics["weight_mode"] == 2
        output = _parse_dopout(tmp_path / "dop.out")
        row, col = np.unravel_index(np.argmax(output.im), output.im.shape)
        recovered = np.array([output.vmap[col], output.vmap[row]])
        map_pixel = output.vmap[1] - output.vmap[0]
        np.testing.assert_allclose(recovered, [vx, vy], rtol=0, atol=map_pixel / 2.)
        np.testing.assert_array_equal(np.sign(recovered), np.sign([vx, vy]))
        current_mtime = (tmp_path / "dopp").stat().st_mtime_ns
        if binary_mtime is None:
            binary_mtime = current_mtime
            baseline_map = output.im.copy()
            baseline_statistic = diagnostics["fit_statistic"]
        else:
            assert current_mtime == binary_mtime  # Same dimensions reuse the executable.
        if flux_scale != 1.:
            assert diagnostics["fit_statistic"] == pytest.approx(baseline_statistic, rel=1e-5)
            np.testing.assert_allclose(output.im / flux_scale, baseline_map,
                                       rtol=2e-4, atol=np.max(np.abs(baseline_map)) * 2e-4)
