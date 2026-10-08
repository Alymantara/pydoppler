import numpy as np
import pytest

import pydoppler
from pydoppler.preprocessing import (
    SpectralCovariance,
    estimate_continuum_sigma,
    interpolation_operator,
    subtract_continuum_and_resample,
)


def prepared(tmp_path, scale=1.0):
    dop = pydoppler.spruit(workdir=tmp_path)
    dop.verbose = False
    wave = np.linspace(6500., 6600., 101)
    dop.wave = [wave]
    dop.flux = np.array([(10. + np.exp(-((wave - 6550.) / 3.)**2)) * scale])
    dop.pha = dop.input_phase = np.array([.2])
    dop.lam0, dop.delw = 6550., 20.
    return dop


def test_interpolated_errors_and_neighbor_covariance():
    sampling = interpolation_operator([0., 2.], [.5, 1., 1.5])
    cov = SpectralCovariance(np.ones(2), sampling)
    assert cov.diagonal() == pytest.approx([.625, .5, .625])
    assert cov.to_dense()[0, 1] == pytest.approx(.5)
    second = interpolation_operator([.5, 1., 1.5], [.75, 1.25])
    composed = cov.transform(second)
    assert composed.diagonal() == pytest.approx([.53125, .53125])
    assert composed.to_dense()[0, 1] == pytest.approx(.46875)


def test_joint_continuum_and_interpolation_covariance_matches_noise_draws():
    wave = np.arange(21.)
    target = np.array([7.25, 7.75, 10.25, 13.75])
    bands = [-1., 5.5, 14.5, 21.]
    sigma = np.linspace(.5, 1.5, len(wave))
    _, covariance = subtract_continuum_and_resample(
        wave, np.ones(len(wave)), target, bands, 1,
        SpectralCovariance(sigma**2),
    )
    rng = np.random.default_rng(130)
    noise = rng.normal(size=(120000, len(wave))) * sigma
    cont = (wave < 5.5) | (wave > 14.5)
    # An independently constructed weighted straight-line fit to noise draws.
    basis = np.column_stack([np.ones(np.count_nonzero(cont)), wave[cont]])
    coeff = np.linalg.lstsq(basis / sigma[cont, None],
                           (noise[:, cont] / sigma[cont]).T, rcond=None)[0]
    interpolated = np.stack([np.interp(target, wave, row) for row in noise])
    output = interpolated - (np.column_stack([np.ones(len(target)), target]) @ coeff).T
    empirical = np.cov(output, rowvar=False)
    expected = covariance.to_dense()
    assert np.allclose(empirical, expected, rtol=.035, atol=.005)
    assert np.allclose(covariance.diagonal(), np.diag(expected), rtol=1e-12)
    assert expected[0, 2] > 0  # common continuum-fit uncertainty couples pixels


def test_scientific_notation_preserves_physical_flux_and_errors(tmp_path):
    dop = prepared(tmp_path, scale=1e-15)
    dop.iw = 1
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6500, 6520, 6580, 6600],
              flux_err=np.full((1, 101), 1e-16), estimate_errors=False)
    lines = (tmp_path / "dopin").read_text().splitlines()
    body = np.fromstring(" ".join(lines[2:]), sep=" ")
    begin = 3 + len(dop.vell)  # phase, iph, width, velocities
    flux = body[begin:begin + len(dop.vell)]
    sigma = body[begin + len(dop.vell):]
    assert np.allclose(flux, dop.normalised_flux[0], rtol=1e-14, atol=0.)
    assert np.allclose(sigma, dop.normalised_flux_err[0], rtol=1e-14, atol=0.)
    assert np.max(flux) > 9e-16
    assert np.all(sigma > 0)


def test_continuum_windows_stay_in_observed_frame_at_nonzero_gamma(tmp_path):
    dop = prepared(tmp_path)
    dop.gama = 200.
    dop.delw = 30.
    dop.flux[0, 6] = 1000.  # outside observed window, inside the old shifted window
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6500, 6504, 6596, 6600])
    factor = np.sqrt((1. + dop.gama / 299792.458) / (1. - dop.gama / 299792.458))
    expected = np.interp(dop.normalised_wave * factor, dop.wave[0], dop.flux[0]) - 10.
    assert np.allclose(dop.normalised_flux[0], expected, atol=1e-12)


def test_redshifted_coverage_need_not_contain_the_rest_wavelength(tmp_path):
    dop = prepared(tmp_path)
    dop.gama = 4000.
    factor = np.sqrt((1. + dop.gama / 299792.458) / (1. - dop.gama / 299792.458))
    dop.wave = [dop.wave[0] * factor]
    assert dop.wave[0][0] > dop.lam0
    dop.Dopin(plot=False, poly_degree=1,
              continuum_band=np.array([6500, 6520, 6580, 6600]) * factor)
    assert np.max(dop.normalised_flux) > .99


def test_invalid_errors_are_rejected_when_estimation_disabled(tmp_path):
    dop = prepared(tmp_path)
    dop.iw = 1
    with pytest.raises(ValueError, match="finite and positive"):
        dop.Dopin(plot=False, flux_err=np.full((1, 101), np.nan), estimate_errors=False)


def test_noiseless_continuum_does_not_invent_unit_errors(tmp_path):
    dop = prepared(tmp_path)
    dop.iw = 1
    with pytest.raises(ValueError, match="no measurable noise"):
        dop.Dopin(plot=False)


def test_continuum_noise_estimator_uses_residual_degrees_of_freedom():
    wave = np.arange(10.)
    flux = np.array([0., 1., -1., 1., -1., 1., -1., 1., -1., 0.])
    sigma = estimate_continuum_sigma(wave, flux, [-1, 4.5, 4.5, 10], 0)
    assert sigma == pytest.approx(np.sqrt(8. / 9.))


def test_native_resampling_and_explicit_flux_mutation(tmp_path):
    coarse = np.arange(6500., 6601., 2.)
    fine = np.arange(6500., 6600., .7)
    coarse_flux = 10. + (coarse == 6550.)
    np.savetxt(tmp_path / "a", np.c_[coarse, coarse_flux, np.full(len(coarse), .1)])
    np.savetxt(tmp_path / "b", np.c_[fine, np.full(len(fine), 10.), np.full(len(fine), .1)])
    (tmp_path / "phases").write_text("a .1\nb .6\n")
    dop = pydoppler.spruit(workdir=tmp_path / "work")
    dop.verbose = False
    dop.base_dir, dop.list = tmp_path, "phases"
    dop.lam0, dop.delw = 6550., 15.
    dop.Foldspec()
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6500, 6520, 6580, 6599])
    assert dop.preprocessing_used_native
    expected = np.interp(dop.normalised_wave, coarse, coarse_flux) - 10.
    assert np.allclose(dop.normalised_flux[0], expected, atol=1e-12)
    before = dop.normalised_flux.copy()
    dop.flux[0] = np.asarray(dop.flux[0]) * 2.
    dop.Dopin(plot=False, poly_degree=1, continuum_band=[6500, 6520, 6580, 6599])
    assert not dop.preprocessing_used_native
    assert np.max(dop.normalised_flux[0]) > 1.5 * np.max(before[0])
    covariance = dop.normalised_flux_covariance[0]
    assert np.allclose(covariance.diagonal(), np.diag(covariance.to_dense()))


@pytest.mark.parametrize("field,value", [("gama", np.nan), ("delw", -1.),
                                         ("delta_phase", 0.), ("lam0", np.inf)])
def test_invalid_scientific_parameters_raise(tmp_path, field, value):
    dop = prepared(tmp_path)
    setattr(dop, field, value)
    with pytest.raises(ValueError):
        dop.Dopin(plot=False)


def test_dopin_can_repeat_without_stale_uncertainties(tmp_path):
    dop = prepared(tmp_path)
    dop.Dopin(plot=False, flux_err=np.full((1, 101), .1))
    assert dop.normalised_flux_err is not None
    dop.Dopin(plot=False)
    assert dop.normalised_flux_err is None
    assert dop.normalised_flux_covariance == [None]
