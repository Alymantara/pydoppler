"""Independent statistical and optimizer checks of the numerical core."""
import time

import numpy as np
import pytest
from scipy.special import ndtr

from pydoppler import solver

PHASES = np.array([0.0, 0.25, 0.5, 0.75])
VELOCITY = np.linspace(-1000.0, 1000.0, 41)


def _swave(phases):
    centre = 300*np.sin(2*np.pi*np.asarray(phases)[:, None])
    return 1 + 4*np.exp(-0.5*((VELOCITY[None, :] - centre)/160)**2)


def _run(*, scale=1.0, iw=1, precision=None, norm=0, phases=None, flux=None, nim=1000,
         ac=8e-4, alpha=.003, mask=(0.0, 0.0), ih=0, nal=0, clim=99.0):
    """Invert four spectra of an analytic S-wave; arrays have shape (nph, nvp)."""
    phases = PHASES if phases is None else np.asarray(phases)
    flux = (_swave(phases) if flux is None else np.asarray(flux))*scale
    sigma = np.full(flux.shape, 0.3*scale)
    ranks = None
    if iw == 2:
        if precision is None:
            precision = np.repeat((np.eye(41)/0.3**2)[None], 4, axis=0)
        precision = np.asarray(precision)/scale**2
        ranks = np.full(4, 41)
    result = solver.solve(
        phases, np.full(4, .001), VELOCITY*1e5, flux, sigma if iw else None,
        precision if iw == 2 else None, ranks, ih=ih, iw=iw, pb0=mask[0], pb1=mask[1],
        ns=3, ac=ac, nim=nim, al0=alpha, alf=1.7, nal=nal, clim=clim, norm=norm,
        wid=1.0, af=0.0, map_size=20)
    return result, flux, sigma


def test_gaussian_statistic_matches_independent_residual_calculation():
    result, flux, sigma = _run()
    assert result.converged
    expected = np.sqrt(np.mean(((flux - result.model)/sigma)**2))
    assert result.status["fit_statistic"] == pytest.approx(expected, rel=1e-10)
    assert result.status["likelihood"] == 1  # Supplied sigma implies Gaussian data.
    assert result.status["stationarity"] <= 8e-4


def test_flux_units_do_not_change_solution_or_statistic():
    reference, _, _ = _run()
    scaled, _, _ = _run(scale=1e-15)
    assert reference.converged and scaled.converged
    assert scaled.status["fit_statistic"] == pytest.approx(
        reference.status["fit_statistic"], rel=1e-6)
    np.testing.assert_allclose(scaled.image/1e-15, reference.image,
                               rtol=2e-4, atol=np.max(np.abs(reference.image))*2e-4)


def test_full_precision_matches_generalized_least_squares():
    covariance = 0.3**2 * 0.35**np.abs(np.arange(41)[:, None]-np.arange(41)[None, :])
    precision = np.repeat(np.linalg.inv(covariance)[None], 4, axis=0)
    result, flux, _ = _run(iw=2, precision=precision, alpha=.03)
    assert result.converged
    residual = flux - result.model
    expected = np.sqrt(np.einsum("pi,pij,pj->", residual, precision, residual)/flux.size)
    assert result.status["fit_statistic"] == pytest.approx(expected, rel=1e-10)


def test_excluded_measurements_do_not_affect_map():
    base, flux, _ = _run(mask=(0.2, 0.3))
    changed_flux = flux.copy()
    changed_flux[1] = -1e6*np.arange(1, 42)
    altered, _, _ = _run(mask=(0.2, 0.3), flux=changed_flux)
    assert base.converged and altered.converged
    assert base.status["statistic_rank"] == 3*41
    assert altered.status["fit_statistic"] == pytest.approx(
        base.status["fit_statistic"], abs=1e-12)
    np.testing.assert_array_equal(base.image, altered.image)


def test_iteration_limit_is_explicit_failure():
    result, _, _ = _run(nim=1)
    assert not result.converged
    assert result.status["optimizer_status"] != 0
    assert result.image is None and result.model is None


def test_normalization_preserves_raw_diagnostic_data_and_gaussian_statistic():
    raw = _swave(PHASES)*np.array([0.6, 1.0, 1.3, 1.8])[:, None]
    result, _, sigma = _run(norm=1, flux=raw)
    assert result.converged
    np.testing.assert_allclose(result.data, raw, rtol=1e-12)
    expected = np.sqrt(np.mean(((raw - result.model)/sigma)**2))
    assert result.status["fit_statistic"] == pytest.approx(expected, rel=1e-10)
    np.testing.assert_allclose(result.normalization*raw.sum(axis=1), raw.sum(axis=1).mean())


def test_default_relative_rms_is_phase_permutation_invariant():
    base, flux, _ = _run(iw=0)
    order = np.array([2, 0, 3, 1])
    permuted, _, _ = _run(iw=0, phases=PHASES[order], flux=flux[order])
    assert base.converged and permuted.converged
    assert permuted.status["fit_statistic"] == pytest.approx(
        base.status["fit_statistic"], rel=1e-6)
    expected = np.sqrt(np.mean((flux - base.model)**2)/np.mean(flux**2))
    assert base.status["fit_statistic"] == pytest.approx(expected, rel=1e-10)
    np.testing.assert_allclose(base.image, permuted.image, rtol=1e-4, atol=1e-8)


def test_signed_samples_are_preserved_including_window_edges():
    flux = _swave(PHASES) - 1.2
    result, raw, _ = _run(flux=flux, alpha=0.0)
    assert result.converged
    # Gaussian fits keep the requested baseline: negative data are not allowed
    # to lift it (which formerly removed the positivity constraint).
    assert result.status["baseline"] == pytest.approx(0.1)
    assert result.status["baseline_raised"] == 0
    assert np.any(result.data < 0)
    np.testing.assert_allclose(result.data, raw, rtol=1e-12, atol=1e-12)


def test_phase_mask_endpoints_remain_open_across_phase_wrap():
    result, _, _ = _run(phases=np.array([.05, .2, .6, .95]), mask=(.95, 1.05), nim=1)
    assert result.status["statistic_rank"] == 4*41


def test_objective_derivatives_and_projection_adjoint():
    """Check derivatives against independent finite differences, including covariance."""
    nph, nvp = 4, 41
    precision = np.repeat(
        (2*np.eye(nvp) - 0.25*(np.eye(nvp, k=1) + np.eye(nvp, k=-1)))[None], nph, axis=0)
    problem = solver.Tomography(
        np.arange(nph)/nph, np.full(nph, 0.01/solver.PI2), -1000 + 50.0*np.arange(nvp),
        np.ones((nph, nvp)), np.ones((nph, nvp)), precision, np.full(nph, nvp),
        iw=2, pb0=0.2, pb1=0.3, ns=3, map_size=20)
    # Fix the statistical model directly, independently of the input scaling.
    problem.precision = precision.copy()
    problem.likelihood_scale = 0.7
    problem.build_projection()
    assert problem.phase_valid.tolist() == [True, False, True, True]

    pixel = np.arange(1, problem.nol + 1)
    datum = np.arange(1, problem.nil + 1)
    ps = 2 + 0.5*np.sin(pixel)
    direction = 0.2*np.cos(0.3*pixel)
    pht = 1 + 0.2*np.cos(datum)
    testdata = np.sin(0.7*datum)
    assert problem.fproj(direction) @ testdata == pytest.approx(
        direction @ problem.bproj(testdata), rel=1e-12)

    def evaluate(x):
        h, dha = problem.like(pht, problem.fproj(x))
        s, dsa, chi, tem = problem.entr(x)
        return h, s, dha, dsa, chi, tem

    al = 0.7
    h, s, dha, dsa, chi, tem = evaluate(ps)
    _, fmu = problem.lamun(pht, ps, direction, direction, problem.fproj(ps), chi, tem, al)
    hessian = -(direction @ dha + al*(direction @ dsa))/fmu
    eps = 1.0e-3
    hplus, splus = evaluate(ps + eps*direction)[:2]
    hminus, sminus = evaluate(ps - eps*direction)[:2]
    np.testing.assert_allclose(
        [direction @ dha, direction @ dsa, hessian],
        [(hplus - hminus)/(2*eps), (splus - sminus)/(2*eps),
         (hplus + al*splus - 2*(h + al*s) + hminus + al*sminus)/eps**2],
        rtol=2e-6, atol=1e-7)


def _projection_profile(phase, width, sigma):
    """One bright map pixel projected through a full exposure, on irregular bins."""
    nph, nvp = 4, 41
    x = -1 + 2*np.arange(nvp)/(nvp - 1)
    centers = 1000*x + 100*np.sin(np.pi*x)
    matrix, indx, indy, _, _ = solver.build_projection(
        np.full(nph, phase*solver.PI2), np.full(nph, width*solver.PI2), centers, 20,
        instrument_sigma=sigma, phase_subsamples=1024)
    spot = np.where((indx == 10) & (indy == 5), float(nph), 0.0)
    assert spot.sum() == nph
    return centers, (matrix @ spot).reshape(nph, nvp)[0]


def _independent_square_pixel(centers, phase, width, sigma, vx=0.5, vy=-4.5):
    """Gauss-Legendre average over exposure phase and a uniform square pixel."""
    edges = np.r_[centers[0]-(centers[1]-centers[0])/2, (centers[:-1]+centers[1:])/2,
                  centers[-1]+(centers[-1]-centers[-2])/2]
    pixel = 2000/19
    phase_nodes, phase_weight = np.polynomial.legendre.leggauss(256)
    pix_nodes, pix_weight = np.polynomial.legendre.leggauss(64)
    angle = 2*np.pi*(phase + width*phase_nodes/2)
    ux = pixel/2*pix_nodes
    # Internal velocity components of the pixel with zero-based grid indices
    # (10, 5) on a 20-pixel grid, in pixel units: index - 9.5.
    centre_x, centre_y = pixel*vx, pixel*vy
    radial = ((centre_x + ux)[None, :, None]*np.sin(angle)[:, None, None]
              - (centre_y + ux)[None, None, :]*np.cos(angle)[:, None, None])
    if sigma > 0:
        cdf = ndtr((edges[:, None, None, None] - radial[None])/sigma)
    else:
        cdf = (edges[:, None, None, None] >= radial[None]).astype(float)
    weights = (phase_weight[:, None, None]/2 * pix_weight[None, :, None]/2
               * pix_weight[None, None, :]/2)
    return np.einsum("vpab,pab->v", np.diff(cdf, axis=0), weights)


def test_exposure_quadrature_lsf_and_irregular_velocity_bins():
    """Compare the forward operator with independent phase/pixel quadrature."""
    centers, actual = _projection_profile(0.0, 0.5, 40.0)
    expected = _independent_square_pixel(centers, 0.0, 0.5, 40.0)
    np.testing.assert_allclose(actual, expected, atol=5e-7, rtol=2e-5)
    assert actual.sum() == pytest.approx(1.0, abs=1e-7)
    assert centers @ actual > 250  # Endpoint-only smearing formerly gave zero.


def test_square_pixels_project_without_ripple_at_cardinal_phases():
    """A uniform pixel seen at phase 0 projects as a box; the old triangle did not."""
    centers, actual = _projection_profile(0.0, 1e-12, 0.0)
    edges = np.r_[centers[0]-(centers[1]-centers[0])/2, (centers[:-1]+centers[1:])/2,
                  centers[-1]+(centers[-1]-centers[-2])/2]
    pixel = 2000/19
    lo, hi = 4.0*pixel, 5.0*pixel  # v = -Vy_internal at phase 0, Vy_internal = -4.5 px
    expected = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0, None)/pixel
    np.testing.assert_allclose(actual, expected, atol=1e-9)
    assert actual.sum() == pytest.approx(1.0, abs=1e-12)


def test_alpha_schedule_brackets_target_and_reports_plateau():
    """The schedule refines a crossing by bisection and detects infeasible targets."""
    noisy = _swave(PHASES) + 0.3*np.random.default_rng(4).normal(size=(4, 41))
    result, _, _ = _run(flux=noisy, alpha=30.0, nal=40, clim=1.1)
    assert result.converged
    assert result.status["target_bracketed"] == 1
    assert result.status["alpha_steps"] > 2
    assert 1.1*(1 - 0.006) <= result.status["fit_statistic"] <= 1.1
    result, _, _ = _run(flux=noisy, alpha=30.0, nal=40, clim=0.2)
    assert not result.converged
    assert result.status["optimizer_status"] == 1
    assert result.status["rms_plateau"] == 1
    assert result.status["rms_floor_estimate"] > 0.2


def test_fast_pixel_response_equals_the_general_formula():
    """The in-place trapezoid used without an LSF matches square_cdf exactly."""
    rng = np.random.default_rng(1)
    ha = rng.uniform(0.0, 2.0, 200)
    hb = rng.uniform(0.0, 2.0, 200)
    ha[:3], hb[:3] = [1.0, 0.0, 1.0], [0.0, 1.0, 1e-9]  # box limits at cardinal phases
    x = rng.uniform(-5.0, 5.0, (200, 7))
    expected = solver.square_cdf(x, ha[:, None], hb[:, None])
    np.testing.assert_array_equal(solver._trapezoid_cdf_inplace(x.copy(), ha, hb), expected)
    assert expected.min() == 0.0 and expected.max() == 1.0


def test_projection_does_not_depend_on_threads_or_work_blocks(monkeypatch):
    phases = np.array([0.03, 0.31, 0.62, 0.9])*solver.PI2
    widths = np.array([0.08, 0.002, 0.04, 0.12])*solver.PI2
    arguments = (phases, widths, VELOCITY, 24)
    reference, indx, indy, va, pixel = solver.build_projection(*arguments, threads=1)
    assert reference.nnz > 0
    monkeypatch.setattr(solver, "_BLOCK_ELEMENTS", 2000)
    for threads in (1, 3):
        matrix = solver.build_projection(*arguments, threads=threads)[0]
        np.testing.assert_array_equal(matrix.indptr, reference.indptr)
        np.testing.assert_array_equal(matrix.indices, reference.indices)
        np.testing.assert_array_equal(matrix.data, reference.data)
    # The storage bound only saves copies: a far too small one changes nothing.
    monkeypatch.setattr(solver, "_capacity_bound", lambda *arguments: 1)
    regrown = solver.build_projection(*arguments)[0]
    np.testing.assert_array_equal(regrown.indices, reference.indices)
    np.testing.assert_array_equal(regrown.data, reference.data)
    # A pixel projected wholly inside the spectral window distributes unit
    # weight over the spectra; pixels at the rim lose part of theirs.
    weight = np.asarray(reference.sum(axis=0)).ravel()
    radius = pixel*np.hypot(indx - 11.5, indy - 11.5)
    assert np.any(radius < 0.8*va)
    np.testing.assert_allclose(weight[radius < 0.8*va], 1.0, rtol=1e-12)
    assert weight.max() <= 1 + 1e-12 and weight.min() < 0.99


def test_expired_deadline_stops_the_run():
    arguments = (PHASES, np.full(4, .001), VELOCITY*1e5, _swave(PHASES))
    with pytest.raises(solver.SolverTimeout):
        solver.solve(*arguments, iw=0, map_size=20, deadline=time.monotonic() - 1)
    problem = solver.Tomography(*arguments, iw=0, map_size=20)
    problem.build_projection()
    problem.set_weights()
    problem.deadline = time.monotonic() - 1
    with pytest.raises(solver.SolverTimeout):
        problem.optimize(problem.dat.ravel() + 0.1, 0.003, np.ones(problem.nol))
