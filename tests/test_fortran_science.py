"""Independent statistical and optimizer checks of the compiled scientific core."""
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from pydoppler.pydoppler import _parse_dopout
from pydoppler.runtime import write_precision_block, write_precision_header


@pytest.fixture(scope="module")
def scientific_core(tmp_path_factory):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the independent core checks")
    work = tmp_path_factory.mktemp("fortran_science_build")
    source = Path(__file__).resolve().parents[1] / "pydoppler" / "fortran_code"
    shutil.copy(source / "dop_core.f90", work / "dop_core.f90")
    shutil.copy(source / "dop.f90", work / "dop.f90")
    result = subprocess.run(
        [compiler, "-O0", "-g", "-fcheck=all", "-Wall",
         "dop_core.f90", "dop.f90", "-o", "dopp"],
        cwd=work, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Warning:" not in result.stderr.replace(
        "Warning: Equality comparison for REAL(8)", "")
    return work / "dopp"


def _compile_with_core(binary, tmp_path, name, program):
    """Compile a test program against the solver module, independently of the driver."""
    shutil.copy(binary.parent / "dop_core.f90", tmp_path / "dop_core.f90")
    (tmp_path / f"{name}.f90").write_text(program)
    result = subprocess.run(
        [shutil.which("gfortran"), "-O0", "-fcheck=all", "dop_core.f90", f"{name}.f90",
         "-o", name], cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    result = subprocess.run([f"./{name}"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def _run(binary, work, *, scale=1.0, iw=1, precision=None, norm=0,
         phases=None, flux=None, nim=1000, ac=8e-4, alpha=.003, mask=(0.0, 0.0),
         ih=0, nal=-1, clim=99.0):
    work.mkdir()
    phases = np.array([0.0, 0.25, 0.5, 0.75]) if phases is None else np.asarray(phases)
    velocity = np.linspace(-1000.0, 1000.0, 41)
    if flux is None:
        flux = 1 + 4*np.exp(-0.5*((velocity[None, :] -
                         300*np.sin(2*np.pi*phases[:, None]))/160)**2)
    flux = np.asarray(flux)*scale
    sigma = np.full(flux.shape, 0.3*scale)
    values = ["4 41 6562.8", "0 0 0 independent-science-test",
              " ".join(map(str, phases)), "1", ".001 .001 .001 .001",
              " ".join(f"{v:.17e}" for v in velocity*1e5),
              " ".join(f"{f:.17e}" for f in flux.T.ravel())]
    if iw:
        values.append(" ".join(f"{s:.17e}" for s in sigma.T.ravel()))
    (work / "dopin").write_text("\n".join(values)+"\n")
    (work / "dop.in").write_text(
        f"{ih}\n{iw}\n{mask[0]} {mask[1]}\n3\n{ac}\n{nim}\n"
        f"{alpha} 1.7 {nal}\n{clim}\n1\n{norm}\n1 0\n0 0 .1\n20\n"
    )
    if iw == 2:
        if precision is None:
            precision = np.repeat((np.eye(41)/0.3**2)[None], 4, axis=0)
        with (work / "dop.precision").open("wb") as handle:
            write_precision_header(handle, 4, 41)
            for matrix in precision:
                write_precision_block(handle, 41, matrix/scale**2)
    result = subprocess.run([str(binary)], cwd=work, capture_output=True, text=True)
    status = {}
    if (work / "dop.status").exists():
        for line in (work / "dop.status").read_text().splitlines():
            key, value = line.split()
            status[key] = float(value)
    parsed = _parse_dopout(work / "dop.out") if (work / "dop.out").exists() else None
    return result, status, parsed, flux, sigma


def test_gaussian_statistic_matches_independent_residual_calculation(scientific_core, tmp_path):
    result, status, output, flux, sigma = _run(scientific_core, tmp_path/"diagonal")
    assert result.returncode == 0, result.stdout+result.stderr
    residual = flux-output.dmr.T
    expected = np.sqrt(np.mean((residual/sigma)**2))
    assert status["fit_statistic"] == pytest.approx(expected, rel=2e-5)
    assert status["likelihood"] == 1  # Supplied sigma implies Gaussian data.
    assert status["stationarity"] <= 8e-4


def test_flux_units_do_not_change_solution_or_statistic(scientific_core, tmp_path):
    reference = _run(scientific_core, tmp_path/"base")
    scaled = _run(scientific_core, tmp_path/"scaled", scale=1e-15)
    assert reference[0].returncode == scaled[0].returncode == 0
    assert scaled[1]["fit_statistic"] == pytest.approx(reference[1]["fit_statistic"], rel=1e-6)
    np.testing.assert_allclose(scaled[2].im/1e-15, reference[2].im,
                               rtol=2e-4, atol=np.max(np.abs(reference[2].im))*2e-4)


def test_full_precision_matches_generalized_least_squares(scientific_core, tmp_path):
    covariance = 0.3**2 * 0.35**np.abs(np.arange(41)[:, None]-np.arange(41)[None, :])
    precision = np.repeat(np.linalg.inv(covariance)[None], 4, axis=0)
    result, status, output, flux, _ = _run(
        scientific_core, tmp_path/"correlated", iw=2, precision=precision, alpha=.03,
    )
    assert result.returncode == 0, result.stdout+result.stderr
    residual = flux-output.dmr.T
    expected = np.sqrt(np.einsum("pi,pij,pj->", residual, precision, residual)/flux.size)
    assert status["fit_statistic"] == pytest.approx(expected, rel=3e-5)


def test_excluded_measurements_do_not_affect_map(scientific_core, tmp_path):
    base = _run(scientific_core, tmp_path/"masked_base", mask=(0.2, 0.3))
    changed_flux = base[3].copy()
    changed_flux[1] = -1e6*np.arange(1, 42)
    altered = _run(scientific_core, tmp_path/"masked_changed", mask=(0.2, 0.3), flux=changed_flux)
    assert base[0].returncode == altered[0].returncode == 0
    assert base[1]["statistic_rank"] == 3*41
    assert altered[1]["fit_statistic"] == pytest.approx(base[1]["fit_statistic"], abs=1e-12)
    np.testing.assert_array_equal(base[2].im, altered[2].im)


def test_iteration_limit_is_explicit_failure(scientific_core, tmp_path):
    result, status, output, _, _ = _run(scientific_core, tmp_path/"limited", nim=1)
    assert result.returncode != 0
    assert status["converged"] == 0
    assert status["optimizer_status"] != 0
    assert output is None


def test_normalization_preserves_raw_diagnostic_data_and_gaussian_statistic(scientific_core, tmp_path):
    seed = _run(scientific_core, tmp_path/"normalization_seed")
    flux = seed[3]*np.array([0.6, 1.0, 1.3, 1.8])[:, None]
    result, status, output, raw, sigma = _run(
        scientific_core, tmp_path/"normalized", norm=1, flux=flux,
    )
    assert result.returncode == 0, result.stdout+result.stderr
    np.testing.assert_allclose(output.dm.T, raw, rtol=5e-5)
    expected = np.sqrt(np.mean(((raw-output.dmr.T)/sigma)**2))
    assert status["fit_statistic"] == pytest.approx(expected, rel=3e-5)


def test_default_relative_rms_is_phase_permutation_invariant(scientific_core, tmp_path):
    base = _run(scientific_core, tmp_path/"ordered", iw=0)
    order = np.array([2, 0, 3, 1])
    permuted = _run(scientific_core, tmp_path/"permuted", iw=0,
                    phases=np.array([0., .25, .5, .75])[order], flux=base[3][order])
    assert base[0].returncode == permuted[0].returncode == 0
    assert permuted[1]["fit_statistic"] == pytest.approx(base[1]["fit_statistic"], rel=1e-6)
    expected = np.sqrt(np.mean((base[3]-base[2].dmr.T)**2)/np.mean(base[3]**2))
    assert base[1]["fit_statistic"] == pytest.approx(expected, rel=2e-5)
    np.testing.assert_allclose(base[2].im, permuted[2].im, rtol=1e-4, atol=1e-8)


def test_signed_samples_are_preserved_including_window_edges(scientific_core, tmp_path):
    phases = np.array([0., .25, .5, .75])
    velocity = np.linspace(-1000, 1000, 41)
    flux = -.2 + 4*np.exp(-.5*((velocity[None, :]-
                             300*np.sin(2*np.pi*phases[:, None]))/160)**2)
    result, status, output, raw, _ = _run(scientific_core, tmp_path/"negative", flux=flux,
                                          alpha=0.0)
    assert result.returncode == 0, result.stdout+result.stderr
    # Gaussian fits keep the requested baseline: negative data are not allowed
    # to lift it (which formerly removed the positivity constraint).
    assert status["baseline"] == pytest.approx(0.1)
    assert status["baseline_raised"] == 0
    assert np.any(output.dm < 0)
    np.testing.assert_allclose(output.dm.T, raw, rtol=5e-5, atol=1e-7)


def test_phase_mask_endpoints_remain_open_across_phase_wrap(scientific_core, tmp_path):
    result, status, _, _, _ = _run(scientific_core, tmp_path/"endpoints",
        phases=np.array([.05, .2, .6, .95]), mask=(.95, 1.05), nim=1)
    assert status, result.stdout+result.stderr
    assert status["statistic_rank"] == 4*41


def test_objective_derivatives_and_projection_adjoint(scientific_core, tmp_path):
    """Check derivatives against independent finite differences, including covariance."""
    output = _compile_with_core(scientific_core, tmp_path, "derivatives", r"""
program derivatives
  use dop_science
  implicit none
  real(dp), allocatable :: pha(:), vp(:), ps(:), direction(:), trial(:), ph(:), pht(:)
  real(dp), allocatable :: dha(:), dsa(:), chi(:), tem(:), testdata(:), back(:), forward(:)
  real(dp) :: h, s, al, q0, derivh, derivs, hessian, eps, hplus, hminus, splus, sminus, fla, fmu
  integer :: ip, iv, j, i
  nph = 4; nvp = 41; nv = 20; nil = nph*nvp
  stat_iw = 2; smoothing_width = 3; likelihood_scale = 0.7_dp
  allocate(pha(nph), dpha(nph), vp(nvp), precision(nvp, nvp, nph), stat_weights(nil))
  allocate(phase_valid(nph))
  phase_valid = .true.; phase_valid(2) = .false.
  stat_weights = 1; precision = 0
  do ip = 1, nph
    pha(ip) = (ip - 1)*pi2/nph
    dpha(ip) = 0.01_dp
    do iv = 1, nvp
      precision(iv, iv, ip) = 2
      if (iv < nvp) then
        precision(iv, iv + 1, ip) = -0.25_dp
        precision(iv + 1, iv, ip) = -0.25_dp
      end if
    end do
  end do
  do iv = 1, nvp
    vp(iv) = -1000 + (iv - 1)*50.0_dp
  end do
  call build_projection(pha, vp)
  allocate(ps(nol), direction(nol), trial(nol), ph(nil), pht(nil), dha(nol), dsa(nol))
  allocate(chi(nol), tem(nol), testdata(nil), back(nol), forward(nil))
  do j = 1, nol
    ps(j) = 2 + 0.5_dp*sin(real(j, dp))
    direction(j) = 0.2_dp*cos(j*0.3_dp)
  end do
  do i = 1, nil
    pht(i) = 1 + 0.2_dp*cos(real(i, dp))
    testdata(i) = sin(i*0.7_dp)
  end do
  call fproj(direction, forward)
  call bproj(testdata, back)
  write(*, *) dot_product(forward, testdata), dot_product(direction, back)
  al = 0.7_dp
  call fproj(ps, ph)
  call like(pht, ph, 1, h, dha)
  call entr(ps, s, dsa, chi, tem)
  derivh = dot_product(direction, dha)
  derivs = dot_product(direction, dsa)
  q0 = h + al*s
  call lamun(pht, ps, direction, direction, ph, chi, tem, al, 1, fla, fmu)
  hessian = -(derivh + al*derivs)/fmu
  eps = 1.0e-3_dp
  trial = ps + eps*direction
  call fproj(trial, ph); call like(pht, ph, 1, hplus, dha); call entr(trial, splus, dsa, chi, tem)
  trial = ps - eps*direction
  call fproj(trial, ph); call like(pht, ph, 1, hminus, dha); call entr(trial, sminus, dsa, chi, tem)
  write(*, *) derivh, (hplus - hminus)/(2*eps)
  write(*, *) derivs, (splus - sminus)/(2*eps)
  write(*, *) hessian, (hplus + al*splus - 2*q0 + hminus + al*sminus)/eps**2
end program derivatives
""")
    pairs = np.array([[float(x) for x in line.split()] for line in output.splitlines()])
    np.testing.assert_allclose(pairs[:, 0], pairs[:, 1], rtol=2e-6, atol=1e-7)


def _projection_profile(binary, tmp_path, phase, width, sigma):
    """One bright map pixel projected through a full exposure, on irregular bins."""
    output = _compile_with_core(binary, tmp_path, "projection", f"""
program projection
  use dop_science
  implicit none
  real(dp), allocatable :: pha(:), vp(:), spot(:), trail(:)
  real(dp) :: x, pi
  integer :: iv, j
  nph = 4; nvp = 41; nv = 20; nil = nph*nvp
  pi = 4*atan(1.0_dp)
  allocate(pha(nph), dpha(nph), vp(nvp))
  pha = {phase!r}_dp*pi2
  dpha = {width!r}_dp*pi2
  do iv = 1, nvp
    x = -1 + 2*real(iv - 1, dp)/(nvp - 1)
    vp(iv) = 1000*x + 100*sin(pi*x)
  end do
  instrument_sigma = {sigma!r}_dp
  phase_subsamples = 1024
  call build_projection(pha, vp)
  allocate(spot(nol), trail(nil))
  spot = 0
  do j = 1, nol
    if (indx(j) == 11 .and. indy(j) == 6) spot(j) = nph
  end do
  call fproj(spot, trail)
  do iv = 1, nvp
    write(*, *) vp(iv), trail(1 + (iv - 1)*nph)
  end do
end program projection
""")
    return np.array([[float(x) for x in line.split()] for line in output.splitlines()]).T


def _independent_square_pixel(centers, phase, width, sigma, vx=0.5, vy=-4.5):
    """Gauss-Legendre average over exposure phase and a uniform square pixel."""
    from scipy.special import ndtr

    edges = np.r_[centers[0]-(centers[1]-centers[0])/2, (centers[:-1]+centers[1:])/2,
                  centers[-1]+(centers[-1]-centers[-2])/2]
    pixel = 2000/19
    phase_nodes, phase_weight = np.polynomial.legendre.leggauss(256)
    pix_nodes, pix_weight = np.polynomial.legendre.leggauss(64)
    angle = 2*np.pi*(phase + width*phase_nodes/2)
    ux = pixel/2*pix_nodes
    # Internal velocity components of pixel (ix=11, iy=6) on a 20-pixel grid,
    # in pixel units: (ix - 10.5, iy - 10.5).
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


def test_exposure_quadrature_lsf_and_irregular_velocity_bins(scientific_core, tmp_path):
    """Compare the forward operator with independent phase/pixel quadrature."""
    centers, actual = _projection_profile(scientific_core, tmp_path, 0.0, 0.5, 40.0)
    expected = _independent_square_pixel(centers, 0.0, 0.5, 40.0)
    np.testing.assert_allclose(actual, expected, atol=5e-7, rtol=2e-5)
    assert actual.sum() == pytest.approx(1.0, abs=1e-7)
    assert centers @ actual > 250  # Endpoint-only smearing formerly gave zero.


def test_square_pixels_project_without_ripple_at_cardinal_phases(scientific_core, tmp_path):
    """A uniform pixel seen at phase 0 projects as a box; the old triangle did not."""
    centers, actual = _projection_profile(scientific_core, tmp_path, 0.0, 1e-12, 0.0)
    edges = np.r_[centers[0]-(centers[1]-centers[0])/2, (centers[:-1]+centers[1:])/2,
                  centers[-1]+(centers[-1]-centers[-2])/2]
    pixel = 2000/19
    lo, hi = 4.0*pixel, 5.0*pixel  # v = -Vy_internal at phase 0, Vy_internal = -4.5 px
    expected = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0, None)/pixel
    np.testing.assert_allclose(actual, expected, atol=1e-9)
    assert actual.sum() == pytest.approx(1.0, abs=1e-12)


def test_alpha_schedule_brackets_target_and_reports_plateau(scientific_core, tmp_path):
    """The schedule refines a crossing by bisection and detects infeasible targets."""
    rng = np.random.default_rng(4)
    phases = np.array([0.0, 0.25, 0.5, 0.75])
    velocity = np.linspace(-1000.0, 1000.0, 41)
    clean = 1 + 4*np.exp(-0.5*((velocity[None, :] - 300*np.sin(2*np.pi*phases[:, None]))/160)**2)
    noisy = clean + 0.3*rng.normal(size=clean.shape)
    result, status, _, _, _ = _run(scientific_core, tmp_path/"bracket", flux=noisy,
                                   alpha=30.0, nal=40, clim=1.1)
    assert result.returncode == 0, result.stdout + result.stderr
    assert status["target_bracketed"] == 1
    assert status["alpha_steps"] > 2
    assert 1.1*(1 - 0.006) <= status["fit_statistic"] <= 1.1
    result, status, _, _, _ = _run(scientific_core, tmp_path/"plateau", flux=noisy,
                                   alpha=30.0, nal=40, clim=0.2)
    assert result.returncode != 0
    assert status["optimizer_status"] == 1
    assert status["rms_plateau"] == 1
