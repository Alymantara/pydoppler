"""Independent statistical and optimizer checks of the compiled scientific core."""
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from pydoppler.pydoppler import _parse_dopout


@pytest.fixture(scope="module")
def scientific_core(tmp_path_factory):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the independent core checks")
    work = tmp_path_factory.mktemp("fortran_science_build")
    source = Path(__file__).resolve().parents[1] / "pydoppler" / "fortran_code"
    (work / "emap.par").write_text(
        "      parameter (npm=4,nvpm=41,nvm=20)\n"
        "      parameter (nd=npm*nvpm,nr=nvm*nvm)\n"
        "      parameter (nri=nr,ndi=nd)\n"
        "      parameter (nf=nd,nfp=ndi,nrp=nri)\n"
        "      parameter (ni=nvm,nj=nvm)\n"
        "      parameter (nch=1,nsu=1)\n"
    )
    shutil.copy(source / "dop.f", work / "dop.f")
    shutil.copy(source / "cclock.f", work / "clock.f")
    result = subprocess.run(
        [compiler, "-O0", "-g", "-fcheck=all", "-ffixed-line-length-none",
         "dop.f", "clock.f", "-o", "dopp"],
        cwd=work, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Warning:" not in result.stderr
    return work / "dopp"


def _run(binary, work, *, scale=1.0, iw=1, precision=None, norm=0,
         phases=None, flux=None, nim=1000, ac=8e-4, alpha=.003, mask=(0.0, 0.0)):
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
        f"0\n{iw}\n{mask[0]} {mask[1]}\n3\n{ac}\n{nim}\n"
        f"{alpha} 1.7 -1\n99\n1\n{norm}\n1 0\n0 0 .1\n"
    )
    if iw == 2:
        if precision is None:
            precision = np.repeat((np.eye(41)/0.3**2)[None], 4, axis=0)
        with (work / "dop.precision").open("w") as handle:
            handle.write("4 41\n")
            for matrix in precision:
                handle.write("41\n")
                np.savetxt(handle, matrix/scale**2, fmt="%.17e")
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
    result, _, output, raw, _ = _run(scientific_core, tmp_path/"negative", flux=flux)
    assert result.returncode == 0, result.stdout+result.stderr
    assert np.any(output.dm < 0)
    np.testing.assert_allclose(output.dm.T, raw, rtol=5e-5, atol=1e-7)


def test_phase_mask_endpoints_remain_open_across_phase_wrap(scientific_core, tmp_path):
    result, status, _, _, _ = _run(scientific_core, tmp_path/"endpoints",
        phases=np.array([.05, .2, .6, .95]), mask=(.95, 1.05), nim=1)
    assert status, result.stdout+result.stderr
    assert status["statistic_rank"] == 4*41


def test_objective_derivatives_and_projection_adjoint(scientific_core, tmp_path):
    """Check derivatives against independent finite differences, including covariance."""
    build = scientific_core.parent
    source = (build / "dop.f").read_text()
    core = source[:source.index("      program dopmap")]+source[source.index("      subroutine statistical_input"):]
    (tmp_path / "core.f").write_text(core)
    shutil.copy(build / "emap.par", tmp_path / "emap.par")
    shutil.copy(build / "clock.f", tmp_path / "clock.f")
    (tmp_path / "derivatives.f").write_text(r'''
      program derivatives
      use dop_science
      implicit real*8(a-h,o-z)
      include 'emap.par'
      common/phases/dpha(npm)
      common/foldn/nv,indx(nr),indy(nr),nol
      common/cabsf/wid,af
      real*4 p(ndi,nr),pt(nri,nd)
      dimension pha(npm),vp(nvpm),ps(nr),direction(nr),trial(nr)
      dimension ph(nd),pht(nd),dha(nr),dsa(nr),chi(nr,1),tem(nr),ex(nr)
      dimension epc(1),testdata(nd),back(nr),forward(nd)
      integer inc(1,1),nsm(1,1),ico(1)
      stat_iw=2
      nphase=npm
      nvelocity=nvpm
      allocate(precision(nvpm,nvpm,npm),stat_weights(nd),phase_valid(npm))
      phase_valid=.true.
      phase_valid(2)=.false.
      stat_weights=1.d0
      precision=0.d0
      do ip=1,npm
        pha(ip)=(ip-1)*8.d0*atan(1.d0)/npm
        dpha(ip)=0.01d0
        do iv=1,nvpm
          precision(iv,iv,ip)=2.d0
          if(iv.lt.nvpm) then
            precision(iv,iv+1,ip)=-0.25d0
            precision(iv+1,iv,ip)=-0.25d0
          endif
        enddo
      enddo
      do iv=1,nvpm
        vp(iv)=-1000.d0+(iv-1)*50.d0
      enddo
      nv=nvm
      wid=1.d0
      af=0.d0
      call prom(va,nv,pha,npm,vp,nvpm,pt,p)
      do j=1,nol
        ps(j)=2.d0+0.5d0*sin(dble(j))
        direction(j)=0.2d0*cos(dble(j)*0.3d0)
      enddo
      do i=1,nd
        pht(i)=1.d0+0.2d0*cos(dble(i))
        testdata(i)=sin(dble(i)*0.7d0)
      enddo
      call fproj(pt,nd,nol,direction,forward)
      call bproj(p,nd,nol,testdata,back)
      write(*,*)dot_product(forward,testdata),
     *          dot_product(direction(1:nol),back(1:nol))
      inc=1
      nsm=3
      ico=1
      epc=1.d0
      al=0.7d0
      call fproj(pt,nd,nol,ps,ph)
      call like(pht,ph,p,nd,nol,1,h,dha)
      call entr(ps,nd,nol,s,dsa,chi,tem,ex,1,inc,nsm,epc,ico,0.d0,1,0.d0)
      derivh=dot_product(direction(1:nol),dha(1:nol))
      derivs=dot_product(direction(1:nol),dsa(1:nol))
      q0=h+al*s
      call lamun(pht,ps,direction,direction,pt,ph,chi,tem,ex,al,nd,nol,
     *           fla,fmu,1,inc,nsm,epc,ico,0.d0,1,1,0.d0)
      hessian=-(derivh+al*derivs)/fmu
      eps=1.d-3
      trial(1:nol)=ps(1:nol)+eps*direction(1:nol)
      call fproj(pt,nd,nol,trial,ph)
      call like(pht,ph,p,nd,nol,1,hplus,dha)
      call entr(trial,nd,nol,splus,dsa,chi,tem,ex,1,inc,nsm,epc,ico,0.d0,1,0.d0)
      trial(1:nol)=ps(1:nol)-eps*direction(1:nol)
      call fproj(pt,nd,nol,trial,ph)
      call like(pht,ph,p,nd,nol,1,hminus,dha)
      call entr(trial,nd,nol,sminus,dsa,chi,tem,ex,1,inc,nsm,epc,ico,0.d0,1,0.d0)
      write(*,*)derivh,(hplus-hminus)/(2.d0*eps)
      write(*,*)derivs,(splus-sminus)/(2.d0*eps)
      write(*,*)hessian,(hplus+al*splus-2.d0*q0+hminus+al*sminus)/eps**2
      end
''')
    result = subprocess.run(
        [shutil.which("gfortran"), "-O0", "-fcheck=all", "-ffixed-line-length-none",
         "core.f", "derivatives.f", "clock.f", "-o", "derivatives"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    result = subprocess.run(["./derivatives"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    pairs = np.array([[float(x) for x in line.split()] for line in result.stdout.splitlines()])
    np.testing.assert_allclose(pairs[:, 0], pairs[:, 1], rtol=2e-6, atol=1e-7)


def test_exposure_quadrature_lsf_and_irregular_velocity_bins(scientific_core, tmp_path):
    """Compare the forward operator with independent phase/pixel quadrature."""
    from scipy.special import ndtr

    build = scientific_core.parent
    source = (build / "dop.f").read_text()
    core = source[:source.index("      program dopmap")]+source[source.index("      subroutine statistical_input"):]
    (tmp_path / "core.f").write_text(core)
    shutil.copy(build / "emap.par", tmp_path / "emap.par")
    shutil.copy(build / "clock.f", tmp_path / "clock.f")
    (tmp_path / "projection.f").write_text(r'''
      program projection
      use dop_science
      implicit real*8(a-h,o-z)
      include 'emap.par'
      common/phases/dpha(npm)
      common/foldn/nv,indx(nr),indy(nr),nol
      common/cabsf/wid,af
      real*4 p(ndi,nr),pt(nri,nd)
      dimension pha(npm),vp(nvpm),spot(nr),trail(nd)
      pi=4.d0*atan(1.d0)
      pha=0.d0
      dpha=pi
      do iv=1,nvpm
        x=-1.d0+2.d0*(iv-1)/(nvpm-1)
        vp(iv)=1000.d0*x+100.d0*sin(pi*x)
      enddo
      nv=nvm
      wid=1.d0
      af=0.d0
      instrument_sigma=40.d0
      phase_subsamples=1024
      call prom(va,nv,pha,npm,vp,nvpm,pt,p)
      spot=0.d0
      do j=1,nol
        if(indx(j).eq.11.and.indy(j).eq.6) spot(j)=npm
      enddo
      call fproj(pt,nd,nol,spot,trail)
      do iv=1,nvpm
        write(*,*)vp(iv),trail(1+(iv-1)*npm)
      enddo
      end
''')
    result = subprocess.run(
        [shutil.which("gfortran"), "-O2", "-fcheck=all", "-ffixed-line-length-none",
         "core.f", "projection.f", "clock.f", "-o", "projection"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    result = subprocess.run(["./projection"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    profile = np.array([[float(x) for x in line.split()] for line in result.stdout.splitlines()])
    centers, actual = profile.T
    edges = np.r_[centers[0]-(centers[1]-centers[0])/2,
                  (centers[:-1]+centers[1:])/2,
                  centers[-1]+(centers[-1]-centers[-2])/2]
    # Independent Gaussian bin integrals, averaged over orbital phase and
    # the two halves of a triangular pixel basis with Gauss-Legendre nodes.
    phase_nodes, phase_weight = np.polynomial.legendre.leggauss(256)
    pixel_nodes, pixel_weight = np.polynomial.legendre.leggauss(48)
    pixel_halfwidth = (2000/19)/2
    pixel_left = pixel_halfwidth*(pixel_nodes-1)/2
    pixel_right = pixel_halfwidth*(pixel_nodes+1)/2
    pixel = np.r_[pixel_left, pixel_right]
    pixel_weight = np.tile(pixel_weight, 2)*pixel_halfwidth/2
    pixel_weight *= (pixel_halfwidth-np.abs(pixel))/pixel_halfwidth**2
    phase = phase_nodes*np.pi/2
    radial = (2000/19)*(.5*np.sin(phase)+4.5*np.cos(phase))
    cdf = ndtr((edges[:, None, None]-radial[None, :, None]-pixel[None, None, :])/40)
    expected = np.einsum("vpu,p,u->v", np.diff(cdf, axis=0), phase_weight/2, pixel_weight)
    np.testing.assert_allclose(actual, expected, atol=5e-7, rtol=2e-5)
    assert actual.sum() == pytest.approx(1.0, abs=1e-7)
    assert centers @ actual > 250  # Endpoint-only smearing formerly gave zero.
