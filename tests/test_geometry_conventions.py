"""Physical-coordinate and geometry regressions independent of display defaults."""

from pathlib import Path

import numpy as np
import pytest

from pydoppler import pydoppler as p


@pytest.fixture
def plt():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as pyplot

    yield pyplot
    pyplot.close("all")


def write_output(path: Path, *, nv=4, va=300.0, phases=(0.25, 0.75),
                 dm=None, dmr=None, im=None, iph=1):
    """Write a small ``dop.out`` with independent physical grids."""
    phases = np.asarray(phases)
    velocities = np.array([-200.0, -100.0, 0.0, 100.0, 200.0])
    nph, nvp = len(phases), len(velocities)
    if dm is None:
        dm = np.ones((nvp, nph))
    if dmr is None:
        dmr = np.asarray(dm)
    if im is None:
        im = np.arange(nv * nv).reshape(nv, nv)
    tokens = [*(str(x) for x in phases * 2 * np.pi), str(iph)]
    if iph == 1:
        tokens += [str(0.1 * 2 * np.pi)] * nph
    tokens += [*(str(x * 1e5) for x in velocities),
               *(str(x) for x in np.ravel(dm)),
               f"0 0 0 0 7 1e-4 1 1 0 1e6 0 {nv} {va * 1e5} dopmap",
               *(str(x) for x in np.ravel(im)),
               f"{nvp * nph} reconstructed data",
               *(str(x) for x in np.ravel(dmr)),
               f"{nv} 4 change maps", *("0" for _ in range(nv * nv))]
    path.write_text(f"{nph} {nvp} {nv} 6562.8 0\n0 input\n"
                    + " ".join(tokens) + "\n", encoding="utf-8")


def test_parser_keeps_spectral_centers_and_independent_map_grid(tmp_path):
    write_output(tmp_path / "dop.out")
    result = p._parse_dopout(tmp_path / "dop.out")
    np.testing.assert_array_equal(result.vp, [-200, -100, 0, 100, 200])
    np.testing.assert_array_equal(result.vmap, [-300, -100, 100, 300])
    assert result.va == 300


def test_parser_supports_fortran_implicit_exposure_widths(tmp_path):
    write_output(tmp_path / "dop.out", iph=0)
    result = p._parse_dopout(tmp_path / "dop.out")
    np.testing.assert_array_equal(result.dpha, [0.5, 0.5])
    np.testing.assert_array_equal(result.vp, [-200, -100, 0, 100, 200])


def test_map_fortran_axis_order_and_physical_pixel_centers(tmp_path, plt):
    image = np.zeros((4, 4))
    # Fortran internal_x=-100, internal_y=+300. Public (Vx,Vy)=(+300,-100).
    image[1, 3] = 9
    write_output(tmp_path / "dop.out", im=image)
    dop = p.spruit(workdir=tmp_path)
    _, returned = dop.Dopmap(plot=True, show=False)
    plotted = plt.gca().images[0]
    assert plotted.origin == "lower"
    assert plotted.get_extent() == pytest.approx([-400, 400, -400, 400])
    row, col = np.unravel_index(np.nanargmax(returned), returned.shape)
    vx, vy = -400 + (col + 0.5) * 200, -400 + (row + 0.5) * 200
    assert (vx, vy) == (300, -100)
    # Fortran projection and public convention agree at independent phases.
    phase = np.array([0, 0.125, 0.25, 0.5]) * 2 * np.pi
    np.testing.assert_allclose(
        -100 * np.sin(phase) - 300 * np.cos(phase),
        -vx * np.cos(phase) + vy * np.sin(phase),
    )


def test_reconstruction_preserves_phase_and_amplitude_errors(tmp_path, plt):
    dm = np.tile([5.0, 1.0], (5, 1))
    write_output(tmp_path / "dop.out", dm=dm, dmr=2 * dm)
    dop = p.spruit(workdir=tmp_path)
    dop.nbins = 10
    dop.Reco(plot=True, colorbar=False, show=False)
    data_panel, model_panel = plt.figure("Reconstruction").axes
    data = data_panel.collections[0]
    model = model_panel.collections[0]
    values = data.get_array().reshape(22, 5)
    np.testing.assert_allclose(model.get_array(), 2 * data.get_array())
    coords = data.get_coordinates()
    centers = (coords[:-1, 0, 1] + coords[1:, 0, 1]) / 2
    i = np.flatnonzero(np.isclose(centers, 0.25))[0]
    assert np.all(values[i] == 1)
    assert data_panel.get_ylim()[0] < data_panel.get_ylim()[1]
    np.testing.assert_array_equal(coords[0, :, 0], [-250, -150, -50, 50, 150, 250])


def test_mean_removal_includes_outer_annulus_and_is_returned(tmp_path):
    image = np.full((4, 4), 4.0)
    center = [1.5, 1.5]
    np.testing.assert_allclose(p.create_profile(image, p.radial_profile(image, center), center), image)
    write_output(tmp_path / "dop.out", im=image)
    _, result = p.spruit(workdir=tmp_path).Dopmap(plot=False, remove_mean=True)
    np.testing.assert_array_equal(result[np.isfinite(result)], 0)
    assert np.isnan(result[0, 0])


def test_smoothing_and_negative_modify_returned_map(tmp_path, plt):
    image = np.zeros((8, 8))
    image[3, 4] = 1
    write_output(tmp_path / "dop.out", nv=8, im=image)
    dop = p.spruit(workdir=tmp_path)
    _, raw = dop.Dopmap(plot=False)
    _, smooth = dop.Dopmap(plot=True, show=False, smooth=True)
    assert not np.array_equal(raw, smooth)
    assert smooth[3, 3] > 0
    np.testing.assert_array_equal(plt.gca().images[0].get_array(), smooth)
    _, negative = dop.Dopmap(plot=False, smooth=True, negative=True)
    np.testing.assert_allclose(negative, -smooth)


def test_phase_overlap_keeps_entire_long_exposure():
    trail, phase = p.rebin_trail([1], [[7]], [0.5], 20, 0.4)
    in_exposure = ((phase % 1) > 0.3) & ((phase % 1) < 0.7)
    np.testing.assert_array_equal(np.isfinite(trail[0]), in_exposure)
    np.testing.assert_allclose(trail[:, in_exposure], 7)


def test_phase_overlap_fractional_weights_and_missing_pixels():
    trail, phase = p.rebin_trail([1, 2], [[2, np.nan], [8, 3]], [0.44, 0.49], 10, 0.08)
    i = np.flatnonzero(np.isclose(phase, 0.45))[0]
    assert trail[0, i] == pytest.approx((2 * 0.08 + 8 * 0.05) / 0.13)
    assert trail[1, i] == 3


def test_phase_overlap_wrap_and_complete_orbits():
    trail, phase = p.rebin_trail([1], [[4]], [0.99], 10, 0.1)
    finite = (phase % 1 < 0.1) | (phase % 1 > 0.9)
    np.testing.assert_array_equal(np.isfinite(trail[0]), finite)
    all_orbits, _ = p.rebin_trail([1], [[9]], [0.31], 10, 2.2)
    np.testing.assert_allclose(all_orbits, 9)


@pytest.mark.parametrize("q", [0.01, 0.35, 1.0, 2.0, 10.0])
def test_roche_l1_equation_and_critical_equipotential(q):
    l1 = p.rlq1(q)
    force = q / (1 - l1) ** 2 - 1 / l1 ** 2 + (1 + q) * l1 - q
    assert abs(force) < 1e-11
    assert l1 + p.rlq1(1 / q) == pytest.approx(1, abs=2e-15)
    r, chi, psi = p.surface(q, l1, 21, 7, None, None, None)
    target = p.pot(q, 1 - l1, 0, 0)[0]
    for i, j in [(2, 1), (5, 3), (10, 6), (20, 0)]:
        x = r[i, j] * np.cos(chi[i])
        y = r[i, j] * np.sin(chi[i]) * np.cos(psi[j])
        z = r[i, j] * np.sin(chi[i]) * np.sin(psi[j])
        assert p.pot(q, x, y, z)[0] == pytest.approx(target, abs=1e-10)


def jacobi(result, q):
    _, _, x, y, inertial, _ = result
    cm = q / (1 + q)
    z = x - cm - 1j * y
    w = inertial - 1j * z
    potential = (1 / abs(z + cm) + q / abs(z - (1 - cm))) / (1 + q) + 0.5 * abs(z) ** 2
    return 2 * potential - abs(w) ** 2


def test_stream_tolerance_convergence_jacobi_and_unequal_grids():
    coarse = p.stream_calculate(0.35, ni=20, nj=5, rtol=1e-7, atol=1e-9)
    fine = p.stream_calculate(0.35, ni=20, nj=5, rtol=1e-10, atol=1e-12)
    assert len(fine[2]) < 1000  # Accepted segment only, not 10,000 padded samples.
    assert len(fine[0]) == 4 * 20
    assert np.ptp(jacobi(fine, 0.35)) < 1e-8
    assert np.ptp(jacobi(fine, 0.35)) < np.ptp(jacobi(coarse, 0.35))
    for i in range(2, 6):
        assert coarse[i][-1] == pytest.approx(fine[i][-1], rel=1e-5, abs=1e-6)
    assert np.rad2deg(np.angle(fine[2][-1] - 1j * fine[3][-1])) == pytest.approx(70, abs=1e-9)


def test_stream_uses_k1_and_warns_about_inconsistent_redundant_values(plt):
    with pytest.warns(RuntimeWarning, match="uses the supplied K1"):
        p.stream(0.35, 20, 0.1769061911, 1.2, 70, title=False)
    lines = plt.gca().lines
    assert lines[-2].get_ydata()[0] == pytest.approx(-20)
    assert lines[-1].get_ydata()[0] == pytest.approx(20 / 0.35)
    first = lines[0].get_xydata().copy()
    plt.close("all")
    with pytest.warns(RuntimeWarning, match="uses the supplied K1"):
        p.stream(0.35, 40, 0.1769061911, 1.2, 70, title=False)
    np.testing.assert_allclose(plt.gca().lines[0].get_xydata(), 2 * first)


@pytest.mark.parametrize("inclination", [30, 70, 90])
def test_resonance_matches_independent_projected_kepler_velocity(plt, inclination):
    q, mass, period = 0.35, 1.2, 0.1769061911 * 86400
    gm = 6.67430e-8 * 1.98847e33 * mass
    omega = 2 * np.pi / period
    separation = (gm * (1 + q) / omega ** 2) ** (1 / 3)
    sini = np.sin(np.deg2rad(inclination))
    k1 = omega * separation * q / (1 + q) * sini / 1e5
    radius = (1 / 3) ** (2 / 3) / (1 + q) ** (1 / 3) * separation
    expected_velocity = np.sqrt(gm / radius) * sini / 1e5
    p.resonance(3, 2, k1, q, period / 86400, mass)
    line = plt.gca().lines[0]
    actual_radii = np.hypot(line.get_xdata(), line.get_ydata() + k1)
    np.testing.assert_allclose(actual_radii, expected_velocity, rtol=1e-12)


@pytest.mark.parametrize("entrypoint", ["Dopmap", "Reco"])
def test_public_output_read_rejects_failed_and_modified_managed_results(tmp_path, entrypoint):
    import hashlib
    import json

    path = tmp_path / "dop.out"
    write_output(path)
    dop = p.spruit(workdir=tmp_path)
    def read():
        return getattr(dop, entrypoint)(plot=False, colorbar=False)

    manifest = {"success": False, "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (tmp_path / "run.json").write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError, match="previous, failed or unverified"):
        read()
    manifest["success"] = True
    (tmp_path / "run.json").write_text(json.dumps(manifest))
    read()
    path.write_text(path.read_text() + "\n")
    with pytest.raises(RuntimeError, match="previous, failed or unverified"):
        read()


def test_phase_overlap_touching_edges_do_not_create_ghost_bins():
    trail, centers = p.rebin_trail([1], [[5]], [0.25], 10, 0.1)
    finite_centers = centers[np.isfinite(trail[0])]
    np.testing.assert_allclose(finite_centers, [0.25, 1.25])


def test_dopmap_display_keeps_zero_flux_at_zero(tmp_path):
    image = np.zeros((4, 4))
    image[1, 1], image[2, 2], image[1, 2] = 4.0, -1.0, 2.0
    write_output(tmp_path / "dop.out", im=image)
    _, shown = p.spruit(workdir=tmp_path).Dopmap(plot=False)
    np.testing.assert_allclose(shown[1, 1], 1.0)
    np.testing.assert_allclose(shown[2, 2], -0.25)  # negative stays negative
    np.testing.assert_allclose(shown[1, 2], 0.5)
    assert shown[2, 1] == 0.0  # zero flux displays as zero


def test_resolved_schedule_defaults_and_validation():
    from pydoppler.runtime import resolved_schedule

    dop = p.spruit()
    assert resolved_schedule(dop, 0) == (0.0, 0, 1.0)
    assert resolved_schedule(dop, 2) == (0.0, 25, 1.0)
    dop.al0, dop.nal, dop.clim = 3.0, 4, 1.2
    assert resolved_schedule(dop, 1) == (3.0, 4, 1.2)
    dop.al0 = 0.0
    with pytest.raises(ValueError, match="al0 must be positive"):
        resolved_schedule(dop, 0)


def test_obsolete_projection_factors_warn(tmp_path):
    dop = p.spruit(workdir=tmp_path)
    with pytest.warns(DeprecationWarning, match="obsolete"):
        with pytest.raises(FileNotFoundError):
            dop.Syncdop(nri=0.9)


def test_packaged_sample_script_matches_repository_copy():
    root = Path(__file__).resolve().parents[1]
    assert ((root / "sample_script.py").read_text()
            == (root / "pydoppler" / "test_data" / "sample_script.py").read_text())
