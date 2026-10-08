import shutil
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pytest

import pydoppler
from pydoppler import pydoppler as pydoppler_module


def _expected_overlap(data_dir: Path) -> Tuple[float, float]:
    mins = []
    maxs = []
    for path in sorted(data_dir.glob("txhugem4*")):
        arr = np.loadtxt(path)
        mins.append(float(arr[0, 0]))
        maxs.append(float(arr[-1, 0]))
    return max(mins), min(maxs)


def _write_mock_dopout(
    path: Path,
    im: np.ndarray,
    dm: Optional[np.ndarray] = None,
    dmr: Optional[np.ndarray] = None,
    vp: Optional[np.ndarray] = None,
) -> None:
    im = np.asarray(im, dtype=float)
    nv = int(im.shape[0])
    if im.shape != (nv, nv):
        raise ValueError("im must be square.")

    nph = 2
    nvp = 3
    if dm is None:
        dm = np.zeros((nvp, nph), dtype=float)
    if dmr is None:
        dmr = np.zeros((nvp, nph), dtype=float)
    dpx = np.zeros((nv, nv), dtype=float)

    pha = np.array([0.0, np.pi], dtype=float)
    dpha = np.array([0.1, 0.1], dtype=float) * 2.0 * np.pi
    if vp is None:
        vp = np.array([1.0e5, 2.0e5, 3.0e5], dtype=float)
    vp = np.asarray(vp, dtype=float)
    params = f"0 0 0.0 1.0 7 1e-4 1.0 1.0 1 1.0 0.0 {nv} {np.max(np.abs(vp))} dopmap"

    payload = " ".join(
        [
            *(f"{value:.8e}" for value in pha),
            "1",
            *(f"{value:.8e}" for value in dpha),
            *(f"{value:.8e}" for value in vp),
            *(f"{value:.8e}" for value in np.asarray(dm, dtype=float).ravel()),
            params,
            *(f"{value:.8e}" for value in im.ravel()),
            "0 0 0",
            *(f"{value:.8e}" for value in np.asarray(dmr, dtype=float).ravel()),
            "0 0 0 0",
            *(f"{value:.8e}" for value in dpx.ravel()),
        ]
    )

    path.write_text(
        "\n".join(
            [
                f"{nph} {nvp} {nv} 6562.8 0.0",
                "0.0 0 0 0",
                payload,
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_foldspec_uses_common_overlap_for_bundled_dataset(tmp_path: Path):
    workdir = tmp_path / "workdir"
    pydoppler.copy_test_data(workdir)

    data_dir = workdir / "ugem99"
    expected_min, expected_max = _expected_overlap(data_dir)

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=workdir)
    dop.base_dir = data_dir
    dop.list = "ugem0all.fas"
    dop.Foldspec()

    assert dop.wave
    assert dop.wave[0][0] >= expected_min - 1e-8
    assert dop.wave[0][-1] <= expected_max + 1e-8
    assert dop.trsp.shape[0] == dop.input_phase.size
    assert dop.trsp.shape[1] == dop.wave[0].size


@pytest.mark.skipif(
    shutil.which("make") is None or shutil.which("gfortran") is None,
    reason="Fortran toolchain is required for the bundled workflow regression.",
)
def test_bundled_dataset_runs_full_pipeline_headless(tmp_path: Path):
    workdir = tmp_path / "workdir"
    pydoppler.copy_fortran_code(workdir)
    pydoppler.copy_test_data(workdir)

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=workdir)
    dop.base_dir = workdir / "ugem99"
    dop.list = "ugem0all.fas"
    dop.lam0 = 6562.8
    dop.delta_phase = 0.003
    dop.delw = 35
    dop.overs = 0.3
    dop.gama = 36.0
    dop.nbins = 28

    dop.Foldspec()
    dop.Dopin(plot=False, continuum_band=[6500, 6537, 6591, 6620])
    dop.Syncdop()
    _, dopmap = dop.Dopmap(plot=False)

    assert (workdir / "dop.out").is_file()
    assert dopmap.ndim == 2
    assert np.any(np.isfinite(dopmap))


def test_dopmap_handles_flat_zero_map(tmp_path: Path):
    workdir = tmp_path
    _write_mock_dopout(workdir / "dop.out", im=np.zeros((4, 4), dtype=float))

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=workdir)
    _, data = dop.Dopmap(plot=False)

    assert data.shape == (4, 4)
    assert np.allclose(data[np.isfinite(data)], 0.0)
    assert np.isnan(data[0, 0])


def test_dopmap_plot_extent_uses_map_grid_not_spectral_grid(tmp_path: Path):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    workdir = tmp_path
    _write_mock_dopout(
        workdir / "dop.out",
        im=np.ones((4, 4), dtype=float),
        vp=np.array([1.0e5, 4.0e5, 7.0e5], dtype=float),
    )

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=workdir)
    dop.Dopmap(plot=True, colorbar=False, show=False)

    extent = plt.gca().images[0].get_extent()

    assert extent == pytest.approx([-28 / 3, 28 / 3, -28 / 3, 28 / 3])
    assert plt.gca().images[0].get_cmap().name == "magma"
    plt.close("all")


def test_dopmap_interactive_colorbar_sets_norm_on_mappable(tmp_path: Path):
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib.colors import Normalize
    import matplotlib.pyplot as plt

    workdir = tmp_path
    _write_mock_dopout(
        workdir / "dop.out",
        im=np.arange(1, 17, dtype=float).reshape(4, 4),
    )

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=workdir)
    cbar, data = dop.Dopmap(
        plot=True,
        colorbar=True,
        show=False,
        limits=[1e-3, 1.0],
    )

    assert data.shape == (4, 4)
    assert cbar.mappable.norm is cbar.cbar.norm
    assert type(cbar.cbar.norm) is Normalize
    assert cbar.cbar.norm.vmin == pytest.approx(1e-3)
    assert cbar.cbar.norm.vmax == pytest.approx(1.0)
    plt.close("all")


def test_reco_colorbars_use_decimal_labels_and_full_panel_height(tmp_path: Path):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    dm = np.array([[0.1, 0.2], [0.4, 0.5], [0.8, 0.9]], dtype=float)
    dmr = np.array([[0.2, 0.1], [0.5, 0.4], [0.9, 0.8]], dtype=float)
    _write_mock_dopout(tmp_path / "dop.out", im=np.ones((4, 4), dtype=float), dm=dm, dmr=dmr)

    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=tmp_path)
    cbar_data, cbar_model, _, _ = dop.Reco(
        plot=True,
        colorbar=True,
        show=False,
        limits=[0.05, 0.95],
    )

    fig = plt.figure("Reconstruction")
    fig.canvas.draw()
    for wrapped_cbar in (cbar_data, cbar_model):
        tick_labels = [tick.get_text() for tick in wrapped_cbar.cbar.ax.get_yticklabels()]
        assert tick_labels
        assert all("e" not in label.lower() for label in tick_labels)
        assert all("\\times" not in label for label in tick_labels)

        panel_box = wrapped_cbar.mappable.axes.get_position()
        cbar_box = wrapped_cbar.cbar.ax.get_position()
        assert cbar_box.y0 == pytest.approx(panel_box.y0, abs=1e-2)
        assert cbar_box.y1 == pytest.approx(panel_box.y1, abs=1e-2)

    plt.close("all")


def test_scale_by_absmax_handles_zero_input():
    scaled = pydoppler_module._scale_by_absmax(np.zeros((4, 4), dtype=float))

    assert scaled.shape == (4, 4)
    assert np.allclose(scaled, 0.0)
