import sys
from pathlib import Path

import numpy as np

import pydoppler


def _make_synthetic_doppler(tmp_path: Path, nspec: int = 6):
    dop = pydoppler.spruit(auto_install=False, interactive=False, workdir=tmp_path)
    wave = np.linspace(6500.0, 6600.0, 120)
    center = 6550.0
    profile = np.exp(-0.5 * ((wave - center) / 3.0) ** 2)
    rng = np.random.default_rng(123)
    flux = [10.0 + profile + 0.1 * rng.normal(size=wave.size) for _ in range(nspec)]

    dop.wave = [wave]
    dop.flux = flux
    phases = np.linspace(0.0, 1.0, nspec, endpoint=False)
    dop.pha = phases
    dop.input_phase = phases
    dop.delta_phase = 0.01
    dop.lam0 = center
    dop.delw = 25
    dop.gama = 0.0
    return dop


def test_dopin_plot_false_does_not_import_matplotlib(tmp_path: Path):
    sys.modules.pop("matplotlib", None)

    dop = _make_synthetic_doppler(tmp_path)
    dop.base_dir = tmp_path
    dop.Dopin(plot=False)

    assert "matplotlib" not in sys.modules
    assert (tmp_path / "dopin").is_file()


def test_dopin_writes_error_bars_when_iw_enabled(tmp_path: Path):
    dop = _make_synthetic_doppler(tmp_path)
    dop.iw = 1
    dop.Dopin(plot=False)

    text = (tmp_path / "dopin").read_text(encoding="utf-8")
    tokens = text.split()
    numeric = []
    for token in tokens:
        try:
            numeric.append(float(token))
        except ValueError:
            continue

    nph = int(numeric[0])
    nvp = int(numeric[1])

    expected = (
        3  # header: nph nvp lam0
        + 3  # gamma + 0 + 0 (the filename token is non-numeric)
        + nph  # phases
        + 1  # iph
        + nph  # delta phases
        + nvp  # velocity grid
        + (nvp * nph)  # flux values
        + (nvp * nph)  # error bars
    )
    assert len(numeric) == expected
