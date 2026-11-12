# PyDoppler — Native Python

Pure NumPy/SciPy port of Spruit's Doppler tomography code.
No Fortran compiler required.

## Structure

```
pydoppler/  (repository root)
├── pydoppler/
│   ├── __init__.py        # package entry point
│   ├── pydoppler.py       # spruit class + stream() + Roche lobe physics
│   ├── mem_engine.py      # MEM tomography engine (replaces dop.f)
│   ├── utils.py           # test_data() helper
│   └── test_data/
│       ├── sample_script.py   # copied to cwd by test_data()
│       └── ugem99/            # U Gem Hα spectra (28 spectra)
├── notebooks/PyDoppler_tutorial.ipynb
├── docs/images/           # README figures
├── setup.py / pyproject.toml / MANIFEST.in
└── README.md, CHANGELOG.md, LICENSE
```

## Install

```bash
pip install -e .
```

Requires: numpy, scipy, matplotlib

## Quick start

```python
import pydoppler, matplotlib.pyplot as plt

pydoppler.test_data()           # copies ugem99/ here

dop = pydoppler.spruit()
dop.object      = 'U Gem'
dop.base_dir    = 'ugem99'
dop.list        = 'ugem0all.fas'
dop.lam0        = 6562.8
dop.delta_phase = 0.003
dop.delw        = 35
dop.gama        = 36.0
dop.nbins       = 28

dop.Foldspec()
dop.Dopin(continnum_band=[6500, 6537, 6591, 6620])
dop.Syncdop()                   # pure Python MEM — no gfortran needed

fig, data = dop.Dopmap(colorbar=True)   # linear stretch, magma_r
pydoppler.stream(0.35, 107, 0.1769061911, 1.2, 70,
                 tidal_radius=True, resonance_32=True, resonance_21=True)
plt.show()
```

## Key parameters (spruit class)

| Param | Default | Meaning |
|---|---|---|
| `nim` | 150 | Max RL iterations per alpha |
| `nal` | 0 | Alpha steps (0 = single step; 10+ for production) |
| `clim` | 1.6 | Target reduced chi² |
| `al0` | 0.002 | Starting regularisation alpha |
| `ns` | 7 | Gaussian default-map smoothing width |
| `overs` | 0.3 | Map oversampling (controls nv auto-size) |

Pass `nv=` to `Syncdop()` to override map grid size (e.g. `nv=100`).

## Dopmap options

```python
dop.Dopmap(
    stretch='linear',     # default; or 'sqrt', 'log'
    percentile_lo=5,      # clip low end at this percentile
    percentile_hi=99,     # clip high end at this percentile
    cmaps=plt.cm.magma_r, # default (pale bg); afmhot for dark bg
    colorbar=True,
    smooth=False,
)
```

## stream() options

```python
pydoppler.stream(q, k1, porb, m1, inc,
    colors='black',         # default; use 'white' on dark cmaps
    tidal_radius=True,      # Paczyński tidal truncation radius
    resonance_32=True,      # 3:2 mean-motion resonance (superhumps)
    resonance_21=True,      # 2:1 resonance
    both_lobes=False,       # full donor Roche lobe
)
```

## Notes

- Version 1.0.0. The projection matrix is built vectorised and sparse: the
  full U Gem sample script runs in ~5 s (Syncdop ~0.7 s at nv=144).
- `mem_engine.py` follows Spruit's `dop.f` (prom, chimg, entr, lamun, basel,
  alfs); see CHANGELOG.md for the fixes relative to the first port.
- The original Fortran `dop.f` / `dopp` binary is no longer needed or shipped.
- Tutorial: `notebooks/PyDoppler_tutorial.ipynb`; README images in `docs/images/`.
- Developed with the help of Claude (Anthropic), using Claude Opus 5.5
  (model `claude-opus-5-5`) in Claude Code.
