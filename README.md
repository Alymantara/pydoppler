# PyDoppler

PyDoppler is a Python implementation of Henk Spruit's floating-default
maximum-entropy Doppler tomography. It turns phase-resolved spectra of a
single emission line into trailed spectra, a Doppler tomogram and a
reconstructed trail, with input validation, noise propagation and run
provenance.

The solver descends from Spruit's 1999 Fortran code and has since been
revised:

- **Forward model:** an exact square-pixel projection, with exposure smearing and an optional Gaussian line-spread function.
- **Likelihoods:** Gaussian likelihoods with diagonal or full noise covariance.
- **Regularization:** an alpha schedule that brackets the RMS target and refines it by bisection.
- **Implementation:** NumPy and SciPy only, so no compiler is needed. It reproduces the earlier Fortran solver to rounding accuracy.

<p align="middle">
   <img src="https://raw.githubusercontent.com/fcotizelati/pydoppler/master/docs/readme_images/Doppler_Map.svg" alt="Doppler tomogram of U Gem with Roche-lobe and stream overlays" width="520" />
</p>

## Contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Input data](#input-data)
4. [Tutorial: U Gem](#tutorial-u-gem)
5. [Noise-weighted fits and regularization](#noise-weighted-fits-and-regularization)
6. [Parameters](#parameters)
7. [Conventions in brief](#conventions-in-brief)
8. [Outputs and reproducibility](#outputs-and-reproducibility)
9. [Scope and limitations](#scope-and-limitations)
10. [Citation](#citation)

## Installation

Requirements:

- Python 3.9 or newer.
- NumPy, SciPy and Matplotlib.

Install this version from GitHub:

```bash
pip install "git+https://github.com/fcotizelati/pydoppler"
```

`pip install pydoppler` installs the release on PyPI, which may predate the
solver and conventions described here.

## Quick start

```python
from pathlib import Path

import pydoppler

workdir = Path.cwd() / "pydoppler-workdir"   # solver input and all outputs
pydoppler.copy_test_data(workdir)             # bundled U Gem spectra

dop = pydoppler.spruit(workdir=workdir)
dop.base_dir = workdir / "ugem99"
dop.list = "ugem0all.fas"
dop.lam0, dop.delw, dop.gama = 6562.8, 35, 36.0   # Angstrom, Angstrom, km/s
dop.overs = 0.3

dop.Foldspec()                                     # read spectra and phases
dop.Dopin(continuum_band=[6500, 6537, 6591, 6620],
          use_list_dpha=True, plot=False)          # continuum-subtract, write input
dop.Syncdop()                                      # run the inversion, validate
_, tomogram = dop.Dopmap(plot=False)               # map / peak; NaN outside support
```

The complete tutorial with all figures is in `sample_script.py`. Run it with
`python sample_script.py`, or `%run sample_script.py` in IPython. The packaged
copy installs with `pydoppler.install_sample_script(".")`.

## Input data

**Spectra.** One whitespace-separated text file per spectrum, with columns
`wavelength flux [flux_error]`. Lines starting with `#` are ignored.

- Wavelengths are in Angstrom.
- Fluxes are samples of a spectral density, not integrated pixel counts.
- The optional third column holds 1σ errors and must be strictly positive.

**Phase file.** One line per spectrum: `file_name phase [exposure_width]`.

- Phases are orbital cycles at **mid-exposure**.
- Exposure widths are full widths in cycles. Use them with `Dopin(use_list_dpha=True)`; otherwise the scalar `dop.delta_phase` applies to every spectrum.

```text
txhugem4004  0.7150 0.03930
txhugem4005  0.7794 0.03930
...
```

Keep the spectra and the phase file in one directory (`dop.base_dir`). All
solver files go to `dop.workdir`. Use `dop.make_run_dir()` to start a fresh run
directory for each map.

## Tutorial: U Gem

The bundled data are 28 H-alpha spectra of U Gem in quiescence, kindly
provided by J. Echevarria and published in
[Echevarria et al. 2007, AJ, 134, 262](https://ui.adsabs.harvard.edu/abs/2007AJ....134..262E/abstract).
Each has an effective 600 s exposure, which is 0.0393 cycles. Copy them with
`pydoppler.test_data()` or `pydoppler.copy_test_data(path)`.

### 1. Read the spectra

```python
dop = pydoppler.spruit(workdir="pydoppler-workdir")
dop.object = "U Gem"
dop.base_dir = "ugem99"
dop.list = "ugem0all.fas"
dop.lam0 = 6562.8   # rest wavelength, Angstrom
dop.delw = 35       # half-width of the rest-frame window, Angstrom
dop.gama = 36.0     # systemic velocity, km/s
dop.overs = 0.3     # map pixels per axis / spectral channels
dop.nbins = 28      # phase bins for diagnostic plots only
dop.Foldspec()
```

`Foldspec()` trims all spectra to their common wavelength range and keeps the
native samples, so later resampling happens once.

### 2. Subtract the continuum

```python
dop.Dopin(continuum_band=[6500, 6537, 6591, 6620], poly_degree=2,
          use_list_dpha=True)
```

`continuum_band` gives two continuum windows in the **observed** wavelength
frame.

- The fitted polynomial is subtracted from the spectrum, not divided into it.
- The spectra are resampled onto a uniform velocity grid around `lam0 ± delw`.
- Without `continuum_band`, the outer 5% of the coverage is used. Check that those windows are free of other lines.

<p align="middle">
   <img src="https://raw.githubusercontent.com/fcotizelati/pydoppler/master/docs/readme_images/Trail.svg" alt="Average profile and trailed spectra" width="400" />
</p>

### 3. Run the inversion

```python
diagnostics = dop.Syncdop()
```

Every call runs the inversion. A run counts as successful only if both of
these hold:

- the optimizer converged;
- the RMS statistic met the target `clim`.

Anything else raises an exception, so failed runs can never be plotted. The
bundled spectra have no errors, so this fit is unweighted and its statistic is
the relative RMS residual. See
[noise-weighted fits](#noise-weighted-fits-and-regularization) for fits with errors.

### 4. Plot the tomogram

```python
import matplotlib.pyplot as plt

cbar, tomogram = dop.Dopmap(colorbar=True, cmaps=plt.cm.magma)
pydoppler.stream(q=0.35, k1=107, porb=0.1769061911, m1=1.2, inc=70)
```

How the map is displayed:

- It is divided by its positive peak, so zero flux is zero on the colour bar.
- The default colour limits are 0 and the 99.5th percentile.
- Optional processing: `smooth`, `remove_mean` and `negative`.
- The returned array is exactly what is displayed.

`stream()` overlays the Roche lobes, the ballistic stream and the Keplerian
velocities along it. `k1` is the projected primary semi-amplitude and sets the
velocity scale. The masses, period and inclination are checked against it, with
a warning if they disagree.

### 5. Compare data and reconstruction

```python
cb_data, cb_model, model, data = dop.Reco(colorbar=True, cmaps=plt.cm.magma)
```

`data` and `model` have shape `(n_velocity, n_spectra)` at the original
phases; only the plots are binned in phase. Both panels share one flux scale,
so amplitude errors are visible. Any map feature should be traceable in the
trailed spectra.

<p align="middle">
   <img src="https://raw.githubusercontent.com/fcotizelati/pydoppler/master/docs/readme_images/Reconstruction.svg" alt="Input and reconstructed trailed spectra" width="680" />
</p>

### 6. Residuals

```python
residuals, standardized = dop.Residuals()
# Full-covariance fits only: residuals whitened by the propagated covariance
residuals, whitened = dop.Residuals(whiten=True)
```

`standardized` is `None` when no measured or estimated errors exist. Spectra
excluded by `pb0, pb1` are NaN.

## Noise-weighted fits and regularization

For a statistically interpretable fit, set `dop.iw = 1` **before** `Dopin()`.

- **Error source:** errors come from the spectrum files' third column. Without them, a constant noise level per spectrum is estimated from the continuum-window residuals.
- **Covariance:** interpolation and continuum subtraction correlate the prepared channels. `Dopin()` propagates their full covariance, and with the default `covariance_mode="full"` the solver fits it exactly. `covariance_mode="diagonal"` ignores the correlations, overstates the information in the data and tends to over-fit noise.
- **Statistic:** the fit statistic is the RMS standardized residual. It is not a reduced chi-square.

```python
dop.iw = 1
dop.Dopin(continuum_band=[6500, 6537, 6591, 6620], poly_degree=2,
          use_list_dpha=True)
dop.clim = 1.2          # target RMS of standardized residuals (default 1.0)
dop.Syncdop()
```

**Choosing regularization.** The entropy weight `alpha` trades smoothness
against fit.

- With `al0=None` (default), the solver starts from the historical weight 0.002 per unit of channel information. That is exactly 0.002 in unweighted fits and the equivalent value in weighted fits.
- Unweighted fits use that single alpha.
- Weighted fits decrease alpha by `alf` for up to `nal` steps until the RMS statistic reaches `clim`, then bisect in log(alpha) to land within 0.5% below the target. The result is close to the most regularized map consistent with the target.
- If the statistic levels off above `clim`, `Syncdop()` raises an error that reports the attainable floor.

**When the target is out of reach.** For the U Gem data with
continuum-estimated errors, the floor is about 1.07, so `clim=1.0` is not
attainable. Possible reasons:

- the noise estimate is somewhat low;
- the emission may vary in ways a stationary positive map cannot describe;
- the continuum may be misplaced.

Raise `clim` only for a documented reason. Lower targets fit noise: they give
maps with more pixel-scale structure.

**Positivity.** The solver adds a flat baseline map so that a positive latent
map can describe signed, continuum-subtracted data. In least-squares and
Gaussian fits the baseline stays at `baseline_fraction`, 10% of the RMS flux.
As a result, returned maps can dip only slightly below zero. The legacy log
objective `ih=0` cannot handle negative samples. It therefore raises the
baseline as needed, which effectively removes positivity, and `Syncdop()` logs
a warning when that happens.

## Parameters

Set these attributes on the `spruit` object before `Dopin()` (data
parameters) or `Syncdop()` (solver parameters).

| Attribute | Default | Meaning |
|---|---|---|
| `lam0` | 6562.83 | Rest wavelength (Å), same air/vacuum convention as the data |
| `delw` | 80 | Half-width of the rest-frame window (Å): `lam0 ± delw` |
| `gama` | 0 | Systemic velocity (km/s); positive removes a redshift |
| `delta_phase` | 0.001 | Full exposure width (cycles) when `use_list_dpha=False` |
| `iw` | 0 | 0 unweighted; 1 Gaussian fit using errors (set before `Dopin`) |
| `covariance_mode` | `"full"` | `"full"` propagated covariance or `"diagonal"` approximation |
| `ih` | 1 | Unweighted objective: 1 least squares, 0 legacy log objective |
| `al0` | `None` | Starting entropy weight; `None` selects the automatic scale |
| `alf` | 1.7 | Factor by which alpha decreases per step |
| `nal` | `None` | Alpha steps; `None` gives 0 (single alpha) unweighted, 25 weighted |
| `clim` | `None` | RMS target; `None` is 1.0 (relative RMS for `iw=0`) |
| `ns` | 7 | Default-map smoothing `exp(-r²/ns²)`, in map pixels |
| `ac`, `nim` | 8e-4, 1000 | Stationarity tolerance and iteration limit per alpha |
| `overs` / `map_size` | 0.4 / `None` | Map pixels per axis relative to channels, or an explicit even size |
| `pb0`, `pb1` | 0, 0 | Open phase interval excluded from the fit; equal values exclude nothing |
| `instrument_sigma` | 0 | Gaussian line-spread sigma (km/s, not FWHM) |
| `phase_subsamples` | 0 | Minimum exposure-quadrature samples (0 automatic) |
| `norm` | 0 | 1 rescales spectra to equal line flux (weighted fits also need `conditional_normalization=True`) |
| `baseline_fraction` | 0.1 | Baseline map level relative to the RMS flux |
| `af`, `wid` | 0, 10e5 | Empirical central-absorption prescription (amplitude, cm/s); 0 disables |
| `nbins` | 20 | Phase bins for diagnostic plots only |
| `run_timeout`, `max_memory_bytes` | 600 s, 2 GB | Resource limits per run |
| `threads` | `None` | Threads that build the projection; `None` uses up to four. The result does not depend on it |

## Conventions in brief

| Topic | Convention |
|---|---|
| Map axes | `(Vx, Vy)` in km/s, including `sin(i)` |
| Line-of-sight velocity | `v = gamma - Vx cos(2πφ) + Vy sin(2πφ)` |
| Overlays | Phase zero is the donor's inferior conjunction |
| Wavelength to velocity | Relativistic longitudinal Doppler formula |
| Map units | Input spectral-density units per km/s; the map's velocity-area integral equals the line integral over radial velocity |

PyDoppler does not apply barycentric or air/vacuum corrections. The data,
`lam0` and `gama` must already use consistent conventions.

The [scientific guide](https://github.com/fcotizelati/pydoppler/blob/master/docs/scientific_conventions.md)
gives the full definitions:

- the forward model;
- the likelihoods and covariance propagation;
- the optimizer and alpha schedule;
- units and validation.

## Outputs and reproducibility

All files are written to `dop.workdir`:

| File | Contents |
|---|---|
| `dopin` | Prepared spectra |
| `dopin.json` | Preprocessing settings and array hashes |
| `dop.in` | Solver parameters |
| `dop.out` | Map and data/model trails |
| `dop.status` | Convergence, RMS statistic, alpha, baseline |
| `dop.log` | Iteration log |
| `run.json` | Input, output and solver-source hashes, library versions, and success flag |

`Dopmap()` and `Reco()` refuse outputs that are missing from a successful
`run.json`. A lock file prevents two runs from sharing a directory.

Earlier versions compiled a Fortran solver in the run directory. Scripts
written for them should drop the `pydoppler.copy_fortran_code()` call, which no
longer exists. The Fortran sources, the `dopp` executable and `dop.precision`
left in old run directories are ignored and can be deleted.

To regenerate the README figures, run `python docs/make_readme_figures.py`.

## Scope and limitations

PyDoppler implements classical two-dimensional Doppler tomography of one line.
It assumes emission fixed in the co-rotating frame and a known ephemeris,
systemic velocity and rest wavelength.

It does not model:

- modulated or time-variable emission;
- occultations or out-of-plane velocities;
- several lines jointly;
- wavelength-dependent line-spread functions.

It does not provide posterior uncertainties for map pixels. Check features
against the following:

- the trailed spectra and residuals;
- the choice of continuum windows, map size and `clim`;
- phase coverage;
- injection and recovery tests on simulated data.

## Citation

If you use PyDoppler, please cite the original method and code, and this package:

* Spruit 1998, arXiv:astro-ph/9806141 ([ADS](https://ui.adsabs.harvard.edu/abs/1998astro.ph..6141S/abstract))
* Hernandez Santisteban 2021, ASCL:2106.003 ([ADS](https://ui.adsabs.harvard.edu/abs/2021ascl.soft06003H/abstract))

The original code and IDL interface are at <https://wwwmpa.mpa-garching.mpg.de/~henk/>.
Please report problems through the
[issue tracker](https://github.com/fcotizelati/pydoppler/issues).
