# Changelog

## 1.0.2 (2026-11-13)

### Fixed
- `Reco` scaled the observed and reconstructed trails each by their own
  maximum, so a single cosmic ray in the data squashed the observed panel
  while the reconstruction saturated. Both panels now share one
  outlier-robust scale (99.5th percentile of the observed trail) and are
  directly comparable. Plot only; maps are unchanged.

Found while testing on public ESO X-shooter spectra of GP Com (programme
084.D-0814A; the same data as Kupfer et al. 2016, arXiv:1601.02841).

## 1.0.1 (2026-11-13)

### Fixed
- `Foldspec` crashed with `IndexError` when `nbins` was larger than the number
  of spectra (e.g. after removing eclipse spectra). Its unused phase-binning
  loop is removed; `nbins` only sets the phase grid of the trail plot, and the
  Doppler map is unchanged.
- `Dopin`: the average spectrum was divided by `nbins` instead of the number of
  spectra, mis-scaling the plot when they differ, and the lower panel's y-limit
  subtracted a spurious 1.0, flipping the axis for continuum-normalised spectra.
  Plot only; the per-spectrum normalisation was already correct.

Found while testing on the public IP Peg outburst spectra of Papadaki, Boffin
& Steeghs (2008, The Journal of Astronomical Data 14, 2).

## 1.0.0 (2026-11-12)

First stable release: a pure Python implementation of Spruit's Doppler
tomography code.

### Changed
- The Fortran `dop.f` engine and `gfortran`/`make` build step are replaced by
  `pydoppler/mem_engine.py` (NumPy/SciPy). No compiler is required.
- The projection matrix is built in vectorised form as a sparse matrix. A full
  run on the U Gem test data takes about 5 s, and the tomography step about 0.7 s.
- `Dopmap` defaults to a linear stretch and the `magma_r` colour map, and
  masks outside the |v| < va circle. `stream()` lines default to black.
- `Reco` uses constrained layout (the colour bar no longer overlaps the
  panels), and its default colour limits come from percentiles.
- `Syncdop` reports the sparse projection-matrix size.

### Fixed (relative to the pre-release native port)
- Projection geometry: removed a +1 velocity-channel offset (+6.7 km/s on
  the test data) from mixing 1-based Fortran bin edges with 0-based data indices.
- MEM: entropy value and gradient (including the floating-default term),
  second derivatives in the Lucy acceleration step, default-map kernel width
  and the baseline/normalisation order now follow Spruit's `dop.f`. The
  reconstruction now converges instead of stalling at the first alpha.
- Alpha schedule: no longer diverges when chi stops decreasing at the noise
  floor of the data.
- `Reco` drew the phase axis upside down.
- Plot labels render with `text.usetex = True`.

### Notes
- Developed with the help of Claude (Anthropic), using Claude Opus 5.5
  (model `claude-opus-5-5`) in Claude Code.
