# Changelog

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
