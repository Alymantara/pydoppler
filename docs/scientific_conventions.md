# Scientific conventions and validation

PyDoppler implements classical, two-dimensional, floating-default
maximum-entropy Doppler tomography. The map is a regularized estimate of a
stationary emission pattern in projected velocity space. It is not a spatial
image, a posterior sample, or a significance map.

## Coordinates, wavelengths, and exposures

The public map axes are `(Vx, Vy)` in km/s. Positive `Vx` points from the
primary to the donor and positive `Vy` follows the donor's orbital motion.
They include the projection factor `sin(inclination)`. The standard overlays
assume phase zero at inferior conjunction of the donor, with

```text
v_rest(phase) = -Vx cos(2 pi phase) + Vy sin(2 pi phase).
```

The usual nonrelativistic observed radial velocity is `gamma + v_rest`.
Input wavelengths and `lam0` must use the same air/vacuum convention and
barycentric/heliocentric reference frame. PyDoppler does not perform those
calibrations. Wavelengths are in Angstrom. Its exact wavelength convention is

```text
D(v) = sqrt((1 + v/c) / (1 - v/c)),  c = 299792.458 km/s
lambda_rest = lambda_observed / D(gama)
v_rest = c ((lambda_rest/lam0)^2 - 1) / ((lambda_rest/lam0)^2 + 1).
```

Thus positive `gama` removes a redshift. Relativistic velocity composition
differs slightly from simply subtracting velocities; an ephemeris or velocity
derived with an optical `cz` convention should be converted consistently.
`delw` is the **half-width** of the rest-frame wavelength window.

Phases are in cycles and refer to **mid-exposure**. `delta_phase`, or the third
phase-file column with `use_list_dpha=True`, is the full exposure width, not
the half-width. Every exposure enters the inversion separately. `nbins`
changes diagnostic plots only. Display binning computes periodic interval
overlaps, including exposures longer than one display bin.

The forward model integrates each exposure with midpoint quadrature. Its
automatic sampling places about three samples per smaller of a map pixel and
spectral channel along the largest projected velocity excursion; in tests
this keeps the error below about 0.05% of a pixel's peak response.
`phase_subsamples` can impose a larger minimum. `instrument_sigma` is the
Gaussian line-spread function **standard deviation in km/s**, with zero
meaning no instrumental broadening; `FWHM = sqrt(8 ln 2) * sigma`. Each map
pixel is a uniform square in velocity space. Its exact line-of-sight
projection at phase `phi` is a trapezoid (two boxes of half-widths
`dv |sin 2 pi phi| / 2` and `dv |cos 2 pi phi| / 2` convolved), integrated over
spectral-channel edges and convolved with the Gaussian LSF. Adjacent pixels
therefore tile exactly at every phase. (The original code used a triangular
pixel response of one pixel width, which does not tile and produced up to
50% ripple in projected smooth maps near phases 0, 0.25, 0.5 and 0.75.) The
LSF is constant across wavelength and exposure.

`pb0,pb1` define a periodically repeated, open phase-exclusion interval, with
`0 <= pb1-pb0 < 1`. Equal endpoints exclude nothing; `.95,1.05` excludes the
region around phase zero. Boundary contacts remain included to floating-point
precision. The same exclusion mask controls the objective, convergence
statistic, normalization factors, and residual diagnostics. Exclusion is by
the exposure's midpoint; exposures crossing an eclipse boundary require a
deliberate selection by the user.

## Preprocessing and uncertainty propagation

Spectra are samples of **spectral density**, not integrated pixel counts.
Convert integrated measurements and their uncertainties to densities first.
`Foldspec()` retains native spectra; `Dopin()` composes resampling and continuum
subtraction in one native-to-output operation. Linear interpolation preserves
sampled values; it is not a flux-conserving integration of pixel counts.
There is no clipping of negative continuum-subtracted measurements.

Continuum bands are specified in the **observed wavelength frame**. The
polynomial is subtracted, not divided out. Automatic bands select the outer
5% of the available coverage; they do not identify uncontaminated continuum.
Line wings in these bands will bias the inferred line and map. Choose clean
bands and test the sensitivity to their placement and polynomial order.

For native flux vector `f`, resampling operator `R`, continuum design matrix
`X`, and linear coefficient estimator `B`, the prepared data and covariance are

```text
y = (R - X B) f = A f
C_y = A C_f A^T.
```

For independent native pixels `C_f = diag(sigma_i^2)`. This includes both
interpolation covariance and continuum-fit covariance, including their cross
terms. In particular, interpolation with weights `a,b` has variance
`a^2 sigma_1^2 + b^2 sigma_2^2`, not a linear interpolation of variances.
Compact operators keep native covariance storage small. Dense matrices are
materialized only for the selected inversion window and eigendecomposition.

An explicit `Dopin(flux_err=...)` describes independent errors on the current
input grid; supply `flux_covariance=...` when those inputs are correlated.
Files loaded by `Foldspec()` can contain a third, strictly positive 1-sigma
error column. Nonfinite samples, nonpositive supplied errors, extrapolation,
and unidentifiable continuum fits are rejected.

Without supplied errors, weighted preparation estimates a constant noise
level per spectrum from continuum residuals, accounting for fitted polynomial
degrees of freedom. This assumes independent, homoscedastic native noise and
an adequate continuum model. These noise estimates are then treated as fixed;
their own uncertainty is not propagated. A perfectly noiseless or
underdetermined continuum cannot supply an error estimate and raises an error.
`estimate_errors=False` requires supplied valid errors.

## Likelihood, regularization, and stopping

For a statistical fit, set `iw=1` **before** calling `Dopin()`. With the default
`covariance_mode="full"`, `Syncdop()` uses the prepared full covariance (backend
mode `iw=2`) and a Gaussian objective, regardless of the `ih` setting.
For residual `r_p = y_p - model_p` in unmasked spectrum `p`,

```text
Q_data = sum_p r_p^T C_p^+ r_p
H = -Q_data / 2
fit_statistic = sqrt(Q_data / sum_p rank(C_p)).
```

`C^+` is computed from retained eigenmodes: eigenvalues greater than
`covariance_rcond * largest_eigenvalue`, with default `covariance_rcond=1e-10`.
This defines a fit in the resolved covariance subspace; discarded modes impose
no constraint on the model. It is important to check threshold sensitivity
when interpolation or continuum removal produces near-null modes. Different
exposures are assumed independent. Cross-exposure calibration errors and
uncertainty in gamma, phase, rest wavelength, and the LSF are not included.

`covariance_mode="diagonal"` deliberately approximates prepared channels as
independent and uses `sum((r/sigma)^2)`. It overstates the information in
resampled data: on the bundled U Gem spectra it reaches a given RMS target
only with much weaker regularization and visibly noisier maps than the
full-covariance fit. `Residuals()` returns `r/sigma` when actual
fitted errors are available, otherwise its second return value is `None`.
For a full-covariance run, `Residuals(whiten=True)` returns
`U Lambda^(-1/2) U^T r`. Its squared norm is `Q_data`. For rank-deficient
covariance these channel-shaped values have projector covariance, so the
individual values are not independent unit Gaussians.

The denominator above is measurement rank, **not residual degrees of freedom**:
`fit_statistic` and `clim` are RMS quantities, not reduced chi-square.
Regularization makes an effective number of fitted parameters nontrivial.
An RMS near one can be useful under an adequate Gaussian model, but does not
establish feature significance or model adequacy.

The default `iw=0` has no measurement-noise interpretation. Its statistic is
the relative RMS residual `sqrt(sum(r^2)/sum(y^2))` over unmasked data, so a
zero map scores exactly 1 and the default `clim=1` only requires the map to
explain part of the data. With the default `ih=1` the objective is
`H = -sum(r^2) / (2 nvp)` in units of the RMS flux; legacy `ih=0` uses
`H = sum(y log m - m) / nvp` after a baseline shift. That objective is not a
Poisson likelihood for calibrated or continuum-subtracted flux. The division
by the number of spectral channels `nvp` restores the original code's
normalization, so `alpha` does not depend on spectral sampling and the
historical `al0=0.002` keeps its meaning. Use the weighted Gaussian mode for
quantitative statistical work.

The optimized objective is `H + alpha*S`. For the positive latent map `x`
and its smoothed floating default `chi = W x`,

```text
S = -sum_j [x_j log(x_j/chi_j) - x_j + chi_j].
```

Derivatives include the dependence of `chi` on `x`. The smoothing kernel has
weights proportional to `exp(-distance^2/ns^2)` (Gaussian sigma
`ns/sqrt(2)` map pixels), with finite truncation, zeros outside the active
circular pixels, and reflected conditions at square-array edges. `ns` is a
regularization scale, distinct from `instrument_sigma` and
the optional smoothing of a displayed map.

A projected flat-map baseline is added internally so signed measurements can
be handled with a positive latent map. It is subtracted from both the returned
map and reconstruction. Its level is `baseline_fraction` (default 0.1) of the
internal RMS flux scale at the projected peak. In least-squares and Gaussian
fits this level is fixed, so the returned map is bounded below by minus the
baseline map (about -2% of the peak for the bundled data) and positivity
remains an effective constraint. The legacy log objective `ih=0` needs
positive shifted data, so it raises the level until every unmasked sample is
positive; a single negative noise excursion can then lift the baseline far
above the line flux and the map loses positivity. `dop.status` records the
level used (`baseline`, `baseline_raised`) and `Syncdop()` warns when it was
raised. The baseline affects the entropy prior and is not a physical
continuum model. The legacy `af,wid`
absorption prescription is empirical; the default `af=0` disables it.

The optimizer uses positivity-preserving, preconditioned gradient steps with
Newton acceleration and an Armijo line search. Let `gH` and `gS` be the
likelihood and entropy gradients. Its numerical stationarity criterion is

```text
delta = ||x * (gH + alpha*gS)||_2 /
        max(||x*gH||_2 + ||x*alpha*gS||_2, 1e-12),
```

where multiplication by `x` is elementwise. All calculations use an internal
RMS flux scale so changing physical flux units does not change this test.
It requires `delta <= ac`, or an absolute numerator at most `1e-12`, plus
the requested RMS target. This tests numerical stationarity in the scaled
coordinates, not scientific identifiability or the significance of faint
pixels. It never relaxes `ac` or accepts an iteration limit as convergence.
A failed run raises and cannot be displayed as a successful
managed output. Each accepted line-search point's projection, likelihood and
entropy are reused by the next iteration.

The alpha schedule starts at `al0`. With `al0=None` the solver uses
`0.002 * nvp * I`, where `I` is the mean per-channel information (the
likelihood curvature `1/sigma^2` in RMS-flux units, scaled by `1/nvp` for
unweighted fits): exactly 0.002 for unweighted fits and the equivalent
regularization for Gaussian fits. Alpha is divided by `alf` for up to `nal`
steps (`None`: 0 extra steps unweighted, 25 weighted), each warm-started from
the previous map, until the RMS statistic meets `clim` (`None`: 1). A crossing
is refined by bisection in `log(alpha)`, starting each trial from the
smoother bracketing map, until the statistic lies within 0.5% below `clim`,
so the accepted map is close to the most regularized one meeting the target.
The schedule stops early when the statistic stops improving (less than 0.1%
per step) or when two consecutive shrinking improvements extrapolate
geometrically, using the slower rate, to a floor above `clim`; the error
message and `rms_floor_estimate` report that floor. This is a discrepancy-
principle schedule, not evidence optimization or posterior inference.

In `dop.status`, the legacy-named `optimizer_status` is a combined status:
0 means success, 1 means the optimizer converged but the target was not met
(`rms_plateau=1` when the statistic levelled off), 11 means iteration limit,
12 means line-search failure, 13 means nonfinite objective, and 19 means
nonfinite fit statistic. The overall `converged` flag is one only when both
optimization and target checks succeed. `target_bracketed=0` with a schedule
means the first alpha already met the target; `Syncdop()` then warns that a
larger `al0` could give a smoother map.

`norm=0` preserves phase-dependent line flux and is the default. Opt-in
`norm=1` rescales each included spectrum to their mean line integral. In a
weighted fit this additionally requires `conditional_normalization=True`:
errors are transformed conditional on the fitted factors, whose uncertainty
and induced cross-exposure covariance are omitted. Factors are recorded in
`dop.normalization`; output trails are restored to original input units.

## Map output and binary overlays

Spectral centers and map centers are separate grids. The raw parsed map uses
`output.vmap` for both axes, while trails use `output.vp`. Plotting uses true
bin edges, phases increasing upwards, and a shared data-derived scale for
input/reconstructed trails. This makes an amplitude error visible.
`Dopmap()` returns its processed, normalized display array, including requested
radial-mean removal, Gaussian smoothing, or sign reversal; samples outside
circular support are NaN. Ordinary maps are divided by their positive peak,
so zero flux remains zero and any negative values remain negative (earlier
versions rescaled minimum-to-maximum, which placed zero flux part-way up the
colour scale); default limits are 0 and the 99.5th percentile. Radial-mean
subtracted maps are scaled by their absolute maximum with symmetric limits.
The raw solver map is available through the output parser; use it for
physical map intensities.

Raw map-density units from the updated solver are **input spectral-density
units per km/s**; older `dop.out` files may use the former conversion and must
not be assumed to share these units. The updated map's
velocity-area integral corresponds to the line-profile integral over radial
velocity when the spectral window captures the complete model. It is not
automatically a line flux, luminosity, or `F_nu`. If the input is `F_lambda`,
conversion to flux per unit velocity requires the wavelength-to-velocity
Jacobian; converting physical integrals back to wavelength also requires it.
Do not compare absolute maps without accounting for units, sampling, and
preprocessing. Spectral trails retain the original input density units.

`stream(q, k1, porb, m1, inc)` uses `q=M2/M1` and the supplied **projected**
primary velocity amplitude `k1` to set its velocity scale. Passing `k1=None`
derives it from masses, period, and inclination. Inconsistent redundant inputs
produce a warning. Resonance radii describe Keplerian orbital frequencies
`j/(j-k)` times the binary frequency and their overlays use projected
velocities. Roche geometry handles equal masses and uses bracketed roots;
the ballistic stream uses adaptive integration. These remain idealized
point-mass, circular-binary overlays.

## Reproducibility, tests, and limits

`Syncdop()` **executes every requested inversion** in the calling process,
with the NumPy/SciPy solver in `pydoppler/solver.py`. The projection is sparse
and exactly sized, so memory scales with the number of nonzero projection
weights. It is built on up to four threads (`threads`); the result does not
depend on their number. A work-directory lock
prevents concurrent collisions. `dopin.json` records continuum settings, noise
recipes, and hashes of the actual input arrays. `run.json` includes matching
preprocessing metadata, parameters, input and output hashes, the hash of the
solver source, the Python, NumPy and SciPy versions, diagnostics, and
success/failure;
it explicitly marks preprocessing metadata unavailable when manual input edits
break that association. These records do not replace upstream calibration or
ephemeris provenance. `dop.status` records convergence. Managed
plots reject unsuccessful or modified output. A memory estimate is checked
before projection and dense covariance preparation; `max_memory_bytes` and `run_timeout`
bound individual runs. The time limit is checked between projection work
blocks and optimizer iterations. Fortran sources, executables and
`dop.precision` files left in working directories by earlier versions are
ignored.

The solver follows, step by step, the arithmetic of the Fortran code it
replaced. On the bundled data both take the same number of iterations and
their maps agree to rounding (1e-12 of the peak or better). Rounding
differences accumulate over long optimizations, so results from different
platforms or library versions can differ at the level allowed by the
stationarity tolerance `ac`.

The automated suite includes independent variance/covariance calculations,
Monte Carlo propagation through continuum subtraction, finite-difference
objective derivatives, projection-adjoint checks, exact square-pixel
projection at a cardinal phase, alpha-schedule bracketing and plateau
detection, excluded-phase corruption,
intentional nonconvergence, run provenance, analytic S-wave recovery
in all four map quadrants, and comparisons with stored results of the former
Fortran solver (`tests/data/fortran_reference.npz`). Injection tests cover both signs of gamma and
flux/error scaling by `1e-15`. Geometry checks cover phase wraps, exposure
overlaps, map centers, Roche geometry (including equal masses), and projected
resonances.

The bundled U Gem example uses the phase file's 600 s exposure widths
(0.0393 cycles; Echevarria et al. 2007). The unweighted default fit converges
at alpha 0.002 with a relative RMS of 0.236 and a smooth map whose minimum is
about -2% of the peak. With continuum-estimated errors and full covariance,
targets of 1.3 and 1.2 are reached at alpha 21 and 1.1 (all phases, `norm=0`,
144 x 144 map). The statistic levels off near 1.07, so a target of 1.0 is
correctly reported as unattainable; the continuum-residual noise estimate is
about 8% below the scatter of line-free channels, which accounts for most of
that excess. These results do not establish the adequacy of that noise model
or the stationary-emission assumption for U Gem. Map resolution, physical
phase selection, noise assumptions, regularization, and iteration budget
require data-specific checks. Increasing a target merely to suppress an error
does not validate a reconstruction.

These are regression and numerical validation checks, not a proof that every
dataset is identifiable. For scientific use, vary the continuum bands,
window, map resolution, regularization, phase coverage, and plausible system
parameters; perform data-specific injection/recovery checks. Classical maps
do not model orbital modulation, occultations, out-of-plane velocities,
time-variable emissivity, multi-line data jointly, or wavelength-dependent
LSFs. Full uncertainty inference and those extended tomography models remain
outside this pipeline's scope.
