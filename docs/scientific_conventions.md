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
automatic sampling resolves the largest projected velocity excursion to
approximately one eighth of the smaller of a map pixel and spectral channel;
`phase_subsamples` can impose a larger minimum. Verify quadrature convergence
for especially demanding data. `instrument_sigma` is the Gaussian line-spread
function **standard deviation in km/s**, with zero meaning no instrumental
broadening; `FWHM = sqrt(8 ln 2) * sigma`. The map pixel response is a
normalized triangular kernel integrated over spectral-channel edges and
convolved with this Gaussian. This is a discretization convention, not an
exact projection of a uniform square pixel. The current LSF is constant
across wavelength and exposure.

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
mode `iw=2`) and a Gaussian objective, regardless of the legacy `ih` setting.
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
independent and uses `sum((r/sigma)^2)`. It is cheaper but can overstate the
information in resampled data. `Residuals()` returns `r/sigma` when actual
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

The compatibility default `iw=0` has no measurement-noise interpretation.
Its stopping statistic is
`sqrt(sum(r^2)/sum(y^2))` over unmasked data. With `ih=1` it minimizes unweighted
squared residuals; with legacy `ih=0` it uses a positive log objective after a
baseline shift. That objective is not a Poisson likelihood for calibrated or
continuum-subtracted flux. Use the weighted Gaussian mode for quantitative
statistical work.

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
map and reconstruction. `baseline_fraction` (default 0.1 of internally scaled
flux) sets a minimum; the solver raises it when needed to accommodate negative
samples. This affects the entropy prior and is not a physical continuum model.
Returned map values can consequently be negative. The legacy `af,wid`
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
managed output. Increasing `nim` can help a slow fit; an infeasible `clim`,
poor coverage, an inappropriate model, or excessive regularization requires
scientific reconsideration. `al0` sets the initial entropy weight; positive
`nal` permits that many values, successively divided by `alf`. With `nal=0`
only the first value is tried. The first converged fit meeting `clim` is used;
this schedule is not evidence optimization or posterior inference.
In `dop.status`, the legacy-named `optimizer_status` is a combined status:
0 means success, 1 means the optimizer converged but the alpha schedule missed
the target, 11 means iteration limit, 12 means line-search failure, 13 means
nonfinite objective, and 19 means nonfinite fit statistic. The overall
`converged` flag is one only when both optimization and target checks succeed.

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
circular support are NaN. Ordinary maps are scaled to [0,1]; radial-mean
subtracted maps are scaled by their absolute maximum. The raw solver map is
available through the output parser; use it for physical map intensities.

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

`Syncdop()` hashes sources, compiler identity, and the executable to cache the
build, then **executes every requested inversion**. A work-directory lock
prevents concurrent collisions. `dopin.json` records continuum settings, noise
recipes, and hashes of the actual input arrays. `run.json` includes matching
preprocessing metadata, parameters, hashes, diagnostics, and success/failure;
it explicitly marks preprocessing metadata unavailable when manual input edits
break that association. These records do not replace upstream calibration or
ephemeris provenance. `dop.status` records convergence. Managed
plots reject unsuccessful or modified output. A memory estimate is checked
before dense covariance preparation; `max_memory_bytes` and `run_timeout`
bound individual runs. Refresh old working directories explicitly with
`copy_fortran_code(workdir, overwrite=True)`; the default installer preserves
existing files.

The automated suite includes independent variance/covariance calculations,
Monte Carlo propagation through continuum subtraction, finite-difference
objective derivatives, projection-adjoint checks, excluded-phase corruption,
intentional nonconvergence, run/cache provenance, and analytic S-wave recovery
in all four map quadrants. Injection tests cover both signs of gamma and
flux/error scaling by `1e-15`. Geometry checks cover phase wraps, exposure
overlaps, map centers, equal-mass Roche geometry, and projected resonances.

The bundled unweighted U Gem example converges and its figures have been
regenerated. A separate full-covariance experiment on those real spectra
exposes a remaining practical limitation: a 144-by-144 map with the legacy
weak entropy weight `al0=0.002` did not reach stationarity in 1000 iterations.
Stronger weights 0.2 and 1 also exhausted a 3000-iteration budget. Coarser
32-by-32 and 48-by-48 maps reached the unchanged stationarity tolerance, but
their RMS values were 1.281 and 1.147, respectively, so a requested target of
1.0 was correctly rejected. These tests used continuum-estimated errors,
`norm=0`, and all phases, including eclipse phases. They do not establish the
adequacy of that noise model or the stationary-emission assumption for U Gem.
Map resolution, physical phase selection, noise assumptions, regularization,
and iteration budget require data-specific checks. Increasing a target merely
to suppress an error does not validate a reconstruction.

These are regression and numerical validation checks, not a proof that every
dataset is identifiable. For scientific use, vary the continuum bands,
window, map resolution, regularization, phase coverage, and plausible system
parameters; perform data-specific injection/recovery checks. Classical maps
do not model orbital modulation, occultations, out-of-plane velocities,
time-variable emissivity, multi-line data jointly, or wavelength-dependent
LSFs. Full uncertainty inference and those extended tomography models remain
outside this pipeline's scope.
