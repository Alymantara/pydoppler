"""Floating-default maximum-entropy Doppler tomography: numerical core.

Based on H. C. Spruit's 1999 code (dopmap, version 2.0), with a revised
forward model, likelihoods, optimizer and regularization schedule. This
NumPy/SciPy module replaces the Fortran solver of earlier versions and
follows the same arithmetic step by step.

Units follow the ``dopin`` file: phases and exposure widths are in cycles;
velocities, ``instrument_sigma`` and ``wid`` are in cm/s. Data arrays have
shape ``(nph, nvp)`` with flat index ``ip*nvp + iv``. Map pixels inside the
circular support are stored x-major, and ``image[ix, iy]`` is the full grid.
"""

from __future__ import annotations

import os
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import sparse
from scipy.ndimage import correlate1d
from scipy.special import erfc

PI2 = 8.0 * np.arctan(1.0)
TINY = np.finfo(float).tiny
MAX_HALF_WIDTH = 100
# Exposure-quadrature samples per (smaller of map pixel, channel) of
# radial-velocity travel; 3 keeps errors below ~0.05% of a pixel's peak response.
QUADRATURE_DENSITY = 3.0
MAX_PROJECTION_THREADS = 4
# Array elements per projection work block; bounds the temporary memory.
_BLOCK_ELEMENTS = 400_000


class SolverError(RuntimeError):
    """The data or settings cannot be inverted."""


class SolverTimeout(SolverError):
    """The run exceeded its time limit."""


def phase_mask(phases, start, stop):
    """Valid phases outside the open, periodically repeated exclusion interval."""
    if not np.isfinite([start, stop]).all() or not 0 <= stop - start < 1:
        raise ValueError("Phase exclusion must satisfy 0 <= pb1 - pb0 < 1 cycle.")
    distance = np.mod(np.asarray(phases, dtype=float) - start, 1.0)
    edge = 16 * np.finfo(float).eps
    return ~((distance > edge) & (distance < stop - start - edge))


def _normal_first(x, sigma):
    """First antiderivative of the Gaussian CDF (of the step if sigma=0)."""
    if sigma <= 0:
        return np.maximum(x, 0.0)
    z = x / sigma
    return x * 0.5 * erfc(-z / np.sqrt(2.0)) + sigma * np.exp(-0.5 * z * z) / np.sqrt(PI2)


def _normal_second(x, sigma):
    """Second antiderivative of the Gaussian CDF (of the step if sigma=0)."""
    if sigma <= 0:
        return 0.5 * np.maximum(x, 0.0) ** 2
    z = x / sigma
    cdf = 0.5 * erfc(-z / np.sqrt(2.0))
    pdf = np.exp(-0.5 * z * z) / np.sqrt(PI2)
    return 0.5 * ((x * x + sigma * sigma) * cdf + x * sigma * pdf)


def square_cdf(x, ha, hb, sigma=0.0):
    """CDF of a uniform square map pixel projected onto the line of sight.

    The projection is the convolution of boxes of half-widths ``ha`` and
    ``hb`` (a trapezoid), optionally convolved with a Gaussian line-spread
    function. Adjacent pixels tile exactly at every orbital phase.
    """
    h1 = np.maximum(ha, hb)
    h2 = np.minimum(ha, hb)
    narrow = h2 <= 1.0e-6 * h1
    h2s = np.where(narrow, 1.0, h2)
    a = x + h1
    b = x - h1
    z = (_normal_second(a + h2s, sigma) - _normal_second(a - h2s, sigma)
         - _normal_second(b + h2s, sigma) + _normal_second(b - h2s, sigma)) / (4 * h1 * h2s)
    if np.any(narrow):
        box = (_normal_first(a, sigma) - _normal_first(b, sigma)) / (2 * h1)
        z = np.where(narrow, box, z)
    return np.clip(z, 0.0, 1.0)


def _trapezoid_cdf_inplace(x, ha, hb):
    """``square_cdf(x, ha[:, None], hb[:, None], 0)`` computed in the storage of *x*."""
    h1 = np.maximum(ha, hb)
    h2 = np.minimum(ha, hb)
    narrow = h2 <= 1.0e-6 * h1
    any_narrow = bool(narrow.any())
    if any_narrow:
        saved = x[narrow]
        h2 = np.where(narrow, 1.0, h2)
    h1c, h2c = h1[:, None], h2[:, None]
    a = x + h1c
    np.subtract(x, h1c, out=x)                      # x now holds x - h1
    acc = a + h2c                                   # max(x + h1 + h2, 0)^2
    np.maximum(acc, 0.0, out=acc)
    np.multiply(acc, acc, out=acc)
    np.subtract(a, h2c, out=a)                      # - max(x + h1 - h2, 0)^2
    np.maximum(a, 0.0, out=a)
    np.multiply(a, a, out=a)
    np.subtract(acc, a, out=acc)
    np.add(x, h2c, out=a)                           # - max(x - h1 + h2, 0)^2
    np.maximum(a, 0.0, out=a)
    np.multiply(a, a, out=a)
    np.subtract(acc, a, out=acc)
    np.subtract(x, h2c, out=x)                      # + max(x - h1 - h2, 0)^2
    np.maximum(x, 0.0, out=x)
    np.multiply(x, x, out=x)
    np.add(acc, x, out=acc)
    # The common factor 1/2 of the four terms is exact, so it joins the divisor.
    np.divide(acc, ((8.0 * h1) * h2)[:, None], out=acc)
    if any_narrow:
        hn = h1[narrow][:, None]
        acc[narrow] = (np.maximum(saved + hn, 0.0) - np.maximum(saved - hn, 0.0)) / (2 * hn)
    np.clip(acc, 0.0, 1.0, out=acc)
    return acc


def _capacity_bound(travel, pixel_reach, narrowest, nvp):
    """Upper bound on the number of projection weights.

    An exposure reaches at most the channels between its extreme line-of-sight
    velocities (``travel`` apart at most), widened on both sides by the pixel
    and the line-spread function (``pixel_reach``).
    """
    channels = np.ceil((travel + 2 * pixel_reach) / narrowest) + 3
    return int(np.minimum(nvp, channels).sum())


def build_projection(pha, dpha, vp, nv, *, instrument_sigma=0.0, phase_subsamples=0,
                     wid=1.0, af=0.0, threads=None, deadline=None):
    """Exact square-pixel projection with exposure-phase quadrature.

    ``pha`` are mid-exposure phases and ``dpha`` full exposure widths, both in
    radians; ``vp`` are increasing spectral bin centres. Each exposure is
    integrated by midpoint quadrature with about three samples per bin of
    radial-velocity travel.

    Returns ``(matrix, indx, indy, va, dvmap)``: a sparse ``(nph*nvp, nol)``
    matrix mapping the pixels inside the circular support to the data, the
    zero-based grid indices of those pixels, the map half-width and the pixel
    size.
    """
    pha = np.asarray(pha, dtype=float)
    dpha = np.asarray(dpha, dtype=float)
    vp = np.asarray(vp, dtype=float)
    nph, nvp = pha.size, vp.size
    nil = nph * nvp
    if nv < 4:
        raise SolverError("Map size must be at least four pixels.")
    va = float(max(abs(vp[0]), abs(vp[-1])))
    dvmap = 2 * va / (nv - 1)
    dvi = float(np.min(np.diff(vp)))
    if dvi <= 0:
        raise SolverError("Velocities must increase.")
    edges = np.empty(nvp + 1)
    edges[1:-1] = (vp[:-1] + vp[1:]) / 2
    edges[0] = vp[0] - (vp[1] - vp[0]) / 2
    edges[-1] = vp[-1] + (vp[-1] - vp[-2]) / 2
    narrowest = np.min(np.diff(edges))
    centre = (nv + 1) / 2.0
    ix, iy = np.meshgrid(np.arange(1, nv + 1), np.arange(1, nv + 1), indexing="ij")
    vxg = dvmap * (ix - centre)
    vyg = dvmap * (iy - centre)
    inside = ~(vxg ** 2 + vyg ** 2 >= (va + dvmap / 10) ** 2)
    indx, indy = ix[inside] - 1, iy[inside] - 1
    vx, vy = vxg[inside], vyg[inside]
    nol = vx.size
    radius = np.sqrt(vx * vx + vy * vy)
    nq = np.ceil(np.abs(dpha)[None, :] * radius[:, None] / min(dvmap, dvi) * QUADRATURE_DENSITY)
    nq = np.maximum(np.maximum(nq, 1).astype(np.int64), int(phase_subsamples))
    if nq.max() > 100000:
        raise SolverError("Exposure integration requires too many samples.")

    # Work blocks of whole pixels, bounded in quadrature samples and in the
    # dense (pixel, phase, channel) accumulator.
    half = 0.5 * dvmap
    slots_guess = int(2 * (1.5 * half + 8 * instrument_sigma) / narrowest) + 3
    max_samples = max(1000, _BLOCK_ELEMENTS // (slots_guess + 1))
    max_pixels = max(1, 10 * _BLOCK_ELEMENTS // nil)
    samples = np.cumsum(nq.sum(axis=1))
    bounds = [0]
    while bounds[-1] < nol:
        a = bounds[-1]
        done = samples[a - 1] if a else 0
        b = int(np.searchsorted(samples, done + max_samples, side="right"))
        bounds.append(min(nol, max(a + 1, min(b, a + max_pixels))))

    def block(limits):
        a, b = limits
        if deadline is not None and time.monotonic() > deadline:
            raise SolverTimeout("The projection exceeded the run time limit.")
        npix = b - a
        reps = nq[a:b].ravel()                      # one entry per (pixel, phase)
        npair = reps.size
        first = np.cumsum(reps) - reps
        iq = np.arange(int(first[-1] + reps[-1])) - np.repeat(first, reps)
        nqs = np.repeat(reps, reps)
        ph = (np.repeat(np.tile(pha, npix), reps)
              + np.repeat(np.tile(dpha, npix), reps) * ((iq + 0.5) / nqs - 0.5))
        sp, cp = np.sin(ph), np.cos(ph)
        vr = (np.repeat(np.repeat(vx[a:b], nph), reps) * sp
              - np.repeat(np.repeat(vy[a:b], nph), reps) * cp)
        np.abs(sp, out=sp)
        np.abs(cp, out=cp)
        h1 = half * sp
        h2 = half * cp
        span = h1 + h2 + 8 * instrument_sigma
        upper = vr + span
        # First channel reached by each sample: its lower edge is the last
        # one not above vr - span (the first channel if the sample starts below).
        c0 = np.searchsorted(edges, vr - span, side="right")
        np.maximum(c0, 1, out=c0)
        c0 -= 1
        slots = int(2 * float(span.max()) / narrowest * (1 + 1e-9)) + 2
        while True:
            padded = np.concatenate((edges, np.full(slots + 1, edges[-1])))
            window = np.take(sliding_window_view(padded, slots + 1), c0, axis=0)
            if not np.any((window[:, slots] <= upper) & (c0 + slots < nvp)):
                break
            slots += 2
        x = window - vr[:, None]
        if instrument_sigma > 0:
            cdf = square_cdf(x, h1[:, None], h2[:, None], instrument_sigma)
        else:
            cdf = _trapezoid_cdf_inplace(x, h1, h2)
        weight = cdf[:, 1:] - cdf[:, :-1]
        np.maximum(weight, 0.0, out=weight)
        denominator = (nqs * nph).astype(float)
        if af != 0:
            vfa = 1 + af * np.exp(-np.minimum((vr / wid) ** 2, 50.0))
            if np.any(vfa <= 0):
                raise SolverError("Absorption prescription makes the response nonpositive.")
            denominator = denominator * vfa
        np.divide(weight, denominator[:, None], out=weight)
        # Channels whose lower edge lies above vr + span are not reached.
        # Slots past the last channel repeat the final edge and weigh zero.
        np.putmask(weight, window[:, :slots] > upper[:, None], 0.0)
        target = (np.repeat(np.arange(npair) * nvp, reps) + c0)[:, None] + np.arange(slots)
        total = np.bincount(target.ravel(), weights=weight.ravel(),
                            minlength=npair * nvp)[:npair * nvp]
        nonzero = np.flatnonzero(total > 0)
        return nonzero % nil, total[nonzero], np.bincount(nonzero // nil, minlength=npix)

    # Storage for the weights, sized by an upper bound so that finished blocks
    # are copied once. Pages that stay unused are never touched.
    capacity = _capacity_bound(np.abs(dpha)[None, :] * radius[:, None],
                               1.5 * half + 8 * instrument_sigma, narrowest, nvp)
    index_type = np.int32 if max(capacity, nil) < 2 ** 31 - 1 else np.int64
    indices = np.empty(capacity, dtype=index_type)
    values = np.empty(capacity)
    counts = np.zeros(nol, dtype=np.int64)
    nnz = 0

    def store(results):
        nonlocal indices, values, nnz
        for (a, b), (index, value, count) in zip(spans, results):
            stop = nnz + index.size
            if stop > indices.size:                 # not reached while the bound holds
                indices = np.concatenate((indices[:nnz], np.empty(stop, dtype=index_type)))
                values = np.concatenate((values[:nnz], np.empty(stop)))
            indices[nnz:stop] = index
            values[nnz:stop] = value
            counts[a:b] = count
            nnz = stop

    spans = list(zip(bounds[:-1], bounds[1:]))
    if threads is None:
        threads = min(MAX_PROJECTION_THREADS, os.cpu_count() or 1)
    if threads > 1 and len(spans) > 1:
        pool = ThreadPoolExecutor(int(threads))
        try:
            store(pool.map(block, spans))
        finally:
            pool.shutdown(wait=True, cancel_futures=True)   # drop queued blocks on failure
    else:
        store(map(block, spans))
    if nnz >= 2 ** 31 - 1:
        index_type = np.int64
    indptr = np.zeros(nol + 1, dtype=index_type)
    np.cumsum(counts, out=indptr[1:])
    matrix = sparse.csc_matrix((values[:nnz], indices[:nnz], indptr), shape=(nil, nol))
    return matrix, indx, indy, va, dvmap


@dataclass
class SolverResult:
    """Outcome of one inversion.

    ``status`` holds the entries of ``dop.status``. After a failed run the
    arrays are ``None``. ``image[ix, iy]`` is the map density in input
    spectral-density units per km/s; ``data`` and ``model`` have shape
    ``(nph, nvp)`` in the original input units.
    """

    status: dict
    image: Optional[np.ndarray] = None
    data: Optional[np.ndarray] = None
    model: Optional[np.ndarray] = None
    change_maps: Optional[np.ndarray] = None
    normalization: Optional[np.ndarray] = None
    va: float = 0.0
    identification: float = 0.0

    @property
    def converged(self) -> bool:
        return self.status["converged"] == 1


class Tomography:
    """One inversion problem: data, statistical model, forward model, optimizer."""

    def __init__(self, phases, widths, velocity, flux, sigma=None, precision=None, ranks=None,
                 *, ih=1, iw=0, pb0=0.0, pb1=0.0, ns=7, ac=8e-4, nim=1000, al0=0.0, alf=1.7,
                 nal=0, clim=1.0, norm=0, wid=1.0e6, af=0.0, instrument_sigma=0.0,
                 phase_subsamples=0, baseline_fraction=0.1, map_size=4, threads=None,
                 deadline=None, log: Optional[Callable[[str], None]] = None):
        flux = np.asarray(flux, dtype=float)
        if flux.ndim != 2:
            raise ValueError("flux must have shape (nph, nvp).")
        self.nph, self.nvp = flux.shape
        self.nil = self.nph * self.nvp
        self.pha = np.asarray(phases, dtype=float) * PI2
        self.dpha = np.asarray(widths, dtype=float) * PI2
        self.vp = np.asarray(velocity, dtype=float)
        if self.pha.shape != (self.nph,) or self.dpha.shape != (self.nph,):
            raise ValueError("phases and widths must have one entry per spectrum.")
        if self.vp.shape != (self.nvp,) or self.nvp < 3:
            raise ValueError("velocity must have one entry per channel, at least three.")
        if ih not in (0, 1) or iw not in (0, 1, 2):
            raise ValueError("Invalid ih or iw.")
        if ns < 1 or ac <= 0 or nim < 1 or alf <= 1 or clim <= 0:
            raise ValueError("Invalid ns, ac, nim, alf or clim.")
        if norm not in (0, 1):
            raise ValueError("norm must be 0 or 1.")
        if instrument_sigma < 0 or phase_subsamples < 0 or baseline_fraction < 0:
            raise ValueError("Invalid LSF sigma, phase subsamples or baseline fraction.")
        if map_size < 4 or map_size % 2:
            raise ValueError("Map size must be even and >= 4.")
        if wid <= 0:
            wid, af = 1.0, 0.0
        self.nv = int(map_size)
        self.iw = int(iw)
        self.smoothing_width = int(ns)
        self.ac, self.nim = float(ac), int(nim)
        self.al0, self.alf, self.nal, self.clim = float(al0), float(alf), int(nal), float(clim)
        self.wid, self.af = float(wid), float(af)
        self.instrument_sigma = float(instrument_sigma)
        self.phase_subsamples = int(phase_subsamples)
        self.baseline_fraction = float(baseline_fraction)
        self.threads, self.deadline = threads, deadline
        self._log = log if log is not None else (lambda text: None)
        self.phase_valid = phase_mask(np.asarray(phases, dtype=float), pb0, pb1)
        self._statistical_input(flux, sigma, precision, ranks, int(norm))
        # Supplied uncertainties imply Gaussian data.
        self.ih = 1 if self.iw >= 1 else int(ih)
        # Unweighted objectives are scaled per spectral channel, so alpha does
        # not depend on spectral sampling and keeps its historical meaning.
        self.likelihood_scale = 1.0 / self.nvp if self.iw == 0 else 1.0
        self.matrix = None
        self.objective = 0.0
        self.stationarity = 0.0
        self.baseline_level = 0.0
        self.baseline_raised = False
        self.rms_floor_estimate = -1.0
        self._kernel = None
        self._check_deadline()

    def _check_deadline(self):
        if self.deadline is not None and time.monotonic() > self.deadline:
            raise SolverTimeout("The inversion exceeded the run time limit.")

    # -- statistical model -------------------------------------------------
    def _statistical_input(self, flux, sigma, precision, ranks, norm):
        nph, nvp, valid = self.nph, self.nvp, self.phase_valid
        if not valid.any():
            raise SolverError("No unmasked orbital phases remain.")
        dat = np.array(flux, dtype=float)
        if not np.isfinite(dat).all():
            raise SolverError("Non-finite input data.")
        if self.iw >= 1:
            if sigma is None:
                raise ValueError("Weighted fits require the standard deviations.")
            eb = np.array(sigma, dtype=float)
            if eb.shape != dat.shape or not np.isfinite(eb).all() or np.any(eb <= 0):
                raise SolverError("Invalid error bars.")
        else:
            eb = np.ones_like(dat)
        normalization = np.ones(nph)
        rank = np.full(nph, nvp, dtype=int)
        self.precision = None
        if self.iw == 2:
            if precision is None or ranks is None:
                raise ValueError("Full-covariance fits require precision matrices and ranks.")
            self.precision = np.array(precision, dtype=float)
            rank = np.asarray(ranks, dtype=int)
            if self.precision.shape != (nph, nvp, nvp) or rank.shape != (nph,):
                raise SolverError("The precision matrices have the wrong shape.")
            if np.any(rank < 1) or np.any(rank > nvp):
                raise SolverError("Invalid covariance rank.")
            if not np.isfinite(self.precision).all():
                raise SolverError("Non-finite precision.")
        # Optional normalization to the mean continuum-subtracted line integral.
        if norm == 1:
            line = dat[valid].sum(axis=1)
            if np.any(line <= 0):
                raise SolverError("Phase normalization requires positive line flux.")
            normalization[valid] = (line.sum() / valid.sum()) / line
        dat *= normalization[:, None]
        eb *= normalization[:, None]
        # An internal RMS flux scale makes all numerical choices unit invariant.
        vmax = np.abs(dat[valid]).max()
        if vmax <= 0:
            raise SolverError("No nonzero unmasked line flux remains.")
        self.flux_scale = float(vmax * np.sqrt(np.sum((dat[valid] / vmax) ** 2)
                                               / (valid.sum() * nvp)))
        dat /= self.flux_scale
        eb /= self.flux_scale
        self.statistic_rank = int(rank[valid].sum())
        if self.precision is not None:
            self.precision *= ((self.flux_scale / normalization) ** 2)[:, None, None]
        self.dat, self.eb, self.normalization = dat, eb, normalization

    def set_weights(self):
        """Diagonal weights, and the mean per-channel information.

        The information is the curvature of the likelihood and sets the
        automatic alpha scale. It uses the marginal precision 1/sigma^2 also
        for full covariance: the diagonal of an inverted near-singular
        covariance would overstate it.
        """
        valid = self.phase_valid
        weights = np.zeros((self.nph, self.nvp))
        weights[valid] = 1.0 / self.eb[valid] ** 2 if self.iw == 1 else 1.0
        used = int(valid.sum()) * self.nvp
        total = float(np.sum(1.0 / self.eb[valid] ** 2)) if self.iw >= 1 else float(used)
        self.stat_weights = weights
        self.information = self.likelihood_scale * total / used

    def apply_precision(self, residual):
        """Return ``(W r, r^T W r)``; masked spectra contribute nothing."""
        r = residual.reshape(self.nph, self.nvp)
        if self.iw == 2:
            z = np.matmul(self.precision, r[:, :, None])[:, :, 0]
            z[~self.phase_valid] = 0.0
        else:
            z = self.stat_weights * r
        return z.ravel(), float(np.sum(r * z))

    # -- forward model -----------------------------------------------------
    def build_projection(self):
        self.matrix, self.indx, self.indy, self.va, self.dvmap = build_projection(
            self.pha, self.dpha, self.vp, self.nv, instrument_sigma=self.instrument_sigma,
            phase_subsamples=self.phase_subsamples, wid=self.wid, af=self.af,
            threads=self.threads, deadline=self.deadline)
        self._adjoint = self.matrix.T               # CSR view of the same storage
        self.nol = self.matrix.shape[1]
        self.change_maps = np.zeros((4, self.nol), dtype=np.float32)

    def fproj(self, x):
        """Forward projection: map (nol) -> data (nil)."""
        return self.matrix @ x

    def bproj(self, d):
        """Adjoint (back) projection: data (nil) -> map (nol)."""
        return self._adjoint @ d

    def baseline(self, dat):
        """Projected flat map that lets a positive latent map describe signed data.

        Its level is ``baseline_fraction`` of the RMS flux at the projected
        peak. Gaussian objectives (ih=1) use this fixed level, so the returned
        map is bounded below by minus the baseline map. The legacy log
        objective (ih=0) needs positive shifted data and raises the level when
        necessary.
        """
        datb = self.fproj(np.ones(self.nol))
        d0 = datb.max()
        if d0 <= 0:
            raise SolverError("The flat map has no projected flux.")
        level = max(self.baseline_fraction, 1.0e-10)
        self.baseline_raised = False
        if self.ih == 0:
            usable = np.repeat(self.phase_valid, self.nvp) & (datb > 0)
            if usable.any():
                needed = float(np.max((-dat[usable] + 1.0e-10) * d0 / datb[usable]))
                if needed > level:
                    level, self.baseline_raised = needed, True
        self.baseline_level = level
        return level * datb / d0, np.full(self.nol, level / d0)

    # -- entropy -----------------------------------------------------------
    def chims(self, ps):
        """Default map: Gaussian-smoothed map, zero outside the circular support.

        The kernel is separable, ``exp(-r^2/ns^2)`` truncated at 1.5 ns, with
        half-sample reflection at the edges of the square array.
        """
        nv = self.nv
        if self._kernel is None:
            half = min(int(1.5 * self.smoothing_width), MAX_HALF_WIDTH, nv - 1)
            offsets = np.arange(-half, half + 1, dtype=float)
            kernel = np.exp(-offsets * offsets / float(self.smoothing_width) ** 2)
            self._kernel = kernel / kernel.sum()
        grid = np.zeros((nv, nv))
        grid[self.indx, self.indy] = ps
        grid = correlate1d(grid, self._kernel, axis=1, mode="reflect")
        grid = correlate1d(grid, self._kernel, axis=0, mode="reflect")
        return grid[self.indx, self.indy]

    def entr(self, ps):
        """Entropy S = -sum[x log(x/chi) - x + chi] with chi = W x.

        Returns ``(S, dS/dx, chi, log(x/chi))``. The gradient includes the
        dependence of chi on x (W is symmetric).
        """
        chi = self.chims(ps)
        tem = np.log(ps) - np.log(chi)
        s = float(np.sum(-ps * (tem - 1) - chi))
        sp = -tem + self.chims((ps - chi) / chi)
        return s, sp, chi, tem

    # -- likelihood --------------------------------------------------------
    def like(self, pht, ph):
        """Likelihood H and its gradient with respect to the map."""
        if self.ih == 0 and self.iw == 0:
            used = self.stat_weights.ravel() != 0
            pp = np.maximum(ph, TINY)
            h = float(np.sum((pht * np.log(pp) - pp)[used]))
            phh = np.where(used, pht / pp - 1, 0.0)
        else:
            phh, quadratic = self.apply_precision(pht - ph)
            h = -0.5 * quadratic
        return h * self.likelihood_scale, self.bproj(phh * self.likelihood_scale)

    def lamun(self, pht, ps, d1, d2, ph, chi, tem, al):
        """Newton step in the plane of two directions.

        Uses exact first and second directional derivatives of Q = H + alpha S
        and returns the coefficients ``(fla, fmu)`` of ``d1`` and ``d2``.
        """
        dph1, dph2 = self.fproj(d1), self.fproj(d2)
        if self.ih == 0 and self.iw == 0:
            used = self.stat_weights.ravel() != 0
            target = pht[used]
            pp = np.maximum(ph[used], TINY)
            dpp = target / pp - 1
            xl, xm = dph1[used] / pp, dph2[used] / pp
            hl, hm = np.sum(dph1[used] * dpp), np.sum(dph2[used] * dpp)
            hll, hlm, hmm = (-np.sum(xl ** 2 * target), -np.sum(xl * xm * target),
                             -np.sum(xm ** 2 * target))
        else:
            rw, _ = self.apply_precision(pht - ph)
            _, q1 = self.apply_precision(dph1)
            z2, q2 = self.apply_precision(dph2)
            hl, hm, hlm, hll, hmm = dph1 @ rw, dph2 @ rw, -(dph1 @ z2), -q1, -q2
        hl, hll, hm, hlm, hmm = (value * self.likelihood_scale
                                 for value in (hl, hll, hm, hlm, hmm))
        sl = sll = sm = slm = smm = 0.0
        if al > 0:
            dch1, dch2 = self.chims(d1), self.chims(d2)
            pl, pm = d1 / ps, d2 / ps
            tl, tm = dch1 / chi, dch2 / chi
            pz = ps - chi
            sl, sm = np.sum(-d1 * tem + pz * tl), np.sum(-d2 * tem + pz * tm)
            sll = -np.sum(ps * (pl - tl) ** 2)
            slm = -np.sum(ps * (pl - tl) * (pm - tm))
            smm = -np.sum(ps * (pm - tm) ** 2)
        ql, qll, qlm = hl + al * sl, hll + al * sll, hlm + al * slm
        qm, qmm = hm + al * sm, hmm + al * smm
        dd = qmm * qll - qlm ** 2
        denominator = qll ** 2 + qmm ** 2
        fla, fmu = 0.0, 1.0
        if denominator <= TINY:
            return fla, fmu
        if abs(dd) / denominator < 1.0e-8 or dd <= 0 or qll >= 0:
            if qmm < -TINY:
                fmu = -qm / qmm
        else:
            fmu = (qlm * ql - qll * qm) / dd
            fla = (qlm * qm - qmm * ql) / dd
        return float(fla), float(fmu)

    # -- optimizer ---------------------------------------------------------
    def optimize(self, pht, al, ps):
        """Positivity-preserving preconditioned ascent of Q = H + alpha S.

        Combines two-direction Newton acceleration with an Armijo line search.
        Converged when ||x*grad|| / (||x*gH|| + ||x*alpha*gS||) <= ac.

        Returns ``(map, iterations, Q, S, status)`` with status 0 converged,
        1 iteration limit, 2 line-search failure, 3 nonfinite objective.
        """
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return self._optimize(pht, al, ps)

    def _optimize(self, pht, al, ps):
        ac, nim = self.ac, self.nim
        self.change_maps[:] = 0
        s, dsa, chi, tem = 0.0, np.zeros(self.nol), None, None
        ph = self.fproj(ps)
        h, dha = self.like(pht, ph)
        if al > 0:
            s, dsa, chi, tem = self.entr(ps)
        dps0 = None
        for nit in range(nim + 1):
            self._check_deadline()
            hs = h + al * s
            self.objective = hs
            if not np.isfinite(hs):
                return ps, nit, hs, s, 3
            grad = dha + al * dsa
            dps = ps * grad
            norm = float(np.sqrt(np.sum(dps ** 2)))
            denominator = max(float(np.sqrt(np.sum((ps * dha) ** 2))
                                    + np.sqrt(np.sum((ps * al * dsa) ** 2))), 1.0e-12)
            delta = norm / denominator
            self.stationarity = delta
            if nit == 0 or nit % 10 == 0 or delta <= ac:
                self._log(f" iteration objective stationarity {nit} {hs!r} {delta!r}")
            if delta <= ac or norm <= 1.0e-12:
                return ps, nit, hs, s, 0
            if nit == nim:
                return ps, nit, hs, s, 1
            if nit == 0:
                dps0 = dps
            fla, fmu = self.lamun(pht, ps, dps0, dps, ph, chi, tem, al)
            up = fla * dps0 + fmu * dps
            slope = float(grad @ up)
            if not np.isfinite(slope) or slope <= 0:
                up = dps
                slope = float(grad @ up)
            step = 1.0
            shrinking = up < 0
            if shrinking.any():
                step = min(step, float(np.min(-0.99 * ps[shrinking] / up[shrinking])))
            accepted = False
            for _ in range(60):
                self._check_deadline()
                trial = ps + step * up
                ph_t = self.fproj(trial)
                hnew, dha_t = self.like(pht, ph_t)
                snew, dsa_t, chi_t, tem_t = 0.0, dsa, chi, tem
                if al > 0:
                    snew, dsa_t, chi_t, tem_t = self.entr(trial)
                qnew = hnew + al * snew
                if np.isfinite(qnew) and qnew >= hs + 1.0e-4 * step * slope:
                    accepted = True
                    break
                step /= 2
            if not accepted:
                return ps, nit, hs, s, 2
            # Reuse the accepted trial's projection, likelihood and entropy.
            ps, ph, dha, dsa, h, s, chi, tem = trial, ph_t, dha_t, dsa_t, hnew, snew, chi_t, tem_t
            dps0 = step * up
            self.change_maps[nit % 4] = grad
        raise AssertionError("The optimizer loop always returns.")

    def rrms(self, dat, model):
        """RMS statistic: sqrt(r^T W r / rank) over unmasked spectra."""
        _, quadratic = self.apply_precision(dat - model)
        return float(np.sqrt(max(quadratic, 0.0) / self.statistic_rank))

    # -- regularization schedule ------------------------------------------
    def alpha_schedule(self, dat, al0, nalpha, vmo):
        """Decrease alpha from ``al0`` by factors ``alf`` until the RMS meets ``clim``.

        A bracketed crossing is refined by bisection in log(alpha), so the
        accepted map is close to the most regularized one meeting the target.
        The schedule stops early (ier=1, plateau=1) when the statistic no
        longer improves. ``ier``: 0 success; 1 target not met; 10+k optimizer
        failure k; 19 nonfinite statistic.
        """
        alf, clim = self.alf, self.clim
        out = dict(ier=1, bracketed=0, plateau=0, nit_total=0, nit_last=0, steps=0,
                   rr=np.inf, alx=al0, sx=0.0, vmo=vmo, dato=None)
        self.rms_floor_estimate = -1.0
        trial = dict(vm=vmo.copy(), model=None, rr=np.inf, nit=0, s=0.0)

        def solve(alpha):
            trial["vm"], trial["nit"], hs, trial["s"], status = self.optimize(dat, alpha, trial["vm"])
            out["steps"] += 1
            out["nit_total"] += trial["nit"]
            trial["model"] = self.fproj(trial["vm"])
            trial["rr"] = self.rrms(dat, trial["model"])
            self._log(f" ni, al, hs, rr: {trial['nit']} {alpha!r} {hs!r} {trial['rr']!r} {clim!r}")
            if status != 0 or not np.isfinite(trial["rr"]):
                out["ier"] = 10 + status if status != 0 else 19
                self._log(f" optimizer failure {out['ier']}")
                status = max(status, 1)
            return status

        def accept(alpha):
            out.update(vmo=trial["vm"].copy(), dato=trial["model"].copy(), alx=alpha,
                       rr=trial["rr"], sx=trial["s"], nit_last=trial["nit"])

        rr_prev, drop_prev, ratio_prev = np.inf, -1.0, -1.0
        al_hi, vm_hi = al0, trial["vm"].copy()
        for ial in range(1, nalpha + 1):
            al = al0 / alf ** (ial - 1)
            if solve(al) != 0:
                # Report the failing solve; no map is accepted from it.
                out.update(rr=trial["rr"], alx=al, nit_last=trial["nit"])
                return out
            if trial["rr"] <= clim:
                accept(al)
                out["ier"] = 0
                if ial > 1:
                    out["bracketed"] = 1
                    al_lo = al
                    for _ in range(8):
                        if clim - out["rr"] <= 0.005 * clim:
                            break
                        al_mid = float(np.sqrt(al_hi * al_lo))
                        trial["vm"] = vm_hi.copy()
                        if solve(al_mid) != 0:
                            out["ier"] = 0
                            break
                        if trial["rr"] <= clim:
                            accept(al_mid)
                            al_lo = al_mid
                        else:
                            vm_hi, al_hi = trial["vm"].copy(), al_mid
                return out
            if ial > 1:
                drop = rr_prev - trial["rr"]
                # Stop when the statistic has stopped improving, or when geometric
                # extrapolation of steadily shrinking improvements (a lower bound in
                # practice) shows the remaining schedule cannot reach the target.
                if drop < 1.0e-3 * rr_prev:
                    self.rms_floor_estimate = trial["rr"]
                    accept(al)
                    out.update(plateau=1, ier=1)
                    return out
                if drop_prev > 0:
                    ratio = drop / drop_prev
                    # Two consecutive shrinking steps are required, and the slower of
                    # the two decay rates is used, so the extrapolated floor errs low.
                    if 0 < ratio < 0.8 and 0 < ratio_prev < 0.8:
                        slower = max(ratio, ratio_prev)
                        self.rms_floor_estimate = trial["rr"] - drop * slower / (1 - slower)
                        if self.rms_floor_estimate > clim * (1 + 1.0e-3):
                            accept(al)
                            out.update(plateau=1, ier=1)
                            return out
                    ratio_prev = ratio
                drop_prev = drop
            accept(al)
            out["ier"] = 1
            vm_hi, al_hi, rr_prev = trial["vm"].copy(), al, trial["rr"]
        return out

    # -- driver ------------------------------------------------------------
    def solve(self) -> SolverResult:
        """Build the forward model, run the alpha schedule and scale the outputs."""
        nph, nvp = self.nph, self.nvp
        self._log(" Floating-default maximum entropy Doppler tomography")
        self._log(f" ih, iw     {self.ih} {self.iw}")
        self._log(f" ns         {self.smoothing_width}")
        self._log(f" ac, nim    {self.ac!r} {self.nim}")
        self._log(f" al0 alf nal {self.al0!r} {self.alf!r} {self.nal}")
        self._log(f" clim       {self.clim!r}")
        self._log(f" wid, af    {self.wid!r} {self.af!r}")
        self._log(f" map size   {self.nv}")
        started = time.perf_counter()
        self.build_projection()
        self._log(f" time for geometry {time.perf_counter() - started:10.2f} s"
                  f"   nonzeros {self.matrix.nnz:12d}")
        dat = self.dat.ravel().copy()
        datb, vmb = self.baseline(dat)
        dat += datb
        nvalid = int(self.phase_valid.sum())
        sdat = float(dat.reshape(nph, nvp)[self.phase_valid].sum()) * nph / nvalid
        self.set_weights()
        # Automatic start: the legacy alpha 0.002 per unit channel information.
        alpha_start = self.al0 if self.al0 > 0 else 0.002 * nvp * self.information
        self._log(f" alpha start, information, baseline {alpha_start!r} "
                  f"{self.information!r} {self.baseline_level!r}")
        # Flat initial map with the data's total flux.
        vmo = np.maximum(sdat / self.nol, vmb)
        out = self.alpha_schedule(dat, alpha_start, max(1, self.nal), vmo)
        converged = int(out["ier"] == 0 and np.isfinite(out["rr"]))
        status = dict(
            converged=converged, optimizer_status=int(out["ier"]),
            fit_statistic=float(out["rr"]), fit_target=self.clim,
            iterations=int(out["nit_last"]), total_iterations=int(out["nit_total"]),
            alpha=float(out["alx"]), alpha_start=float(alpha_start),
            alpha_steps=int(out["steps"]), target_bracketed=int(out["bracketed"]),
            rms_plateau=int(out["plateau"]), rms_floor_estimate=float(self.rms_floor_estimate),
            entropy=float(out["sx"]), objective=float(self.objective),
            stationarity=float(self.stationarity), flux_scale=self.flux_scale,
            likelihood=self.ih, weight_mode=self.iw, statistic_rank=self.statistic_rank,
            information=float(self.information), baseline=float(self.baseline_level),
            baseline_raised=int(self.baseline_raised),
        )
        result = SolverResult(status=status, normalization=self.normalization.copy(), va=self.va)
        if not converged:
            return result
        vmo = out["vmo"]
        if not np.isfinite(vmo).all():
            raise SolverError("The converged map contains nonfinite values.")
        # Identification number (historical), in the original sampling order.
        fortran_order = dat.reshape(nph, nvp).T.ravel()
        s1 = sum(vmo[::max(1, self.nol // 30)].tolist())
        s2 = sum(fortran_order[::max(1, self.nil // 30)].tolist())
        result.identification = float(np.sin(100 * (s1 / s2))) if abs(s2) > TINY else 0.0
        # Map density in input spectral-density units per km/s: its velocity-area
        # integral equals the line integral over radial velocity in km/s.
        channel = (self.vp[-1] - self.vp[0]) / (nvp - 1)
        image = np.zeros((self.nv, self.nv))
        image[self.indx, self.indy] = ((vmo - vmb) * self.flux_scale
                                       / self.dvmap ** 2 * (channel * 1.0e5) / nph)
        # Restore spectra to their original physical flux scale for diagnostics.
        scale = (self.flux_scale / self.normalization)[:, None]
        change_maps = np.zeros((4, self.nv, self.nv))
        change_maps[:, self.indx, self.indy] = self.change_maps
        result.image = image
        result.data = (dat - datb).reshape(nph, nvp) * scale
        result.model = (out["dato"] - datb).reshape(nph, nvp) * scale
        result.change_maps = change_maps
        return result


def solve(phases, widths, velocity, flux, sigma=None, precision=None, ranks=None,
          **settings) -> SolverResult:
    """Invert one data set; see :class:`Tomography` for the settings."""
    return Tomography(phases, widths, velocity, flux, sigma, precision, ranks,
                      **settings).solve()
