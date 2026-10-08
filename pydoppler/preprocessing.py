"""Linear spectral preprocessing with compact, complete noise covariance.

Spectra are treated as samples of flux density, not integrated pixel counts.
Linear interpolation therefore preserves the sampled function; it is not a
flux-conserving rebinning operation for integrated-count inputs.
"""

from dataclasses import dataclass
import hashlib

import numpy as np


def array_digest(values):
    """SHA-256 of shape and canonical little-endian float64 numerical values."""
    array = np.ascontiguousarray(values, dtype="<f8")
    digest = hashlib.sha256(str(array.shape).encode("ascii"))
    if array.size:
        digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def covariance_digest(covariance):
    """Hash the stored covariance representation without materializing it."""
    if covariance is None:
        return None
    if isinstance(covariance, SpectralCovariance):
        fields = (covariance.variance, covariance.sampling.indices,
                  covariance.sampling.weights, covariance.design, covariance.estimator)
    else:
        fields = (covariance.matrix,)
    content = type(covariance).__name__ + ":" + ":".join(array_digest(x) for x in fields)
    return hashlib.sha256(content.encode("ascii")).hexdigest()


@dataclass
class SamplingOperator:
    """Sparse row-wise linear mapping, with only a few inputs per output."""

    indices: np.ndarray
    weights: np.ndarray
    n_input: int

    def apply(self, values):
        values = np.asarray(values)
        if values.ndim == 1:
            return np.sum(values[self.indices] * self.weights, axis=1)
        return np.sum(values[self.indices] * self.weights[:, :, None], axis=1)

    def left_multiply(self, matrix):
        """Return ``matrix @ self`` without constructing a dense sampling map."""
        matrix = np.asarray(matrix, dtype=float)
        result = np.zeros((matrix.shape[0], self.n_input))
        for row in range(matrix.shape[0]):
            np.add.at(
                result[row], self.indices.ravel(),
                (matrix[row, :, None] * self.weights).ravel(),
            )
        return result

    def compose(self, previous):
        indices = previous.indices[self.indices].reshape(len(self.indices), -1)
        weights = (
            self.weights[:, :, None] * previous.weights[self.indices]
        ).reshape(len(self.indices), -1)
        return SamplingOperator(indices, weights, previous.n_input)


def interpolation_operator(source, target):
    """Linear interpolation inside a finite, strictly increasing native grid."""
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if (source.ndim != 1 or len(source) < 2 or
            not np.all(np.isfinite(source)) or np.any(np.diff(source) <= 0)):
        raise ValueError("Wavelength grid must be finite and strictly increasing.")
    if target.ndim != 1 or not np.all(np.isfinite(target)):
        raise ValueError("Target wavelength grid must be finite and one-dimensional.")
    tol = 16 * np.finfo(float).eps * max(1., np.max(np.abs(source)))
    if np.any(target < source[0] - tol) or np.any(target > source[-1] + tol):
        raise ValueError("Spectral interpolation would extrapolate beyond input data.")
    target = np.clip(target, source[0], source[-1])
    left = np.clip(np.searchsorted(source, target, side="right") - 1, 0, len(source) - 2)
    fraction = (target - source[left]) / (source[left + 1] - source[left])
    return SamplingOperator(
        np.column_stack((left, left + 1)),
        np.column_stack((1. - fraction, fraction)), len(source),
    )


class SpectralCovariance:
    """Covariance of ``(R - X B) f`` for independent native-pixel noise.

    ``R`` is a compact sampling operator. ``X B`` contains the continuum-fit
    correction. Storage is linear in spectrum length for fixed polynomial
    degree. Use :meth:`to_dense` only when a dense covariance is needed.
    """

    def __init__(self, variance, sampling=None, design=None, estimator=None):
        self.variance = np.asarray(variance, dtype=float).copy()
        if (self.variance.ndim != 1 or not np.all(np.isfinite(self.variance))
                or np.any(self.variance <= 0)):
            raise ValueError("Native noise variances must be finite and positive.")
        size = len(self.variance)
        self.sampling = sampling or SamplingOperator(
            np.arange(size)[:, None], np.ones((size, 1)), size,
        )
        self.design = (np.zeros((len(self.sampling.indices), 0)) if design is None
                       else np.asarray(design, dtype=float))
        self.estimator = (np.zeros((0, size)) if estimator is None
                          else np.asarray(estimator, dtype=float))
        self.shape = (len(self.sampling.indices),) * 2

    def _terms(self, rows):
        indices = self.sampling.indices[rows]
        weights = self.sampling.weights[rows]
        design = self.design[rows]
        cross = np.sum(
            weights[:, :, None] * self.variance[indices, None]
            * self.estimator[:, indices].transpose(1, 2, 0), axis=1,
        )
        coefficient_cov = (self.estimator * self.variance) @ self.estimator.T
        return indices, weights, design, cross, coefficient_cov

    def diagonal(self):
        indices, weights, design, cross, coefficient_cov = self._terms(slice(None))
        # Composed resampling can refer to the same native pixel more than once.
        diagonal = np.zeros(len(indices))
        for j in range(indices.shape[1]):
            for k in range(indices.shape[1]):
                diagonal += (indices[:, j] == indices[:, k]) * (
                    weights[:, j] * weights[:, k] * self.variance[indices[:, j]]
                )
        diagonal -= 2. * np.sum(cross * design, axis=1)
        diagonal += np.sum((design @ coefficient_cov) * design, axis=1)
        return np.maximum(diagonal, 0.)

    def submatrix(self, rows):
        """Materialize covariance for selected output-pixel indices."""
        indices, weights, design, cross, coefficient_cov = self._terms(rows)
        result = np.zeros((len(indices), len(indices)))
        for j in range(indices.shape[1]):
            for k in range(indices.shape[1]):
                result += (indices[:, j, None] == indices[None, :, k]) * (
                    (weights[:, j] * self.variance[indices[:, j]])[:, None]
                    * weights[None, :, k]
                )
        result -= cross @ design.T + design @ cross.T
        result += design @ coefficient_cov @ design.T
        return (result + result.T) / 2.

    def to_dense(self):
        return self.submatrix(slice(None))

    def __array__(self, dtype=None, copy=None):
        return np.asarray(self.to_dense(), dtype=dtype)

    def transform(self, sampling, design=None, estimator=None):
        """Propagate another resampling and optional continuum subtraction."""
        composed = sampling.compose(self.sampling)
        old_design = sampling.apply(self.design)
        if design is None:
            return SpectralCovariance(
                self.variance, composed, old_design, self.estimator,
            )
        mapped_estimator = self.sampling.left_multiply(estimator)
        mapped_estimator -= (estimator @ self.design) @ self.estimator
        return SpectralCovariance(
            self.variance, composed, np.column_stack((old_design, design)),
            np.vstack((self.estimator, mapped_estimator)),
        )


class DenseSpectralCovariance:
    """Explicit user-supplied covariance on the current input wavelength grid."""

    def __init__(self, matrix):
        matrix = np.asarray(matrix, dtype=float)
        if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
                or not np.all(np.isfinite(matrix)) or np.any(np.diag(matrix) <= 0)):
            raise ValueError("Covariance must be finite, square, and have positive diagonal.")
        scale = float(np.max(np.diag(matrix)))
        if not np.allclose(matrix / scale, matrix.T / scale, atol=1e-12, rtol=1e-10):
            raise ValueError("Covariance must be symmetric.")
        if np.min(np.linalg.eigvalsh(matrix / scale)) < -1e-10:
            raise ValueError("Covariance must be positive semidefinite.")
        self.matrix = (matrix + matrix.T) / 2.
        self.shape = matrix.shape

    def diagonal(self):
        return np.diag(self.matrix).copy()

    def submatrix(self, rows):
        indices = np.arange(self.shape[0])[rows]
        return self.matrix[np.ix_(indices, indices)]

    def to_dense(self):
        return self.matrix.copy()

    def __array__(self, dtype=None, copy=None):
        return np.array(self.matrix, dtype=dtype, copy=True)

    def transform(self, sampling, design=None, estimator=None):
        operator = sampling.left_multiply(np.eye(len(sampling.indices)))
        if design is not None:
            operator -= design @ estimator
        return DenseSpectralCovariance(operator @ self.matrix @ operator.T)


def continuum_mask(wave, bands):
    return (((wave > bands[0]) & (wave < bands[1])) |
            ((wave > bands[2]) & (wave < bands[3])))


def _continuum_design(wave, target, bands, degree):
    mask = continuum_mask(wave, bands)
    if (not np.any((wave > bands[0]) & (wave < bands[1])) or
            not np.any((wave > bands[2]) & (wave < bands[3]))):
        raise ValueError("Each observed-frame continuum window needs input samples.")
    if np.count_nonzero(mask) < max(degree + 1, 3):
        raise ValueError("Continuum windows contain too few samples for the polynomial fit.")
    origin = (bands[0] + bands[3]) / 2.
    scale = (bands[3] - bands[0]) / 2.
    fit_design = np.polynomial.polynomial.polyvander((wave[mask] - origin) / scale, degree)
    out_design = np.polynomial.polynomial.polyvander((target - origin) / scale, degree)
    return mask, fit_design, out_design


def estimate_continuum_sigma(wave, flux, bands, degree):
    """Estimate homoscedastic noise with the continuum fit's residual DOF."""
    mask, design, _ = _continuum_design(wave, wave, bands, degree)
    coefficients, _, rank, _ = np.linalg.lstsq(design, flux[mask], rcond=None)
    dof = int(np.count_nonzero(mask)) - rank
    if dof <= 0:
        raise ValueError("Cannot estimate noise without continuum residual degrees of freedom.")
    residual = flux[mask] - design @ coefficients
    sigma = float(np.linalg.norm(residual) / np.sqrt(dof))
    rounding = 32 * np.finfo(float).eps * np.max(np.abs(flux[mask]))
    if not np.isfinite(sigma) or sigma <= rounding:
        raise ValueError(
            "Continuum residuals contain no measurable noise; provide positive flux errors."
        )
    return sigma


def subtract_continuum_and_resample(wave, flux, target, bands, degree, covariance=None):
    """Fit in observed wavelength coordinates and propagate the whole operator."""
    mask, fit_design, out_design = _continuum_design(wave, target, bands, degree)
    rows = np.flatnonzero(mask)
    independent = (
        isinstance(covariance, SpectralCovariance) and covariance.design.shape[1] == 0
        and covariance.sampling.indices.shape == (len(wave), 1)
        and np.array_equal(covariance.sampling.indices[:, 0], np.arange(len(wave)))
    )
    if covariance is None or independent:
        weights = (np.ones(len(rows)) if covariance is None
                   else 1. / np.sqrt(covariance.diagonal()[rows]))
        weights /= np.max(weights)
        weighted_design = fit_design * weights[:, None]
        fit_operator = np.linalg.pinv(weighted_design) * weights[None, :]
    else:
        continuum_cov = covariance.submatrix(rows)
        scale = float(np.max(np.diag(continuum_cov)))
        eigenvalues, eigenvectors = np.linalg.eigh(continuum_cov / scale)
        retained = eigenvalues > max(1e-12, np.max(eigenvalues) * 1e-12)
        whitening = eigenvectors[:, retained].T / np.sqrt(eigenvalues[retained, None])
        weighted_design = whitening @ fit_design
        fit_operator = np.linalg.pinv(weighted_design) @ whitening
    if np.linalg.matrix_rank(weighted_design) < degree + 1:
        raise ValueError("Continuum fit is not identifiable with the supplied noise covariance.")
    coefficients = fit_operator @ flux[mask]
    sampling = interpolation_operator(wave, target)
    output = sampling.apply(flux) - out_design @ coefficients
    if covariance is None:
        return output, None
    estimator = np.zeros((degree + 1, len(wave)))
    estimator[:, mask] = fit_operator
    return output, covariance.transform(sampling, out_design, estimator)
