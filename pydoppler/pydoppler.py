#! Python PyDoppler

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np

from importlib.resources import files as resource_files


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


_MATPLOTLIB_CONFIGURED = False


def _lazy_import_matplotlib_pyplot():
    """Import Matplotlib lazily so importing :mod:`pydoppler` stays headless-safe."""

    global _MATPLOTLIB_CONFIGURED
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - depends on local environment
        raise ImportError(
            "Matplotlib is required for plotting functionality. "
            "Install it with `pip install matplotlib`."
        ) from exc

    if not _MATPLOTLIB_CONFIGURED:
        mpl.rcParams.update({"font.size": 12})
        _MATPLOTLIB_CONFIGURED = True

    return mpl, plt


def _copy_tree(source, destination: Path, overwrite: bool = False) -> List[Path]:
    """Copy the files contained in *source* into *destination*.

    Parameters
    ----------
    source:
        An importlib ``Traversable`` pointing to a directory.
    destination:
        Filesystem directory that will receive the copied files.
    overwrite:
        When ``True`` existing files will be replaced, otherwise they are kept
        untouched.
    """

    dest_path = Path(destination)
    dest_path.mkdir(parents=True, exist_ok=True)

    copied: List[Path] = []

    def _walk(node, relative: Path = Path(".")) -> None:
        for item in node.iterdir():
            item_rel = relative / item.name
            if item.is_dir():
                _walk(item, item_rel)
                continue

            target = dest_path / item_rel
            if target.exists() and not overwrite:
                continue

            target.parent.mkdir(parents=True, exist_ok=True)
            with item.open("rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            copied.append(target)

    _walk(source)

    return copied


@dataclass
class DopplerOutput:
    """Raw solver arrays and grids.

    ``vp``, ``va`` and ``vmap`` are in km/s, phases and exposure widths in
    cycles; ``gamma`` retains the input header's cm/s units. ``im`` uses
    displayed [Vy, Vx] array order after the documented Fortran-axis exchange.
    """
    path: Path
    nph: int
    nvp: int
    nv: int
    w0: float
    aa: float
    gamma: float
    pha: np.ndarray
    dpha: np.ndarray
    vp: np.ndarray
    va: float
    vmap: np.ndarray
    dm: np.ndarray
    im: np.ndarray
    dmr: np.ndarray
    dpx: np.ndarray


def _scale_by_absmax(data: np.ndarray) -> np.ndarray:
    """Scale *data* by its finite absolute maximum, preserving NaNs."""

    arr = np.asarray(data, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.full_like(arr, np.nan, dtype=float)

    scale = float(np.nanmax(np.abs(arr[finite])))
    out = np.full_like(arr, np.nan, dtype=float)
    if not np.isfinite(scale) or scale <= 0:
        out[finite] = 0.0
        return out

    out[finite] = arr[finite] / scale
    return out


def _scale_by_peak(data: np.ndarray) -> np.ndarray:
    """Divide *data* by its positive peak, so zero flux stays at zero.

    Negative values remain negative. Maps without positive values are scaled
    by their absolute maximum instead.
    """

    arr = np.asarray(data, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.full_like(arr, np.nan, dtype=float)
    peak = float(np.max(arr[finite]))
    if not np.isfinite(peak) or peak <= 0:
        return _scale_by_absmax(arr)
    out = np.full_like(arr, np.nan, dtype=float)
    out[finite] = arr[finite] / peak
    return out


def _auto_display_limits(data: np.ndarray) -> List[float]:
    """Robust 5th-95th percentile display limits from finite samples."""

    arr = np.asarray(data, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return [0.0, 1.0]
    lo = float(np.percentile(finite, 5.0))
    hi = float(np.percentile(finite, 95.0))
    if hi <= lo:
        vmax = float(np.max(np.abs(finite)))
        if vmax <= 0:
            return [0.0, 1.0]
        lo, hi = -vmax, vmax
    return [lo, hi]


def _set_mappable_norm(mappable, norm) -> None:
    """Apply a Matplotlib norm to a mappable across Matplotlib versions."""

    if hasattr(mappable, "set_norm"):
        mappable.set_norm(norm)
    else:  # pragma: no cover - defensive fallback for unusual mappables
        mappable.norm = norm


def _bin_edges(centers: np.ndarray) -> np.ndarray:
    """Return midpoint bin edges from finite, strictly increasing centers."""
    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 1 or centers.size < 2 or not np.all(np.isfinite(centers)):
        raise ValueError("At least two finite one-dimensional bin centers are required.")
    if np.any(np.diff(centers) <= 0):
        raise ValueError("Bin centers must be strictly increasing.")
    return np.concatenate(([centers[0] - (centers[1] - centers[0]) / 2],
                           (centers[:-1] + centers[1:]) / 2,
                           [centers[-1] + (centers[-1] - centers[-2]) / 2]))


def _read_verified_dopout(path: Union[Path, str]) -> DopplerOutput:
    """Reject stale/failed managed outputs before displaying their diagnostics."""
    path = Path(path)
    manifest_path = path.parent / "run.json"
    if manifest_path.exists():
        import hashlib
        import json
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
        except (OSError, ValueError) as exc:
            raise RuntimeError(f"Cannot verify managed output {path}.") from exc
        if (not isinstance(manifest, dict) or manifest.get("success") is not True
                or manifest.get("output_sha256") != digest):
            raise RuntimeError(
                f"Refusing previous, failed or unverified output {path}; "
                "rerun Syncdop successfully before plotting."
            )
    return _parse_dopout(path)


def _parse_dopout(dopout_path: Union[Path, str]) -> DopplerOutput:
    """Parse Spruit's ``dop.out`` output into structured NumPy arrays."""

    path = Path(dopout_path)
    with path.open(encoding="utf-8") as handle:
        lines = handle.readlines()

    if len(lines) < 2:
        raise RuntimeError(f"{path} is incomplete: expected at least two header lines.")

    try:
        head0 = lines[0].split()
        nph = int(head0[0])
        nvp = int(head0[1])
        nv_header = int(head0[2])
        w0 = float(head0[3])
        aa = float(head0[4])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"Failed to parse the first header line of {path}.") from exc

    try:
        head1 = lines[1].split()
        gamma = float(head1[0])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"Failed to parse the second header line of {path}.") from exc

    tokens = [
        token.replace("D", "e").replace("E", "e")
        for token in "".join(lines[2:]).split()
    ]
    index = 0

    def take(count: int, label: str) -> List[str]:
        nonlocal index
        end = index + count
        if end > len(tokens):
            raise RuntimeError(f"{path} ended unexpectedly while reading {label}.")
        chunk = tokens[index:end]
        index = end
        return chunk

    def take_float_array(count: int, label: str) -> np.ndarray:
        try:
            return np.asarray([float(value) for value in take(count, label)], dtype=float)
        except ValueError as exc:
            raise RuntimeError(f"Failed to parse {label} in {path}.") from exc

    if nph < 1 or nvp < 2 or nv_header < 2:
        raise RuntimeError(f"Invalid array dimensions in {path}.")
    pha = take_float_array(nph, "phase grid") / (2.0 * np.pi)
    iph = int(take(1, "exposure-width flag")[0])
    dpha = (take_float_array(nph, "phase widths") / (2.0 * np.pi)
            if iph == 1 else np.full(nph, 1.0 / nph))
    # Fortran writes spectral bin centers, already corrected for gamma.
    vp = take_float_array(nvp, "velocity grid") / 1e5
    if not np.all(np.isfinite(vp)) or np.any(np.diff(vp) <= 0):
        raise RuntimeError(f"Velocity centers must be finite and increasing in {path}.")

    dm = take_float_array(nvp * nph, "input trail").reshape(nvp, nph)

    try:
        params = take(14, "inversion parameters")
        nv = int(params[11])
        va = float(params[12]) / 1e5
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"Failed to parse inversion parameters in {path}.") from exc

    if nv < 2 or not np.isfinite(va) or va <= 0:
        raise RuntimeError(f"Invalid Doppler map grid in {path}.")
    vmap = np.linspace(-va, va, nv)

    im = take_float_array(nv * nv, "Doppler map").reshape(nv, nv)
    take(3, "reconstruction header")
    dmr = take_float_array(nvp * nph, "reconstructed trail").reshape(nvp, nph)
    take(4, "projection header")
    dpx = take_float_array(nv * nv, "projection map").reshape(nv, nv)

    if nv != nv_header:
        LOGGER.debug(
            "dop.out header nv=%s differs from payload nv=%s in %s", nv_header, nv, path
        )

    return DopplerOutput(
        path=path,
        nph=nph,
        nvp=nvp,
        nv=nv,
        w0=w0,
        aa=aa,
        gamma=gamma,
        pha=pha,
        dpha=dpha,
        vp=vp,
        va=va,
        vmap=vmap,
        dm=dm,
        im=im,
        dmr=dmr,
        dpx=dpx,
    )


def get_fortran_code_path():
    """Return a Traversable pointing to the bundled Fortran sources."""

    return resource_files(__package__).joinpath("fortran_code")


def get_test_data_path():
    """Return a Traversable pointing to the bundled test data."""

    return resource_files(__package__).joinpath("test_data")


def copy_fortran_code(
    destination: Union[Path, str], overwrite: bool = False
) -> List[Path]:
    """Copy the bundled Fortran assets into *destination*."""

    return _copy_tree(get_fortran_code_path(), Path(destination), overwrite)


def copy_test_data(
    destination: Union[Path, str], overwrite: bool = False
) -> List[Path]:
    """Copy the bundled sample data into *destination*."""

    return _copy_tree(get_test_data_path(), Path(destination), overwrite)


def install_sample_script(
    destination: Union[Path, str],
    overwrite: bool = False,
) -> Path:
    """Install the example ``sample_script.py`` in *destination*.

    When ``overwrite`` is ``False`` the script will not clobber an existing
    file. Instead it will fall back to the ``sample_script-<n>.py`` naming
    convention previously used by the project.
    """

    dest_dir = Path(destination)
    dest_dir.mkdir(parents=True, exist_ok=True)

    target = dest_dir / "sample_script.py"
    if target.exists() and not overwrite:
        index = 1
        while target.exists():
            target = dest_dir / f"sample_script-{index}.py"
            index += 1

    source = get_test_data_path().joinpath("sample_script.py")
    with source.open("rb") as src, target.open("wb") as dst:
        shutil.copyfileobj(src, dst)

    return target


class spruit:
    """Wrapper around Henk Spruit's Doppler tomography code.

    Typical workflow
    ----------------
    1. Set ``base_dir`` and ``list`` to point at the spectra and phase file.
    2. Call :meth:`Foldspec` to load and interpolate spectra to a common grid.
    3. Call :meth:`Dopin` to continuum-subtract and write ``dopin`` in ``workdir``.
    4. Call :meth:`Syncdop` to compile/run the Fortran code in ``workdir``.
    5. Use :meth:`Dopmap`, :meth:`Reco`, and :meth:`Residuals` for diagnostics.

    Notes
    -----
    - ``workdir`` isolates all Fortran build products and outputs (defaults to
      ``./pydoppler-workdir``).
    - ``nbins`` controls the *diagnostic* phase binning used in trail plots; it
      does not change the inversion, which uses all input spectra.
    """
    def __init__(
        self,
        force_install: bool = False,
        auto_install: bool = False,
        interactive: bool = False,
        workdir: Optional[Union[Path, str]] = None,
    ):
        default_workdir = Path.cwd() / "pydoppler-workdir"
        self.workdir = default_workdir if workdir is None else Path(workdir)
        self.workdir = self.workdir.expanduser().resolve()
        self._workdir_base = self.workdir
        self.object = 'disc'
        self.wave = 0.0
        self.flux = 0.0
        self.pha = 0.0
        self.input_files = 0.0
        self.input_phase = 0.0
        self.input_dpha = None
        self.trsp = 0.0
        self.nbins = 20
        self.normalised_flux = 0.0
        self.normalised_wave = 0.0
        self.base_dir = '.'
        self.lam0 = 6562.83
        self.delw = 80
        self.list = 'phases.txt'
        self.gama = 0.0
        self.delta_phase = 0.001
        self.flux_err = None
        self.normalised_flux_err = None
        self.normalised_flux_covariance = None

        self.verbose = True
        self.interactive = interactive
        self.logger = logging.getLogger(f"{__name__}.spruit")
        if not self.logger.handlers:
            self.logger.addHandler(logging.NullHandler())

        # Plotting and map-grid parameters
        self.plot = True  # Produce diagnostic figures by default
        self.overs = 0.4  # Map pixels per axis relative to spectral channels

        ####### Dop.in parameters
        self.ih = 1  # 1: least squares; 0: legacy positive log objective
        self.iw = 0
        self.pb0 = 0.0
        self.pb1 = 0.0
        self.ns = 7
        self.ac = 8e-4
        self.nim = 1000
        self.al0 = None  # None: automatic (0.002 per unit channel information)
        self.alf = 1.7
        self.nal = None  # None: 0 for unweighted fits, 25 for Gaussian fits
        self.clim = None  # None: 1.0
        self.ipri = 2
        self.norm = 0
        self.conditional_normalization = False
        self.wid = 10e5
        self.af = 0.0
        self.instrument_sigma = 0.0  # Gaussian LSF sigma in km/s, excluding pixel integration
        self.phase_subsamples = 0  # 0 selects adaptive exposure quadrature
        self.baseline_fraction = 0.1
        self.covariance_mode = "full"
        self.covariance_rcond = 1e-10
        self.covariance_ranks = None
        self.map_size = None  # otherwise use overs times the spectral channel count
        self.run_timeout = 600.0
        self.max_memory_bytes = 2_000_000_000
        self.fit_diagnostics = None
        self.valid_phase_mask = None


        if auto_install:
            self.install_assets(self.workdir, force_install)

    # ------------------------------------------------------------------
    # helper methods
    def _log(self, level: int, message: str) -> None:
        """Emit *message* through :mod:`logging` and optionally echo it."""

        self.logger.log(level, message)
        if self.verbose and level >= logging.INFO:
            print(message)

    def install_assets(
        self, destination: Union[Path, str], force: bool = False
    ) -> None:
        """Install bundled assets into *destination*.

        Parameters
        ----------
        destination:
            Directory where the Fortran sources and helper script should be
            installed.
        force:
            When ``True`` existing files are overwritten.
        """

        dest_dir = Path(destination)
        dest_dir.mkdir(parents=True, exist_ok=True)

        copied_fortran = copy_fortran_code(dest_dir, overwrite=force)
        script_path: Optional[Path] = None
        if force or copied_fortran:
            script_path = install_sample_script(dest_dir, overwrite=force)

        if copied_fortran:
            self._log(logging.INFO, f"Installed Fortran assets into {dest_dir}.")
        else:
            self._log(logging.DEBUG, "Fortran assets already present; skipping copy.")

        if script_path:
            self._log(logging.INFO, f"Sample script available at {script_path}.")
        else:
            self._log(logging.DEBUG, "Sample script already present; not copied.")


    def set_workdir(self, destination: Union[Path, str]) -> Path:
        """Set the working directory used for Fortran execution and outputs.

        Parameters
        ----------
        destination:
            Path where ``dopin``, ``dop.in``, ``dop.out``, and intermediate build
            products should be created.

        Returns
        -------
        pathlib.Path
            The resolved working directory.
        """

        self.workdir = Path(destination).expanduser().resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self._workdir_base = self.workdir
        return self.workdir

    def make_run_dir(self, prefix: str = "pydoppler-run") -> Path:
        """Create a new run directory and switch ``workdir`` to it."""

        base = Path(self._workdir_base)
        base.mkdir(parents=True, exist_ok=True)
        index = 0
        while True:
            name = prefix if index == 0 else f"{prefix}-{index}"
            candidate = base / name
            try:
                candidate.mkdir(parents=True, exist_ok=False)
            except FileExistsError:
                index += 1
                continue
            break
        copy_fortran_code(candidate, overwrite=False)
        self.workdir = candidate.expanduser().resolve()
        return self.workdir


    def Foldspec(self):
        """Load spectra and orbital phases and interpolate onto a common grid.

        The phase file (``base_dir/list``) must contain two whitespace-separated
        columns: ``spectrum_filename  orbital_phase``. An optional third column
        is interpreted as the exposure width in orbital phase and is stored in
        ``input_dpha`` (use :meth:`Dopin` with ``use_list_dpha=True`` to write it).

        Individual spectrum files must contain at least two columns (wavelength,
        flux). A third column is interpreted as 1σ flux uncertainty and will be
        propagated through the interpolation step. ``flux_covariance`` retains
        the induced pixel covariance as compact objects with ``to_dense()``.
        Native spectra are retained so :meth:`Dopin` can resample only once.
        """
        from .preprocessing import SpectralCovariance, interpolation_operator
        list_path = Path(self.base_dir) / self.list
        if not list_path.exists():
            raise FileNotFoundError(
                f"Phase file '{list_path}' is not accessible. Check 'base_dir' and 'list'. "
                "If you are trying the bundled example dataset, run `pydoppler.test_data()` "
                "or `pydoppler.copy_test_data(path)` to copy it into your working directory."
            )

        list_dir = list_path.parent

        rows = []
        with open(list_path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"Invalid phase-file row (expected >=2 columns): {line!r}")
                fname = parts[0]
                try:
                    phase = float(parts[1])
                except ValueError as exc:
                    raise ValueError(f"Invalid phase value in row: {line!r}") from exc
                if not np.isfinite(phase):
                    raise ValueError(f"Orbital phase must be finite: {line!r}")
                dpha = None
                if len(parts) >= 3:
                    try:
                        dpha = float(parts[2])
                    except ValueError as exc:
                        raise ValueError(f"Invalid exposure width value in row: {line!r}") from exc
                    if not np.isfinite(dpha) or not 0 < dpha <= 1:
                        raise ValueError(f"Exposure width must be in (0, 1] cycles: {line!r}")
                rows.append((fname, phase, dpha))

        if not rows:
            raise ValueError("No entries found in the list file.")

        files = np.asarray([row[0] for row in rows], dtype=str)
        phases = np.asarray([row[1] for row in rows], dtype=float)
        dpha_vals = [row[2] for row in rows]
        if any(val is not None for val in dpha_vals):
            input_dpha = np.asarray(
                [np.nan if val is None else float(val) for val in dpha_vals],
                dtype=float,
            )
        else:
            input_dpha = None

        def _read_spectrum(path: Path):
            arr = np.loadtxt(path, comments="#")
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)
            if arr.ndim != 2 or arr.shape[1] < 2:
                raise ValueError("Spectrum file must contain at least two columns (wavelength, flux).")
            wave = np.asarray(arr[:, 0], dtype=float)
            flux = np.asarray(arr[:, 1], dtype=float)
            err = None
            if arr.shape[1] >= 3:
                err = np.asarray(arr[:, 2], dtype=float)

            # Keep the interpolation grid well-defined for np.interp:
            # sort by wavelength and collapse duplicate coordinates.
            order = np.argsort(wave)
            wave = wave[order]
            flux = flux[order]
            if err is not None:
                err = err[order]

            if not np.all(np.isfinite(wave)) or np.any(wave <= 0):
                raise ValueError("Wavelength samples must be finite and positive.")
            if not np.all(np.isfinite(flux)):
                raise ValueError("Flux samples must be finite; mask bad pixels explicitly before loading.")
            if err is not None and np.any((~np.isfinite(err)) | (err <= 0)):
                raise ValueError("Supplied flux uncertainties must be finite and positive.")
            if wave.size < 2:
                raise ValueError("Spectrum must contain at least two finite wavelength/flux samples.")

            if np.any(np.diff(wave) <= 0):
                uniq_wave, inv = np.unique(wave, return_inverse=True)
                counts = np.bincount(inv).astype(float)
                flux_sum = np.bincount(inv, weights=flux)
                flux = flux_sum / counts
                if err is not None:
                    err2_sum = np.bincount(inv, weights=np.square(err))
                    err = np.sqrt(err2_sum) / counts
                wave = uniq_wave
            if wave.size < 2:
                raise ValueError("Spectrum needs at least two distinct wavelengths.")
            return wave, flux, err

        loaded = []
        flux_err: List[Optional[np.ndarray]] = []

        if getattr(self, "nbins", None) is None:
            phases_sorted = np.sort(np.asarray(phases, dtype=float))
            dphi = np.diff(phases_sorted)
            dphi = dphi[dphi > 0]
            if dphi.size == 0:
                self.nbins = 20
            else:
                self.nbins = int(max(4, round(1.5 / float(np.median(dphi)))))

        try:
            rep_dphi = float(np.median(np.diff(np.sort(phases))))
        except Exception:  # pragma: no cover - extremely small data set
            rep_dphi = float("nan")

        self._log(
            logging.INFO,
            f"Number of bins: {self.nbins} (median phase spacing {rep_dphi:.5f})",
        )

        for z, (fname, ph) in enumerate(zip(files, phases)):
            sp_path = list_dir / fname
            try:
                w, f, e = _read_spectrum(sp_path)
            except Exception as exc:
                raise RuntimeError(f"Failed to read spectrum '{sp_path}': {exc}") from exc

            loaded.append((fname, w, f, e))

            self._log(
                logging.INFO,
                f"{str(z + 1).zfill(3)} {fname}  {ph} {w.size}",
            )

        overlap_min = float(max(w[0] for _, w, _, _ in loaded))
        overlap_max = float(min(w[-1] for _, w, _, _ in loaded))
        if overlap_min >= overlap_max:
            raise ValueError(
                "Input spectra do not share a common wavelength overlap. "
                f"Latest starting wavelength: {overlap_min:.6f}; earliest ending "
                f"wavelength: {overlap_max:.6f}."
            )

        overlap_counts = [
            int(np.count_nonzero((w >= overlap_min) & (w <= overlap_max)))
            for _, w, _, _ in loaded
        ]
        ref_index = int(np.argmax(overlap_counts))
        ref_name, ref_wave, _, _ = loaded[ref_index]
        ref_mask = (ref_wave >= overlap_min) & (ref_wave <= overlap_max)
        wo = np.asarray(ref_wave[ref_mask], dtype=float)
        if wo.size < 2:
            raise ValueError(
                "The common wavelength overlap contains fewer than two samples on the "
                f"chosen reference grid ({ref_name})."
            )

        full_min = float(min(w[0] for _, w, _, _ in loaded))
        full_max = float(max(w[-1] for _, w, _, _ in loaded))
        if overlap_min > full_min or overlap_max < full_max:
            self._log(
                logging.INFO,
                "Trimming spectra to common wavelength overlap "
                f"[{overlap_min:.6f}, {overlap_max:.6f}] using reference grid '{ref_name}'.",
            )

        wave: List[np.ndarray] = []
        flux: List[np.ndarray] = []
        flux_covariance = []
        for _, w, f, e in loaded:
            wave.append(wo.copy())
            sampling = interpolation_operator(w, wo)
            flux.append(sampling.apply(f))
            if e is None:
                flux_err.append(None)
                flux_covariance.append(None)
            else:
                covariance = SpectralCovariance(np.square(e), sampling)
                flux_err.append(np.sqrt(covariance.diagonal()))
                flux_covariance.append(covariance)

        self.wave = wave
        self.flux = flux
        self.pha = np.asarray(phases, float)
        self.input_files = np.asarray(files, dtype=str)
        self.input_phase = np.asarray(phases, dtype=float)
        self.input_dpha = input_dpha
        self.trsp = np.vstack(flux)
        if any(err is not None for err in flux_err):
            self.flux_err = flux_err
        else:
            self.flux_err = None
        self.flux_covariance = flux_covariance
        self._foldspec_covariance_snapshot = tuple(flux_covariance)
        self._native_spectra = [(w.copy(), f.copy(), None if e is None else e.copy())
                                for _, w, f, e in loaded]
        self._foldspec_wave_snapshot = [w.copy() for w in self.wave]
        self._foldspec_flux_snapshot = np.asarray(self.flux).copy()
        self._foldspec_error_snapshot = [None if e is None else e.copy() for e in flux_err]


    def _auto_continuum_band(self) -> Sequence[float]:
        """Derive continuum windows for non-interactive use.

        The routine samples the outer 5% (at least four samples) on both
        sides of the line and returns two intervals suitable for the continuum
        polynomial fit used in :meth:`Dopin`.
        """

        if np.ndim(self.wave) == 0 or not len(self.wave) or not np.size(self.wave[0]):
            raise RuntimeError(
                "Wavelength grid is empty. Run 'Foldspec' before calling 'Dopin'."
            )

        wave = np.asarray(self.wave[0], dtype=float)
        if wave.size < 8:
            raise ValueError("Cannot determine continuum bands from fewer than eight samples.")

        window = max(4, int(round(0.05 * wave.size)))
        left = wave[:window]
        right = wave[-window:]

        return (
            float(left.min()),
            float(left.max()),
            float(right.min()),
            float(right.max()),
        )


    def Dopin(
        self,
        poly_degree=2,
        continnum_band=False,
        rebin: bool = True,
        plot_median: bool = False,
        rebin_wave=0.0,
        xlim=None,
        two_orbits=True,
        vel_space=True,
        verbose: Optional[bool] = None,
        flux_err: Optional[np.ndarray] = None,
        estimate_errors: bool = True,
        use_list_dpha: bool = False,
        plot: Optional[bool] = None,
        show: Optional[bool] = None,
        continuum_band: Optional[Sequence[float]] = None,
        flux_covariance=None,
    ):
        """Continuum-subtract spectra and write the Fortran input file ``dopin``.

        The method fits a polynomial continuum in two wavelength windows bracketing
        the emission line, subtracts it from each spectrum, interpolates onto a
        common velocity grid, and writes the result to ``workdir/dopin`` for
        Spruit's Fortran tomography code.

        Parameters
        ----------
        poly_degree:
            Polynomial degree for the continuum fit.
        continuum_band:
            Four observed-frame wavelength limits ``[x1, x2, x3, x4]`` defining two continuum
            windows ``(x1, x2)`` and ``(x3, x4)``. When omitted, windows are
            selected automatically if ``interactive`` is disabled, or chosen
            interactively when ``interactive=True``.
        continnum_band:
            Legacy alias for ``continuum_band`` (kept for backward compatibility).
        flux_err:
            Optional 1σ uncertainties with shape ``(nspec, nwave)`` matching the
            output of :meth:`Foldspec`. Explicitly supplied errors are assumed
            independent unless ``flux_covariance`` is also supplied.
        flux_covariance:
            Optional per-spectrum covariance matrices on the loaded wavelength
            grid. Compact covariance objects returned by :meth:`Foldspec` are
            also accepted. Overrides inferred native-pixel covariance.
        estimate_errors:
            When ``True`` (default) and ``iw==1``, missing/invalid uncertainties
            are estimated from continuum residuals.
        use_list_dpha:
            When ``True``, use the optional third column of the phase file
            (exposure widths in orbital phase) instead of the scalar ``delta_phase``.
        rebin:
            When ``True`` (default), phase-bin the diagnostic trail plot into
            ``nbins``; when ``False`` plot the raw spectra ordered by phase.
        verbose:
            Override ``self.verbose`` for this call. Use ``None`` (default) to
            keep the instance setting.
        plot:
            When ``True`` produce diagnostic plots; when ``False`` no Matplotlib
            is imported and only the ``dopin`` file is written.
        show:
            When ``True`` display figures (only relevant when ``plot=True``).

        Notes
        -----
        ``nbins`` only affects the diagnostic trail plot binning; the inversion
        uses all spectra/phases written to ``dopin``.
        """
        from .preprocessing import (
            DenseSpectralCovariance, SpectralCovariance, array_digest, covariance_digest,
            estimate_continuum_sigma,
            subtract_continuum_and_resample,
        )

        if np.ndim(self.wave) == 0 or len(self.wave) == 0:
            raise RuntimeError("No spectra loaded. Run 'Foldspec' before normalising.")
        wave0 = np.asarray(self.wave[0], dtype=float)
        input_flux = np.asarray(self.flux, dtype=float)
        if (wave0.ndim != 1 or wave0.size < 2 or not np.all(np.isfinite(wave0))
                or np.any(wave0 <= 0) or np.any(np.diff(wave0) <= 0)):
            raise ValueError("Loaded wavelength grid must be finite, positive, and strictly increasing.")
        if (input_flux.ndim != 2 or input_flux.shape[1] != wave0.size
                or not np.all(np.isfinite(input_flux)) or not len(input_flux)):
            raise ValueError("Loaded flux must be a finite array with shape (nspec, nwave).")
        if any(not np.array_equal(w, wave0) for w in self.wave):
            raise ValueError("Loaded spectra must share a wavelength grid; run Foldspec first.")
        phases = np.asarray(self.input_phase, dtype=float)
        if phases.shape != (len(input_flux),) or not np.all(np.isfinite(phases)):
            raise ValueError("Input phases must be finite with one value per spectrum.")
        self.pha = phases.copy()
        if (isinstance(poly_degree, bool) or not isinstance(poly_degree, (int, np.integer))
                or poly_degree < 0):
            raise ValueError("poly_degree must be a non-negative integer.")
        if not np.isfinite(self.delta_phase) or not 0 < self.delta_phase <= 1:
            raise ValueError("delta_phase must be in (0, 1] orbital cycles.")
        if not np.isfinite(self.lam0) or self.lam0 <= 0:
            raise ValueError("lam0 must be finite and positive.")
        if not np.isfinite(self.delw) or self.delw <= 0:
            raise ValueError("delw must be finite and positive.")
        if not np.isfinite(self.gama) or abs(self.gama) >= 2.99792458e5:
            raise ValueError("Systemic velocity 'gama' must be finite and satisfy |gama| < c.")
        if rebin_wave is not None and rebin_wave != 0:
            raise ValueError("rebin_wave is not implemented; use the native diagnostic wavelength sampling.")

        log = self._log
        if verbose is not None:
            verbose_flag = bool(verbose)

            def log(level: int, message: str) -> None:
                self.logger.log(level, message)
                if verbose_flag and level >= logging.INFO:
                    print(message)

        if continuum_band is not None:
            if continnum_band is not False and continnum_band is not None:
                raise TypeError(
                    "Pass only one of 'continnum_band' (legacy) or 'continuum_band'."
                )
            continnum_band = continuum_band
        has_continuum = continnum_band is not False and continnum_band is not None

        if plot is None:
            plot = bool(self.plot)
        if show is None:
            show = bool(plot)

        plt = None
        if plot or (self.interactive and not has_continuum):
            _, plt = _lazy_import_matplotlib_pyplot()
        elif self.interactive and not has_continuum:
            raise RuntimeError(
                "Interactive continuum selection requires plotting; set plot=True."
            )

        lam=self.lam0
        cl=2.99792458e5
        systemic_factor = np.sqrt((1. + self.gama / cl) / (1. - self.gama / cl))
        observed_line_center = lam * systemic_factor

        cmaps = plt.cm.binary_r if plt is not None else None  # noqa: F841

        wmin = float(np.min(self.wave[0]))
        wmax = float(np.max(self.wave[0]))
        if observed_line_center < wmin or observed_line_center > wmax:
            raise ValueError(
                f"Systemic-velocity shifted line center {observed_line_center} Å is out of "
                f"bounds. Must be between {wmin} and {wmax}."
            )


        avgspec=np.sum(self.flux,axis=0)
        continuum_fig = None
        continuum_ax = None
        if plot and plt is not None and self.interactive and not has_continuum:
            continuum_fig = plt.figure(num="Continuum Selection", figsize=(6.57, 4.0))
            plt.clf()
            continuum_ax = continuum_fig.add_subplot(111)
            continuum_ax.plot(self.wave[0], avgspec / len(self.pha))
            continuum_ax.set_xlabel(r'Wavelength / $\AA$')
            continuum_ax.set_ylabel('Input flux')

        if not has_continuum:
            if self.interactive:
                if plt is None:
                    raise RuntimeError(
                        "Interactive continuum selection requires Matplotlib; install it or pass continuum_band."
                    )
                if not plot:
                    raise RuntimeError(
                        "Interactive continuum selection requires plot=True."
                    )
                log(logging.INFO, "Choose 4 points to define the continuum.")
                xor=[]
                for _ in np.arange(4):
                    selection=plt.ginput(1,timeout=-1)
                    if not selection:
                        raise RuntimeError("Continuum selection aborted by user.")
                    xor.append(float(selection[0][0]))
                    if continuum_ax is not None and continuum_fig is not None:
                        continuum_ax.axvline(x=xor[-1], linestyle='--', color='k')
                        continuum_fig.canvas.draw_idle()
            else:
                xor=list(self._auto_continuum_band())
                message = (
                    "Using automatic continuum bands at "
                    f"[{xor[0]:.2f}, {xor[1]:.2f}] and [{xor[2]:.2f}, {xor[3]:.2f}] Å."
                )
                log(logging.INFO, message)
                if continuum_ax is not None:
                    for idx,val in enumerate(xor):
                        label='Cont Bands' if idx == 0 else ''
                        continuum_ax.axvline(x=val, linestyle='--', color='k', label=label)
        else:
            xor=[float(v) for v in continnum_band]
            if len(xor) != 4:
                raise ValueError("continuum_band must contain four wavelength limits.")
            if not np.all(np.isfinite(xor)) or np.any(np.diff(xor) <= 0):
                raise ValueError("Continuum limits must be finite, increasing, and non-overlapping.")
            if continuum_ax is not None:
                for idx,val in enumerate(xor):
                    label='Cont Bands' if idx == 0 else ''
                    continuum_ax.axvline(x=val, linestyle='--', color='k', label=label)
        lop = ((self.wave[0]>xor[0]) * (self.wave[0]<xor[1])) + ((self.wave[0]>xor[2]) * (self.wave[0]<xor[3]))
        n_cont = int(np.count_nonzero(lop))
        min_cont = max(poly_degree + 1, 3)
        if n_cont < min_cont:
            raise ValueError(
                f"Continuum windows contain {n_cont} samples, but at least {min_cont} "
                f"are required for a degree-{poly_degree} polynomial fit."
            )
        yor=avgspec[lop]/len(self.pha)
        if continuum_ax is not None:
            continuum_ax.set_ylim(avgspec[lop].min()/len(self.pha)*0.8,avgspec.max()/len(self.pha)*1.1)
            z = np.polyfit(self.wave[0][lop], yor, poly_degree)
            pz = np.poly1d(z)
            linfit = pz(self.wave[0])
            continuum_ax.plot(self.wave[0], linfit, 'r', label='Cont Fit')
            continuum_ax.legend(fontsize=10)
            continuum_ax.set_xlim(xor[0]-10, xor[3]+10)
            continuum_fig.tight_layout()

        beta = self.gama / cl
        nufac = np.sqrt((1.0 + beta) / (1.0 - beta))
        rest_wave = wave0 / nufac
        region = (rest_wave > self.lam0 - self.delw) & (rest_wave < self.lam0 + self.delw)
        if np.count_nonzero(region) < 2:
            raise ValueError(
                "No spectral samples fall inside the Doppler-map wavelength window "
                "(lam0 ± delw) after gamma correction."
            )
        ratios = rest_wave[region] / self.lam0
        velocities = cl * (ratios**2 - 1.) / (ratios**2 + 1.)
        self.vell = np.linspace(velocities[0], velocities[-1], len(velocities))
        self.normalised_wave = self.lam0 * np.sqrt((1. + self.vell / cl) / (1. - self.vell / cl))
        observed_target = self.normalised_wave * nufac

        def same_errors(current, previous):
            if current is None:
                return all(value is None for value in previous)
            if not hasattr(current, "__len__") or len(current) != len(previous):
                return False
            return all(
                (x is None and y is None) or
                (x is not None and y is not None and np.array_equal(x, y))
                for x, y in zip(current, previous)
            )

        wave_unchanged = (
            hasattr(self, "_foldspec_wave_snapshot") and
            len(self.wave) == len(self._foldspec_wave_snapshot) and
            all(np.array_equal(x, y) for x, y in zip(self.wave, self._foldspec_wave_snapshot))
        )
        errors_unchanged = (
            hasattr(self, "_foldspec_error_snapshot") and
            same_errors(self.flux_err, self._foldspec_error_snapshot)
        )
        covariance_unchanged = (
            hasattr(self, "_foldspec_covariance_snapshot") and
            len(self.flux_covariance) == len(self._foldspec_covariance_snapshot) and
            all(x is y for x, y in zip(self.flux_covariance, self._foldspec_covariance_snapshot))
        )
        use_native = (
            flux_err is None and flux_covariance is None and wave_unchanged and
            errors_unchanged and covariance_unchanged and
            np.array_equal(input_flux, self._foldspec_flux_snapshot)
        )
        if use_native:
            sources = self._native_spectra
            source_covariances = [None] * len(sources)
        else:
            errors = self.flux_err if flux_err is None else np.asarray(flux_err, dtype=float)
            if errors is not None:
                if not hasattr(errors, "__len__") or len(errors) != len(input_flux) or any(
                    e is not None and np.asarray(e).shape != wave0.shape for e in errors
                ):
                    raise ValueError("flux_err must have shape (nspec, nwave) matching the loaded spectra.")
            sources = [(wave0, f, None if errors is None else errors[i])
                       for i, f in enumerate(input_flux)]
            if flux_covariance is not None:
                if len(flux_covariance) != len(sources):
                    raise ValueError("flux_covariance must contain one matrix per spectrum.")
                source_covariances = [
                    c if isinstance(c, (SpectralCovariance, DenseSpectralCovariance))
                    else DenseSpectralCovariance(c) for c in flux_covariance
                ]
                if any(c.shape != (len(wave0), len(wave0)) for c in source_covariances):
                    raise ValueError("Each flux covariance must match the loaded wavelength grid.")
            elif flux_err is None and errors_unchanged and hasattr(self, "flux_covariance"):
                source_covariances = self.flux_covariance
                if any(c is not None and c.shape != (len(wave0), len(wave0))
                       for c in source_covariances):
                    raise ValueError("Changed wavelength sampling requires matching flux covariance.")
            else:
                source_covariances = [None] * len(sources)
                if flux_err is not None and hasattr(self, "_native_spectra"):
                    log(logging.WARNING,
                        "Explicit common-grid flux_err assumes independent pixels; pass "
                        "flux_covariance to retain correlations from previous resampling.")

        output_flux, output_covariance, source_records = [], [], []
        for index, ((wave, flu, errors), covariance) in enumerate(zip(sources, source_covariances)):
            wave, flu = np.asarray(wave, dtype=float), np.asarray(flu, dtype=float)
            noise_source = ("explicit_covariance" if flux_covariance is not None
                            else "propagated_input_covariance") if covariance is not None else "unavailable"
            if covariance is None:
                if errors is not None:
                    noise_source = "provided_standard_deviations"
                    errors = np.asarray(errors, dtype=float).copy()
                    bad = (~np.isfinite(errors)) | (errors <= 0)
                    if np.any(bad):
                        if not estimate_errors:
                            raise ValueError("Flux uncertainties must be finite and positive when estimate_errors=False.")
                        sigma = estimate_continuum_sigma(wave, flu, xor, poly_degree)
                        errors[bad] = sigma
                        noise_source = "provided_errors_with_continuum_fill"
                        log(logging.INFO, f"Estimated missing uncertainties for spectrum {index} from continuum residuals.")
                elif self.iw in (1, 2):
                    if not estimate_errors:
                        raise RuntimeError("Weighted inversion requires flux errors or covariance; estimate_errors=False.")
                    sigma = estimate_continuum_sigma(wave, flu, xor, poly_degree)
                    errors = np.full(len(wave), sigma)
                    noise_source = "continuum_residual_homoscedastic_estimate"
                    log(logging.INFO, f"Estimated spectrum {index} uncertainties from continuum residuals.")
                if errors is not None:
                    covariance = SpectralCovariance(np.square(errors))
            source_records.append({
                "spectrum_index": index,
                "samples": len(wave),
                "wavelength_sha256": array_digest(wave),
                "flux_sha256": array_digest(flu),
                "noise_source": noise_source,
                "covariance_representation_sha256": covariance_digest(covariance),
            })
            values, propagated = subtract_continuum_and_resample(
                wave, flu, observed_target, xor, poly_degree, covariance,
            )
            output_flux.append(values)
            output_covariance.append(propagated)
        self.normalised_flux = np.asarray(output_flux)
        self.normalised_flux_covariance = output_covariance
        self.normalised_flux_err = (
            np.vstack([np.sqrt(c.diagonal()) if c is not None else np.full(len(self.vell), np.nan)
                       for c in output_covariance])
            if any(c is not None for c in output_covariance) else None
        )
        self.preprocessing_used_native = use_native

        log(
            logging.INFO,
            f">> Min/Max velocities in map: {self.vell.min():.2f} / {self.vell.max():.2f} km/s",
        )


        self.pha = phases % 1.0
        self.trsp = self.normalised_flux
        if use_list_dpha:
            if self.input_dpha is None:
                raise ValueError("use_list_dpha=True requires exposure widths in the phase file.")
            dpha = np.asarray(self.input_dpha, dtype=float)
            if (dpha.shape != (self.pha.size,) or not np.all(np.isfinite(dpha))
                    or np.any(dpha <= 0) or np.any(dpha > 1)):
                raise ValueError("input_dpha must contain one finite width in (0, 1] per spectrum.")
        else:
            dpha = np.full(self.pha.size, float(self.delta_phase))
        self.used_dpha = dpha.copy()
        self.workdir.mkdir(parents=True, exist_ok=True)
        dopin_path = self.workdir / "dopin"
        with dopin_path.open("w", encoding="utf-8") as handle:
            handle.write(f"{self.pha.size:d} {self.vell.size:d} {self.lam0:.16e}\n")
            list_path = Path(self.base_dir) / self.list
            handle.write(f"{self.gama * 1e5:.16e} 0 0 {list_path}\n")

            def write_values(values):
                flat = np.asarray(values).ravel()
                for start in range(0, len(flat), 6):
                    handle.write(" ".join(f"{value:.16e}" for value in flat[start:start + 6]) + "\n")

            write_values(self.pha)
            handle.write("1\n")
            write_values(dpha)
            write_values(self.vell * 1e5)
            write_values(self.trsp.T)
            if self.iw in (1, 2):
                if (self.normalised_flux_err is None or
                        np.any(~np.isfinite(self.normalised_flux_err)) or
                        np.any(self.normalised_flux_err <= 0)):
                    raise ValueError("Weighted inversion requires positive propagated uncertainties.")
                write_values(self.normalised_flux_err.T)

        # This is bounded preprocessing provenance, not a claim to record
        # upstream calibration, ephemeris construction, or observing metadata.
        import hashlib
        import json
        provenance = {
            "schema_version": 1,
            "object": str(self.object),
            "dopin_sha256": hashlib.sha256(dopin_path.read_bytes()).hexdigest(),
            "continuum_frame": "observed_wavelength",
            "continuum_bands_angstrom": [float(value) for value in xor],
            "continuum_polynomial_degree": int(poly_degree),
            "rest_wavelength_angstrom": float(self.lam0),
            "systemic_velocity_km_s": float(self.gama),
            "rest_frame_window_half_width_angstrom": float(self.delw),
            "velocity_definition": "longitudinal_relativistic",
            "interpolation": "linear_flux_density_samples",
            "source_grid": "native" if use_native else "loaded_common",
            "estimate_errors_enabled": bool(estimate_errors),
            "error_estimator": "continuum_residual_standard_deviation_with_residual_dof",
            "covariance_propagation": "composed_sampling_and_weighted_continuum_fit",
            "phase_sha256": array_digest(self.pha),
            "exposure_width_sha256": array_digest(dpha),
            "sources": source_records,
        }
        (self.workdir / "dopin.json").write_text(
            json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8",
        )

        if not plot or plt is None:
            return

        delp = float(self.delta_phase)
        if use_list_dpha and self.input_dpha is not None and np.any(np.isfinite(self.input_dpha)):
            delp = np.asarray(self.input_dpha, dtype=float)

        if rebin:
            trail, phase = rebin_trail(
                self.vell,
                self.normalised_flux,
                self.input_phase,
                self.nbins,
                delp,
                rebin_wave=None,
            )
            phase_pad = 1.0 / self.nbins
        else:
            phase_once, inverse = np.unique(self.pha, return_inverse=True)
            trail_once = np.zeros((len(self.vell), len(phase_once)))
            for index in range(len(phase_once)):
                trail_once[:, index] = np.mean(self.normalised_flux[inverse == index], axis=0)
            if two_orbits:
                phase = np.concatenate([phase_once, phase_once + 1.0])
                trail = np.concatenate([trail_once, trail_once], axis=1)
            else:
                phase = phase_once
                trail = trail_once
            diffs = np.diff(phase)
            diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
            phase_pad = float(np.nanmedian(diffs)) if diffs.size else 1.0 / self.nbins


        if xlim is None:
            rr = np.ones(self.normalised_wave.size,dtype='bool')
        else:
            if len(xlim) != 2 or not np.all(np.isfinite(xlim)) or xlim[0] >= xlim[1]:
                raise ValueError("xlim must contain two finite, increasing rest-frame wavelengths.")
            rr = (self.normalised_wave > xlim[0]) & (self.normalised_wave < xlim[1])
            if not np.any(rr):
                raise ValueError("xlim contains no output wavelengths.")

        if vel_space:
            x1_lim = float(np.nanmin(self.vell))
            x2_lim = float(np.nanmax(self.vell))
            profile_x = self.vell
            display_xlim = (x1_lim, x2_lim)
            xlabel = 'Velocity / km s$^{-1}$'
            center = 0.0
        else:
            x1_lim = float(np.nanmin(self.normalised_wave))
            x2_lim = float(np.nanmax(self.normalised_wave))
            profile_x = self.normalised_wave
            display_xlim = (self.lam0 - self.delw, self.lam0 + self.delw)
            xlabel = 'Wavelength / $\\AA$'
            center = self.lam0
        if xlim is not None:
            display_xlim = (float(profile_x[rr].min()), float(profile_x[rr].max()))

        fig = plt.figure('Trail',figsize=(6.57,8.57))
        plt.clf()
        ax1 = plt.subplot2grid((6, 1), (0, 0), rowspan=2)
        ax1.minorticks_on()
        profile = (
            np.nanmedian(self.normalised_flux, axis=0)
            if plot_median
            else np.nanmean(self.normalised_flux, axis=0)
        )
        profile_mask = rr if rr.shape == profile.shape else np.ones(profile.shape, dtype=bool)
        ax1.plot(profile_x[profile_mask], profile[profile_mask], color='k')
        ax1.axhline(y=0, ls='--', color='k', alpha=0.7)
        ax1.axvline(x=center, ls='--', color='DarkOrange')
        ax1.set_xlim(*display_xlim)
        ax1.set_ylabel('Continuum-subtracted flux')
        ax1.set_xticklabels([])

        ax2 = plt.subplot2grid((6, 1), (2, 0), rowspan=4)
        ax2.minorticks_on()
        trail_cmap = plt.cm.binary.copy() if hasattr(plt.cm.binary, "copy") else plt.cm.binary
        if hasattr(trail_cmap, "set_bad"):
            trail_cmap.set_bad(color="0.85")
        if len(phase) == 1:
            phase_edges = np.array([phase[0] - phase_pad / 2., phase[0] + phase_pad / 2.])
        else:
            phase_edges = _bin_edges(phase)
        ax2.pcolormesh(_bin_edges(profile_x), phase_edges, trail.T,
                       cmap=trail_cmap, shading='flat')
        ax2.set_xlim(*display_xlim)
        ax2.set_xlabel(xlabel)
        ax2.axvline(x=center, ls='--', color='DarkOrange')

        if rebin:
            if two_orbits:
                lim_two = 2
            else:
                lim_two = 1
            ax2.set_ylim(0., lim_two)
        else:
            ax2.set_ylim(phase_edges[0], phase_edges[-1])
        ax2.set_ylabel('Orbital Phase')
        fig.tight_layout(h_pad=0)
        if show:
            plt.show()



    def Syncdop(self, nri=None, ndi=None):
        """Build and execute the solver, returning validated fit diagnostics.

        Every call executes the solver; the executable is compiled once per
        source and compiler version and serves every data and map size.
        A run succeeds only if the optimizer converges, the RMS target is met,
        and the finite output matches the current input. ``run.json`` records
        provenance. With supplied errors, ``covariance_mode='full'`` (default)
        uses the covariance propagated by :meth:`Dopin`; ``'diagonal'`` is an
        explicit independence approximation that overstates the information in
        resampled data. The fit statistic is RMS, not reduced chi².
        ``nri``/``ndi`` are obsolete and ignored.
        """
        from .runtime import run

        if nri is not None or ndi is not None:
            import warnings
            warnings.warn("Syncdop(nri, ndi) is obsolete: projection storage is now sized "
                          "exactly at run time.", DeprecationWarning, stacklevel=2)
        return run(self)

    def Dopmap(
        self, dopout: str = "dop.out", cmaps=None, limits=None,
        colorbar: bool = False, negative: bool = False,
        remove_mean: bool = False, corrx=0, corry=0,
        smooth: Union[bool, float] = False,
        plot: Optional[bool] = None, show: Optional[bool] = None,
    ):
        """Read a map and return ``(colorbar, displayed_array)``.

        The returned array is identical to the displayed array, including
        mean removal, smoothing and sign reversal. Pixels outside the
        reconstruction's circular support are NaN. Without mean removal the
        map is divided by its positive peak, so zero flux stays at zero and
        negative values stay negative; default limits are 0 and the 99.5th
        percentile. Mean-subtracted maps use absolute-maximum normalization
        and symmetric default limits. ``smooth=True`` applies a Gaussian of sigma one pixel;
        a positive float specifies another sigma. Smoothing precedes radial
        mean removal and normalization. ``negative`` reverses the final sign;
        explicit ``limits`` describe the values before that sign reversal.

        Axes are projected velocities in km/s with the convention
        ``v_rad - gamma = -Vx*cos(2*pi*phase) + Vy*sin(2*pi*phase)``.
        ``corrx/corry`` shift the mean-removal center in pixel units, retaining
        the historical convention of subtracting each correction. Raw map
        values and physical grids are available through ``_parse_dopout``.
        """
        if plot is None:
            plot = bool(self.plot)
        if show is None:
            show = bool(plot)
        if not plot and colorbar:
            raise ValueError("colorbar=True requires plot=True.")
        dopout_path = Path(dopout)
        if not dopout_path.is_absolute():
            dopout_path = self.workdir / dopout_path
        parsed = _read_verified_dopout(dopout_path)
        data = np.asarray(parsed.im, dtype=float).copy()
        if not np.all(np.isfinite(data)):
            raise RuntimeError(f"Non-finite map values found in {dopout_path}.")
        vx, vy = np.meshgrid(parsed.vmap, parsed.vmap)
        map_step = parsed.vmap[1] - parsed.vmap[0]
        support = np.hypot(vx, vy) < parsed.va + map_step / 10
        data[~support] = np.nan
        if not np.any(support):
            raise RuntimeError(f"The map grid has no pixels inside its circular support: {dopout_path}.")
        sigma = float(smooth)
        if not np.isfinite(sigma) or sigma < 0:
            raise ValueError("smooth must be False, True, or a non-negative sigma.")
        if sigma:
            # Separable Gaussian with reflected edges, without a SciPy dependency.
            radius = max(1, int(np.ceil(4 * sigma)))
            kernel = np.exp(-0.5 * (np.arange(-radius, radius + 1) / sigma) ** 2)
            kernel /= kernel.sum()
            numerator = np.where(support, data, 0.0)
            weights = support.astype(float)
            for axis in (0, 1):
                pads = [(0, 0), (0, 0)]
                pads[axis] = (radius, radius)
                def convolve(values):
                    return np.apply_along_axis(
                        lambda row: np.convolve(row, kernel, mode="valid"),
                        axis, np.pad(values, pads, mode="reflect"),
                    )
                numerator, weights = convolve(numerator), convolve(weights)
            data = np.full_like(numerator, np.nan)
            np.divide(numerator, weights, out=data, where=support & (weights > 0))
        if remove_mean:
            center = [(data.shape[1] - 1) / 2 - corrx,
                      (data.shape[0] - 1) / 2 - corry]
            data -= create_profile(data, radial_profile(data, center), center)
            display = _scale_by_absmax(data)
        else:
            display = _scale_by_peak(data)
        if limits is None:
            finite = display[np.isfinite(display)]
            if remove_mean:
                bound = float(np.percentile(np.abs(finite), 99.5)) if finite.size else 1.0
                limits = [-bound, bound] if bound > 0 else [-1.0, 1.0]
            else:
                top = float(np.percentile(finite, 99.5)) if finite.size else 1.0
                limits = [0.0, top if top > 0 else 1.0]
        if negative:
            display = -display
            limits = [-limits[1], -limits[0]]
        if not plot:
            return None, display

        _, plt = _lazy_import_matplotlib_pyplot()
        if cmaps is None:
            cmaps = plt.cm.magma
        if hasattr(cmaps, "with_extremes"):
            cmap_plot = cmaps.with_extremes(bad="0.85")
        else:
            cmap_plot = cmaps.copy() if hasattr(cmaps, "copy") else cmaps
            if hasattr(cmap_plot, "set_bad"):
                cmap_plot.set_bad(color="0.85")
        fig = plt.figure(num="Doppler Map", figsize=(8.57, 8.57))
        fig.clear()
        ax = fig.add_subplot(111)
        ax.minorticks_on()
        edges = _bin_edges(parsed.vmap)
        # Fortran writes x-major. Keeping its array order here exchanges the
        # internal axes, producing the public radial-velocity convention above.
        img = ax.imshow(
            display, interpolation="nearest", cmap=cmap_plot, aspect="equal",
            origin="lower", extent=(edges[0], edges[-1], edges[0], edges[-1]),
            vmin=limits[0], vmax=limits[1],
        )
        ax.set_xlabel('V$_x$ / km s$^{-1}$')
        ax.set_ylabel('V$_y$ / km s$^{-1}$')
        cbar = None
        if colorbar:
            cbar_obj = fig.colorbar(img, ax=ax, format="%.2f", fraction=0.046, pad=0.04)
            cbar_obj.set_label("Flux / map peak" if not remove_mean else "Flux / max |flux|")
            cbar = DraggableColorbar(cbar_obj, img)
            cbar.connect()
        fig.tight_layout()
        if show:
            plt.show()
        return cbar, display

    def Reco(
        self, cmaps=None, limits=None, colorbar: bool = True,
        plot: Optional[bool] = None, show: Optional[bool] = None,
    ):
        """Return ``(data_colorbar, model_colorbar, model, data)``.

        Returned trails contain the input and reconstructed values restored
        to the original flux units and velocity/phase sampling, shape
        ``(nvp, nph)``. Only the plots are phase-binned. Both panels use
        the same normalization (the absolute maximum of the binned data), so
        reconstruction amplitude errors remain visible. Phases increase
        upwards and the velocity coordinates are bin centers in km/s.
        """
        if plot is None:
            plot = bool(self.plot)
        if show is None:
            show = bool(plot)
        if not plot and colorbar:
            raise ValueError("colorbar=True requires plot=True.")
        parsed = _read_verified_dopout(self.workdir / "dop.out")
        dm, dmr = np.asarray(parsed.dm), np.asarray(parsed.dmr)
        if not plot:
            return None, None, dmr, dm
        _, plt = _lazy_import_matplotlib_pyplot()
        if cmaps is None:
            cmaps = plt.cm.binary
        trail_dm, phase = rebin_trail(
            parsed.vp, dm.T, parsed.pha, self.nbins, parsed.dpha)
        trail_dmr, _ = rebin_trail(
            parsed.vp, dmr.T, parsed.pha, self.nbins, parsed.dpha)
        finite_data = np.abs(trail_dm[np.isfinite(trail_dm)])
        scale = float(finite_data.max()) if finite_data.size else 1.0
        if scale == 0:
            scale = 1.0
        trails = [trail_dm / scale, trail_dmr / scale]
        if limits is None:
            limits = _auto_display_limits(np.concatenate([x.ravel() for x in trails]))
        if hasattr(cmaps, "with_extremes"):
            cmap_plot = cmaps.with_extremes(bad="0.85")
        else:
            cmap_plot = cmaps.copy() if hasattr(cmaps, "copy") else cmaps
            if hasattr(cmap_plot, "set_bad"):
                cmap_plot.set_bad(color="0.85")
        xedges, yedges = _bin_edges(parsed.vp), _bin_edges(phase)
        fig = plt.figure("Reconstruction", figsize=(10, 8))
        fig.clear()
        cbars = []
        for index, (trail, title) in enumerate(zip(trails, ("Input", "Reconstruction"))):
            ax = fig.add_subplot(1, 2, index + 1)
            # pcolormesh retains exact bin edges even when output formatting
            # has rounded spectral center coordinates slightly unevenly.
            img = ax.pcolormesh(xedges, yedges, trail.T, shading="flat",
                                cmap=cmap_plot, vmin=limits[0], vmax=limits[1])
            ax.set_xlabel('Velocity / km s$^{-1}$')
            ax.set_title(title)
            if index == 0:
                ax.set_ylabel("Orbital Phase")
            else:
                ax.tick_params(labelleft=False)
            cbar = None
            if colorbar:
                from matplotlib.ticker import FormatStrFormatter
                from mpl_toolkits.axes_grid1 import make_axes_locatable
                cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.08)
                cbar_obj = fig.colorbar(img, cax=cax, format=FormatStrFormatter("%.2f"))
                cbar_obj.set_label("Flux / data peak")
                cbar = DraggableColorbar(cbar_obj, img)
                cbar.connect()
            cbars.append(cbar)
        fig.tight_layout(w_pad=0)
        if show:
            plt.show()
        return cbars[0], cbars[1], dmr, dm

    def Residuals(
        self,
        dm: Optional[np.ndarray] = None,
        dmr: Optional[np.ndarray] = None,
        sigma: Optional[np.ndarray] = None,
        bins: int = 50,
        vmax: float = 6.0,
        plot: Optional[bool] = None,
        show: Optional[bool] = None,
        cmaps=None,
        *,
        whiten: bool = False,
        mask_excluded: bool = True,
    ):
        """Return residuals and uncertainty-scaled residuals on the original phases.

        Supplied ``sigma`` must be positive 1-sigma errors in the data's units.
        Otherwise errors from the matching validated run are used. If no errors
        are available, the second return value is ``None``: flux densities are
        never interpreted as Poisson counts. With ``whiten=True``, require a
        validated full-covariance fit and use its retained covariance snapshot,
        discarding numerically null eigenmodes;
        marginally standardized residuals otherwise remain correlated.
        Excluded phases are NaN by default. Plot arguments are deprecated and
        ignored; use :meth:`Reco` for the input/model comparison and plot the
        returned residual arrays separately.
        """
        from .runtime import _sha, covariance_modes

        if (dm is None) != (dmr is None):
            raise ValueError("Pass both dm and dmr, or neither.")
        if dm is None:
            _, _, dmr, dm = self.Reco(plot=False, colorbar=False)
        dm = np.asarray(dm, dtype=float)
        dmr = np.asarray(dmr, dtype=float)
        if dm.ndim != 2 or dm.shape != dmr.shape:
            raise ValueError("dm and dmr must be matching (nvp, nph) arrays.")
        if not np.isfinite(dm).all() or not np.isfinite(dmr).all():
            raise ValueError("Residual inputs must be finite.")
        residuals = dm - dmr
        matched = False
        output_path = self.workdir / "dop.out"
        if getattr(self, "_output_hash", None) and output_path.exists():
            if self._output_hash == _sha(output_path):
                parsed = _parse_dopout(output_path)
                matched = (parsed.dm.shape == dm.shape and
                           np.array_equal(parsed.dm, dm) and
                           np.array_equal(parsed.dmr, dmr))
        if whiten:
            if sigma is not None:
                raise ValueError("whiten=True uses full covariance; do not supply sigma.")
            covariance = getattr(self, "_fit_covariance", None)
            if not matched or covariance is None:
                raise ValueError("Whitening requires a matching validated full-covariance fit.")
            scaled = np.empty_like(residuals)
            for index, cov in enumerate(covariance):
                values, vectors = covariance_modes(cov, self._fit_rcond)
                scaled[:, index] = vectors @ ((vectors.T @ residuals[:, index]) / np.sqrt(values))
        else:
            if sigma is None and matched:
                sigma = getattr(self, "_fit_sigma", None)
            if sigma is None:
                scaled = None
            else:
                sigma = np.asarray(sigma, dtype=float)
                if sigma.shape != dm.shape or not np.isfinite(sigma).all() or np.any(sigma <= 0):
                    raise ValueError("sigma must be positive, finite and have the data shape.")
                scaled = residuals / sigma
        if mask_excluded and matched:
            valid = np.asarray(self._fit_phase_mask, dtype=bool)
            residuals[:, ~valid] = np.nan
            if scaled is not None:
                scaled[:, ~valid] = np.nan
        return residuals, scaled


def rebin_trail(waver, flux, input_phase, nbins, delp, rebin_wave=None):
    """Exposure-overlap-weighted diagnostic trail on periodic phase bins.

    ``flux`` is (nspec, nwave); ``input_phase`` gives exposure centers in
    cycles and ``delp`` their positive widths (scalar or per spectrum).
    The return values are (nwave, 2*nbins+2) and the corresponding bin
    centers. The first and last columns pad two full orbits. Empty bins are
    NaN. Each finite spectral pixel is weighted by its exposure's exact
    overlap with a bin; missing pixels do not contaminate other exposures.
    This is a diagnostic average, not an inverse-variance statistical binning.
    ``rebin_wave`` must be None or zero; wavelength rebinning is not supported.
    """
    if rebin_wave is not None and rebin_wave != 0:
        raise ValueError("rebin_wave is unsupported; use the native spectral sampling.")
    waver = np.asarray(waver, dtype=float)
    flux = np.asarray(flux, dtype=float)
    input_phase = np.asarray(input_phase, dtype=float)
    if waver.ndim != 1 or not waver.size or not np.all(np.isfinite(waver)):
        raise ValueError("waver must be a non-empty finite 1D array.")
    if input_phase.ndim != 1 or not np.all(np.isfinite(input_phase)):
        raise ValueError("input_phase must be a finite 1D array.")
    if flux.ndim != 2 or flux.shape != (input_phase.size, waver.size):
        raise ValueError("flux must have shape (nspec, nwave).")
    if isinstance(nbins, (bool, np.bool_)) or int(nbins) != nbins or nbins <= 0:
        raise ValueError("nbins must be a positive integer.")
    nbins = int(nbins)
    widths = np.asarray(delp, dtype=float)
    if widths.ndim == 0:
        widths = np.full(input_phase.shape, float(widths))
    if widths.shape != input_phase.shape or not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("delp must contain finite positive exposure widths, one per spectrum.")
    lo = np.arange(nbins, dtype=float) / nbins
    hi = lo + 1.0 / nbins
    weights = np.empty((input_phase.size, nbins), dtype=float)
    for i, (phi, width) in enumerate(zip(input_phase, widths)):
        # Each complete orbit contributes one full bin. Treat the remaining
        # interval in [0,2), including an interval crossing phase zero.
        full_orbits = np.floor(width)
        remainder = width - full_orbits
        start = (phi - width / 2.0) % 1.0
        stop = start + remainder
        overlap = np.maximum(0, np.minimum(hi, stop) - np.maximum(lo, start))
        overlap += np.maximum(0, np.minimum(hi + 1, stop) - np.maximum(lo + 1, start))
        # Roundoff at exactly touching edges must not create a full-valued
        # ghost bin after dividing a tiny numerator by the same tiny weight.
        overlap[overlap <= 8 * np.finfo(float).eps] = 0.0
        weights[i] = full_orbits / nbins + overlap
    valid = np.isfinite(flux)
    numerator = np.where(valid, flux, 0.0).T @ weights
    denominator = valid.astype(float).T @ weights
    trail = np.full_like(numerator, np.nan)
    np.divide(numerator, denominator, out=trail, where=denominator > 0)
    indices = np.arange(-1, 2 * nbins + 1)
    return trail[:, indices % nbins], (indices + 0.5) / nbins


class DraggableColorbar(object):
    def __init__(self, cbar, mappable):
        _, plt = _lazy_import_matplotlib_pyplot()
        self.cbar = cbar
        self.mappable = mappable
        self.press = None
        self.cycle = sorted([i for i in dir(plt.cm) if hasattr(getattr(plt.cm,i),'N')])
        current_cmap = None
        if hasattr(cbar, "get_cmap"):
            current_cmap = cbar.get_cmap()
        elif hasattr(cbar, "cmap"):
            current_cmap = cbar.cmap
        elif hasattr(cbar, "mappable"):
            current_cmap = cbar.mappable.get_cmap()
        else:
            current_cmap = mappable.get_cmap()
        self.index = self.cycle.index(current_cmap.name) if current_cmap else 0

    def _refresh(self) -> None:
        if hasattr(self.cbar, "update_normal"):
            self.cbar.update_normal(self.mappable)
        elif hasattr(self.cbar, "draw_all"):
            self.cbar.draw_all()
        if hasattr(self.cbar, "ax") and self.cbar.ax and self.cbar.ax.figure:
            self.cbar.ax.figure.canvas.draw()

    def connect(self):
        """connect to all the events we need"""
        canvas = self.cbar.ax.figure.canvas
        self.cidpress = canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.keypress = canvas.mpl_connect('key_press_event', self.key_press)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.cbar.ax:
            return
        self.press = event.x, event.y

    def key_press(self, event):
        if event.key=='down':
            self.index += 1
        elif event.key=='up':
            self.index -= 1
        if self.index<0:
            self.index = len(self.cycle) - 1
        elif self.index>=len(self.cycle):
            self.index = 0
        cmap = self.cycle[self.index]
        if hasattr(self.cbar, "set_cmap"):
            self.cbar.set_cmap(cmap)
        self.mappable.set_cmap(cmap)
        #self.mappable.get_axes().set_title(cmap)
        self._refresh()

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None:
            return
        if event.inaxes != self.cbar.ax:
            return
        xprev, yprev = self.press
        dy = event.y - yprev
        self.press = event.x,event.y
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        scale = self.cbar.norm.vmax - self.cbar.norm.vmin
        perc = 0.03
        if event.button==1:
            self.cbar.norm.vmin -= (perc*scale)*np.sign(dy)
            self.cbar.norm.vmax -= (perc*scale)*np.sign(dy)
        elif event.button==3:
            self.cbar.norm.vmin -= (perc*scale)*np.sign(dy)
            self.cbar.norm.vmax += (perc*scale)*np.sign(dy)
        _set_mappable_norm(self.mappable, self.cbar.norm)
        self._refresh()


    def on_release(self, event):
        """on release we reset the press data"""
        self.press = None
        _set_mappable_norm(self.mappable, self.cbar.norm)
        self._refresh()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        canvas = self.cbar.ax.figure.canvas
        canvas.mpl_disconnect(self.cidpress)
        canvas.mpl_disconnect(self.cidrelease)
        canvas.mpl_disconnect(self.cidmotion)


def radial_profile(data, center):
    """Mean of finite pixels in integer-radius annuli about a pixel center."""
    data = np.asarray(data, dtype=float)
    y, x = np.indices(data.shape)
    r = np.hypot(x - center[0], y - center[1]).astype(int)
    valid = np.isfinite(data)
    sums = np.bincount(r[valid], weights=data[valid], minlength=r.max() + 1)
    counts = np.bincount(r[valid], minlength=r.max() + 1)
    result = np.full(sums.shape, np.nan)
    np.divide(sums, counts, out=result, where=counts > 0)
    return result


def create_profile(data, profile, center):
    """Expand an annular profile to every pixel, including the outer annulus."""
    y, x = np.indices(np.shape(data))
    r = np.hypot(x - center[0], y - center[1]).astype(int)
    return np.asarray(profile)[r]


def stream(q, k1, porb, m1, inc, colors='k', both_lobes=False, title=True, label=None):
    """Overlay projected ballistic/Keplerian velocities and Roche outlines.

    ``q=M2/M1``, ``k1`` is the observed primary semi-amplitude in km/s,
    ``porb`` is in days, ``m1`` in solar masses and ``inc`` in degrees.
    K1 sets the velocity scale; the redundant mass/period/inclination values
    are checked and a warning reports discrepancies greater than 2 percent.
    Pass ``k1=None`` to derive the scale from those three quantities instead.
    """
    import warnings
    q = _positive_scalar(q, "q")
    porb = _positive_scalar(porb, "porb")
    m1 = _positive_scalar(m1, "m1")
    inc = float(inc)
    if not np.isfinite(inc) or not 0 <= inc <= 90:
        raise ValueError("inc must be finite and between 0 and 90 degrees.")
    omega = 2 * np.pi / (porb * 86400)
    scale_from_mass = (6.67430e-8 * 1.98847e33 * m1 * (1 + q) * omega) ** (1 / 3)
    scale_from_mass *= np.sin(np.deg2rad(inc)) / 1e5
    expected_k1 = q / (1 + q) * scale_from_mass
    if k1 is None:
        scale = scale_from_mass
    else:
        k1 = _positive_scalar(k1, "k1")
        if inc == 0:
            raise ValueError("A face-on binary cannot have a positive projected K1.")
        if not np.isclose(k1, expected_k1, rtol=0.02, atol=0):
            warnings.warn(
                f"K1={k1:g} km/s differs from the M1/P/q/i prediction "
                f"({expected_k1:.6g} km/s); the overlay uses the supplied K1.",
                RuntimeWarning, stacklevel=2,
            )
        scale = k1 * (1 + q) / q
    xl, yl, _, _, ballistic, keplerian = stream_calculate(q)
    cm = q / (1 + q)
    _, plt = _lazy_import_matplotlib_pyplot()
    plt.figure(num="Doppler Map")
    plt.plot(ballistic.real * scale, ballistic.imag * scale, color=colors)
    plt.plot(keplerian.real * scale, keplerian.imag * scale, color=colors)
    quarter = xl.size // 4
    donor = slice(quarter, 3 * quarter)
    plt.plot(yl[donor] * scale, (xl[donor] - cm) * scale, color=colors)
    if both_lobes:
        ix = np.r_[np.arange(3 * quarter, xl.size), np.arange(quarter)]
        plt.plot(yl[ix] * scale, (xl[ix] - cm) * scale, color=colors, ls="--")
    plt.plot(0, 0, "x", ms=9, c=colors, alpha=0.3)
    plt.plot(0, -cm * scale, "+", ms=10, c=colors, alpha=0.7)
    plt.plot(0, (1 - cm) * scale, "+", ms=10, c=colors, alpha=0.7)
    if title:
        plt.title(rf"$i$={inc:g}, $M_1$={m1:g} $M_\odot$, $q$={q:g}, $P_{{orb}}$={porb * 24:.3g} hr")
    if label is not None:
        plt.text(0.12, 0.1, label, ha="center", va="center", transform=plt.gca().transAxes)
    plt.tight_layout()


def _positive_scalar(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be > 0 and finite.")
    return value


def stream_calculate(qm, ni=100, nj=100, *, rtol=1e-8, atol=1e-10, max_angle=70.0):
    """Roche outlines and a cold ballistic stream launched just inside L1.

    Lengths use separation a=1; velocities use Omega*a=1. Return values are
    ``(lobe_x, lobe_y, stream_x, stream_y, ballistic_v, keplerian_v)``.
    Stream positions retain the historical primary-centered convention with
    ``stream_y=-physical_y``; complex velocities use physical (+x,+y), in the
    inertial frame expressed on the instantaneous binary axes. The stream
    starts at L1-0.001a at rest in the rotating frame and ends at max_angle
    degrees about the primary. This is a test-particle trajectory, not a gas
    hydrodynamics model. Adaptive RK4 step doubling controls local error.
    Returned trajectory arrays contain only accepted samples (no zero pad).
    ``nj`` remains accepted for compatibility; only the orbital-plane outline
    is needed, so a full latitude/longitude surface is not computed.
    """
    qm = _positive_scalar(qm, "Mass ratio q")
    rtol = _positive_scalar(rtol, "rtol")
    atol = _positive_scalar(atol, "atol")
    if int(ni) != ni or ni < 3 or int(nj) != nj or nj < 1:
        raise ValueError("ni must be an integer >= 3 and nj an integer >= 1.")
    max_angle = float(max_angle)
    if not np.isfinite(max_angle) or not 0 < max_angle < 180:
        raise ValueError("max_angle must lie between 0 and 180 degrees.")
    rl1 = rlq1(qm)
    x, y = lobes(qm, rl1, int(ni), 1)
    cm = qm / (1 + qm)
    z1, z2 = -cm, 1 - cm
    z, w = complex(rl1 - cm - 1e-3), 0j
    dt, elapsed = 1e-4, 0.0
    target = np.deg2rad(max_angle)
    positions, velocities = [z], [w]

    def step(z, w, dt):
        dz, dw = intrk(z, w, dt, z1, z2, qm)
        return z + dz, w + dw

    for _ in range(100000):
        coarse_z, coarse_w = step(z, w, dt)
        half_z, half_w = step(z, w, dt / 2)
        fine_z, fine_w = step(half_z, half_w, dt / 2)
        err_z, err_w = (fine_z - coarse_z) / 15, (fine_w - coarse_w) / 15
        error = max(abs(err_z) / (atol + rtol * max(abs(z), abs(fine_z))),
                    abs(err_w) / (atol + rtol * max(abs(w), abs(fine_w))))
        if not np.isfinite(error):
            raise RuntimeError("Non-finite ballistic stream integration step.")
        if error <= 1:
            next_z, next_w = fine_z + err_z, fine_w + err_w
            if np.angle(next_z - z1) >= target:
                # Locate the final angular event, so changing tolerances does
                # not change the requested endpoint of the plotted stream.
                lo, hi = 0.0, dt
                for _ in range(40):
                    mid = (lo + hi) / 2
                    zm, wm = step(z, w, mid / 2)
                    zm, wm = step(zm, wm, mid / 2)
                    if np.angle(zm - z1) < target:
                        lo = mid
                    else:
                        hi = mid
                z, w = zm, wm
                positions.append(z)
                velocities.append(w)
                break
            z, w = next_z, next_w
            positions.append(z)
            velocities.append(w)
            elapsed += dt
            if elapsed > 100:
                raise RuntimeError("Ballistic stream did not reach the requested angle.")
        dt *= min(2.0, max(0.2, 0.9 * error ** (-0.2))) if error > 0 else 2.0
        dt = min(dt, 0.05)
        if dt < 1e-12:
            raise RuntimeError("Ballistic stream integration requires a singular step.")
    else:
        raise RuntimeError("Ballistic stream exceeded its integration step limit.")
    z = np.asarray(positions)
    w = np.asarray(velocities)
    r = np.abs(z - z1)
    inertial = w + 1j * z
    keplerian = 1j * (z - z1) / r / np.sqrt(r * (1 + qm)) - 1j * cm
    return x, y, z.real + cm, -z.imag, inertial, keplerian


def rlq1(q):
    """Primary-to-L1 distance in units of a (not a volume-equivalent radius)."""
    q = _positive_scalar(q, "q")
    if q == 1:
        return 0.5
    # The L1 equation is strictly increasing between the stellar centers.
    # Bracketing avoids Newton guesses outside this physical interval.
    lo, hi = 0.0, 1.0
    for _ in range(60):
        r = (lo + hi) / 2
        f = q / (1 - r) ** 2 - 1 / r ** 2 + (1 + q) * r - q
        if f < 0:
            lo = r
        else:
            hi = r
    return (lo + hi) / 2


def lobes(q,rs,ni,nj):
    '''
    SUBROUTINE
    '''
    r = np.zeros((ni,nj))
    ch = np.zeros(ni)
    ps = np.zeros(nj)
    x  = np.zeros(ni)
    y  = np.zeros(ni)
    x2 = np.zeros(ni)
    y2 = np.zeros(ni)
    nc = ni
    nop = nj

    r,ch,ps = surface(q,rs,nc,nop,r,ch,ps)
    j=0
    for i in np.arange(nc):
        x[i] = 1.0 -r[i,j]*np.cos(ch[i])
        y[i] = -r[i,j] * np.sin(ch[i])

    r,ch,ps = surface(1./q,1.-rs,nc,nop,r,ch,ps)
    j=0
    for i in np.arange(nc):
        x2[i] = r[i,j] * np.cos(ch[i])
        y2[i] = r[i,j] * np.sin(ch[i])
    xt = np.concatenate((x2[::-1],x,x[::-1],x2))
    yt = np.concatenate((y2[::-1],-y,y[::-1],-y2))
    return xt,yt


def pot(q,x,y,z):
    '''
    FUNCTION
    Roche potential. coordinates centered on M2,
    z along rotation axis, x toward M1
    pr is gradient in radius from M2
    first transform to polar coordinates w/r rotation axis
    '''
    r = np.sqrt(x*x+y*y+z*z)
    if (r == 0):
        raise ValueError("Invalid coordinates in pot(): r=0 leads to a singular potential.")
    rh = np.sqrt(x*x+y*y)
    st=rh/r
    if rh == 0:
        cf=1
    else:
        cf=x/rh

    r2 = 1. / (1. + q)
    r1 = np.sqrt(1.0+r**2-2.0*r*cf*st)
    pot=-1.0/r-1.0/q/r1-0.5*(1.0/q+1.0)*(r2**2+(r*st)**2-2.0*r2*r*cf*st)
    pr=1.0/r**2+1.0/q/(r1**3)*(r-cf*st)-0.5*(1.0/q+1)*2.0*(r*st*st-r2*cf*st)
    return pot,pr


def surface(q, rs, nc, nop, r, ch, ps):
    """Critical Roche surface around M2; rs is primary-to-L1 distance.

    Chi is measured from the direction towards M1 and psi from the orbital
    plane. Bracketed radial solves preserve the critical equipotential,
    including the equal-mass limit, without perturbing the mass ratio.
    """
    q = _positive_scalar(q, "q")
    if int(nc) != nc or nc < 3 or int(nop) != nop or nop < 1:
        raise ValueError("nc must be >= 3 and nop >= 1.")
    if not np.isfinite(rs) or not 0 < rs < 1:
        raise ValueError("rs must be between the two stellar centers.")
    ch = np.linspace(0, np.pi, int(nc))
    ps = np.arange(int(nop)) * (2 * np.pi / nop)
    direction_x = np.cos(ch)[:, None]
    direction_y = np.sin(ch)[:, None] * np.cos(ps)[None, :]
    rs1 = 1 - rs
    target, _ = pot(q, rs1, 0.0, 0.0)
    lo = np.zeros((int(nc), int(nop)))
    hi = np.full_like(lo, rs1)
    cm_from_m2 = 1 / (1 + q)
    for _ in range(60):
        radius = (lo + hi) / 2
        distance_m1 = np.sqrt(1 + radius ** 2 - 2 * radius * direction_x)
        potential = (-1 / radius - 1 / (q * distance_m1)
                     - 0.5 * (1 + 1 / q) * (cm_from_m2 ** 2
                       + radius ** 2 * (direction_x ** 2 + direction_y ** 2)
                       - 2 * cm_from_m2 * radius * direction_x))
        inside = potential < target
        lo = np.where(inside, radius, lo)
        hi = np.where(inside, hi, radius)
    r = (lo + hi) / 2
    r[0, :] = rs1
    return r, ch, ps


def eqmot(z,w,z1,z2,qm):
    zr1 = z-z1
    zr2 = z-z2
    ## c change by Guillaume : - sign in Coriolis
    wp = -(qm*zr2/(np.abs(zr2))**3 + zr1/(np.abs(zr1))**3)/(1.0+qm) - complex(0.,2.)*w + z
    zp = w
    return zp,wp


def intrk(z,w,dt,z1,z2,qm):
    zx=z
    wx=w
    zp,wp = eqmot(zx,wx,z1,z2,qm)
    hz0=zp*dt
    hw0=wp*dt
    zx=z+hz0/2.
    wx=w+hw0/2.
    zp,wp = eqmot(zx,wx,z1,z2,qm)
    hz1=zp*dt
    hw1=wp*dt
    zx=z+hz1/2.
    wx=w+hw1/2.
    zp,wp = eqmot(zx,wx,z1,z2,qm)
    hz2=zp*dt
    hw2=wp*dt
    zx=z+hz2
    wx=w+hw2
    zp,wp = eqmot(zx,wx,z1,z2,qm)
    hz3=zp*dt
    hw3=wp*dt
    dz=(hz0+2*hz1+2*hz2+hz3)/6.
    dw=(hw0+2*hw1+2*hw2+hw3)/6.
    return dz,dw

def resonance(j, k, k1, q, porb=None, m1=None):
    """Plot projected j:(j-k) resonance and tidal-truncation velocity circles.

    q=M2/M1; K1 is the observed primary semi-amplitude in km/s. The dashed
    resonance radius is r/a=[(j-k)/j]**(2/3)/(1+q)**(1/3). The solid circle
    uses the approximate tidal radius 0.60*a/(1+q), not the circularization
    radius. Since K1 already includes inclination, no further projection is
    needed. Historical porb/m1 arguments are redundant and accepted only for
    compatibility; positive values are validated when supplied.
    """
    q = _positive_scalar(q, "q")
    k1 = _positive_scalar(k1, "k1")
    if int(j) != j or int(k) != k or not 0 < k < j:
        raise ValueError("Resonance indices must be integers with 0 < k < j.")
    if porb is not None:
        _positive_scalar(porb, "porb")
    if m1 is not None:
        _positive_scalar(m1, "m1")
    resonance_radius = ((j - k) / j) ** (2 / 3) / (1 + q) ** (1 / 3)
    tidal_radius = 0.60 / (1 + q)
    _, plt = _lazy_import_matplotlib_pyplot()
    plt.figure(num="Doppler Map")
    phis = np.linspace(0, 2 * np.pi, 629)
    for radius, style in ((resonance_radius, ":"), (tidal_radius, "-")):
        velocity = k1 * np.sqrt((1 + q) / radius) / q
        xx, yy = velocity * np.cos(phis), velocity * np.sin(phis)
        plt.plot(xx, yy - k1, c="k", ls=style)
    plt.draw()


def test_data(
    destination: Optional[Union[Path, str]] = None,
    overwrite: bool = False,
) -> List[Path]:
    """Copy the bundled test data into *destination* (defaults to CWD)."""

    dest_dir = Path.cwd() if destination is None else Path(destination)
    copied = copy_test_data(dest_dir, overwrite=overwrite)

    LOGGER.info("Copied %d test data files to %s", len(copied), dest_dir)
    return copied

