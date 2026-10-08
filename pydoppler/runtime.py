"""Validated, provenance-tracked execution of the standalone tomography solver."""

from __future__ import annotations

import hashlib
import copy
import json
import logging
import os
import shutil
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np


def phase_mask(phases, start, stop):
    """Valid phases outside the open, periodically repeated exclusion interval."""
    if not np.isfinite([start, stop]).all() or not 0 <= stop - start < 1:
        raise ValueError("Phase exclusion must satisfy 0 <= pb1 - pb0 < 1 cycle.")
    distance = np.mod(np.asarray(phases, dtype=float) - start, 1.0)
    edge = 16 * np.finfo(float).eps
    return ~((distance > edge) & (distance < stop - start - edge))


def covariance_modes(covariance, rcond=1e-10):
    """Eigenmodes of a covariance, evaluated in scale-independent units.

    Null modes from interpolation/continuum subtraction carry no independent
    measurements. They are excluded from the generalized least-squares metric.
    The returned eigenvalues are in the original squared flux units.
    """
    if not np.isfinite(rcond) or not 0 < rcond < 1:
        raise ValueError("covariance_rcond must be strictly between zero and one.")
    c = covariance.to_dense() if hasattr(covariance, "to_dense") else covariance
    c = np.asarray(c, dtype=float)
    if c.ndim != 2 or c.shape[0] != c.shape[1] or not np.isfinite(c).all():
        raise ValueError("Covariance must be a finite square matrix.")
    scale = float(np.max(np.abs(np.diag(c))))
    if scale <= 0:
        raise ValueError("Covariance has no positive variance.")
    if not np.allclose(c / scale, c.T / scale, rtol=1e-10, atol=1e-12):
        raise ValueError("Covariance must be symmetric.")
    values, vectors = np.linalg.eigh((c + c.T) / (2 * scale))
    if values[0] < -1e-9 * max(values[-1], 1.0):
        raise ValueError("Covariance must be positive semidefinite.")
    selected = values > rcond * values[-1]
    if not np.any(selected):
        raise ValueError("Covariance has no resolved positive modes.")
    return values[selected] * scale, vectors[:, selected]


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_changed(path, content):
    path = Path(path)
    if not path.exists() or path.read_text() != content:
        path.write_text(content, encoding="utf-8")


def preprocessing_provenance(workdir, input_hash):
    """Return preprocessing metadata only when its numerical input still matches."""
    path = Path(workdir) / "dopin.json"
    try:
        metadata = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"status": "unavailable", "reason": "No dopin.json sidecar; input may be manual."}
    except (OSError, ValueError):
        return {"status": "unavailable", "reason": "The dopin.json sidecar is unreadable or malformed."}
    if not isinstance(metadata, dict) or metadata.get("dopin_sha256") != input_hash:
        return {"status": "unavailable", "reason": "Input changed since preprocessing metadata was written."}
    return {"status": "matched", "sidecar_sha256": _sha(path), "metadata": metadata}


@contextmanager
def _run_lock(workdir):
    path = workdir / ".pydoppler.lock"
    try:
        fd = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(
            f"{workdir} is already in use. Use a separate run directory. "
            f"If its previous process has exited, remove the stale lock {path}."
        ) from exc
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(str(os.getpid()))
        yield
    finally:
        path.unlink(missing_ok=True)


def read_input(path, weighted=False):
    """Read the numerical input contract independently of the mutable wrapper."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError("dopin is incomplete; run Dopin() first.")
    try:
        head = lines[0].split()
        nph, nvp = int(head[0]), int(head[1])
        lam0 = float(head[2])
        values = np.asarray(
            [float(s.replace("D", "e")) for s in " ".join(lines[2:]).split()]
        )
    except (ValueError, IndexError) as exc:
        raise ValueError("Malformed dopin numerical data.") from exc
    if nph < 1 or nvp < 3 or not np.isfinite(lam0) or lam0 <= 0:
        raise ValueError("dopin requires positive wavelength and at least three channels.")
    if not np.isfinite(values).all():
        raise ValueError("dopin contains NaN or infinite values.")
    phases = values[:nph]
    index = nph
    if values.size <= index or values[index] not in (0, 1):
        raise ValueError("Invalid exposure-width flag in dopin.")
    explicit = int(values[index])
    index += 1
    widths = values[index:index + nph] if explicit else np.full(nph, 1 / nph)
    index += explicit * nph
    vp = values[index:index + nvp]
    index += nvp
    expected = index + nph * nvp * (2 if weighted else 1)
    if values.size < expected:
        raise ValueError("dopin lacks the expected flux/error samples; rerun Dopin().")
    if widths.size != nph or np.any(widths <= 0) or np.any(widths > 1):
        raise ValueError("Exposure widths must lie in (0, 1] orbital cycles.")
    if vp.size != nvp or np.any(np.diff(vp) <= 0):
        raise ValueError("Velocity samples must be strictly increasing.")
    if not np.allclose(np.diff(vp), np.diff(vp)[0], rtol=5e-5, atol=0):
        raise ValueError("The Fortran backend requires uniformly spaced velocities.")
    flux = values[index:index + nph * nvp].reshape(nvp, nph)
    sigma = None
    if weighted:
        sigma = values[index + nph * nvp:expected].reshape(nvp, nph)
        if np.any(sigma <= 0):
            raise ValueError("All supplied standard deviations must be positive.")
    return phases, widths, vp, flux, sigma


def _prepare_precision(dop, flux, sigma):
    mode = dop.covariance_mode
    if mode not in ("full", "diagonal"):
        raise ValueError("covariance_mode must be 'full' or 'diagonal'.")
    covariances = getattr(dop, "normalised_flux_covariance", None)
    if dop.iw == 0 or (dop.iw == 1 and mode == "diagonal"):
        (dop.workdir / "dop.precision").unlink(missing_ok=True)
        return int(dop.iw)
    if covariances is None:
        if dop.iw == 2:
            raise ValueError("iw=2 requires the covariance prepared by Dopin().")
        return 1
    if len(covariances) != flux.shape[1]:
        raise ValueError("Prepared covariance does not match dopin; rerun Dopin().")
    prepared = np.asarray(dop.normalised_flux, dtype=float).T
    if prepared.shape != flux.shape or not np.allclose(prepared, flux, rtol=1e-10, atol=0):
        raise ValueError("dopin differs from the covariance's flux data; rerun Dopin().")
    ranks = []
    path = dop.workdir / "dop.precision"
    with path.open("w", encoding="ascii") as handle:
        handle.write(f"{flux.shape[1]} {flux.shape[0]}\n")
        for index, covariance in enumerate(covariances):
            diagonal = np.asarray(covariance.diagonal(), dtype=float)
            if not np.allclose(diagonal, sigma[:, index] ** 2, rtol=1e-9, atol=0):
                raise ValueError("Prepared covariance and dopin uncertainties differ.")
            values, vectors = covariance_modes(covariance, dop.covariance_rcond)
            precision = (vectors / values) @ vectors.T
            if not np.isfinite(precision).all():
                raise ValueError("Precision matrix overflow; change the flux units.")
            ranks.append(len(values))
            handle.write(f"{len(values)}\n")
            np.savetxt(handle, precision, fmt="%.16e")
    dop.covariance_ranks = np.asarray(ranks)
    return 2


def _parameters(dop, effective_iw):
    numbers = [dop.ac, dop.al0, dop.alf, dop.clim, dop.wid, dop.af,
               dop.instrument_sigma, dop.baseline_fraction, dop.overs]
    if not np.isfinite(numbers).all():
        raise ValueError("Solver parameters must be finite.")
    if dop.ac <= 0 or dop.nim < 1 or dop.clim <= 0 or dop.al0 < 0 or dop.alf <= 1:
        raise ValueError("Require ac>0, nim>=1, clim>0, al0>=0 and alf>1.")
    if dop.norm not in (0, 1) or dop.ih not in (0, 1) or dop.iw not in (0, 1, 2):
        raise ValueError("Invalid likelihood, uncertainty, or normalization mode.")
    if dop.norm == 1 and effective_iw > 0 and not dop.conditional_normalization:
        raise ValueError(
            "norm=1 estimates normalization factors from the data. Use norm=0 for "
            "unconditional noise propagation, or explicitly set conditional_normalization=True "
            "to condition inference on those fixed factors."
        )
    if dop.ns < 1 or dop.instrument_sigma < 0 or dop.baseline_fraction < 0:
        raise ValueError("ns must be positive; LSF sigma and baseline must be nonnegative.")
    if dop.phase_subsamples < 0 or int(dop.phase_subsamples) != dop.phase_subsamples:
        raise ValueError("phase_subsamples must be a nonnegative integer (0=adaptive).")
    if not 0 < dop.overs <= 1:
        raise ValueError("overs must be in (0, 1].")
    if (not np.isfinite([dop.run_timeout, dop.max_memory_bytes, dop.covariance_rcond]).all()
            or dop.run_timeout <= 0 or dop.max_memory_bytes <= 0
            or not 0 < dop.covariance_rcond < 1):
        raise ValueError("Invalid execution limit or covariance eigenvalue threshold.")
    for name in ("ns", "nim", "nal", "ipri"):
        if int(getattr(dop, name)) != getattr(dop, name):
            raise ValueError(f"{name} must be an integer.")
    lines = [
        f"{1 if effective_iw else int(dop.ih)}    ih (0=positive log objective, 1=least squares)",
        f"{effective_iw}    iw (0=unweighted, 1=diagonal, 2=full covariance)",
        f"{dop.pb0:.16e} {dop.pb1:.16e}    excluded open phase interval",
        f"{int(dop.ns)}    default-map smearing pixels",
        f"{dop.ac:.16e}    convergence tolerance",
        f"{int(dop.nim)}    maximum iterations per alpha",
        f"{dop.al0:.16e} {dop.alf:.16e} {int(dop.nal)}    alpha schedule",
        f"{dop.clim:.16e}    RMS target (not reduced chi-square)",
        f"{int(dop.ipri)}    log verbosity",
        f"{int(dop.norm)}    opt-in flat integrated flux normalization",
        f"{dop.wid:.16e} {dop.af:.16e}    absorption prescription",
        f"{dop.instrument_sigma * 1e5:.16e} {int(dop.phase_subsamples)} "
        f"{dop.baseline_fraction:.16e}    LSF sigma cm/s, phase quadrature, baseline",
    ]
    return "\n".join(lines) + "\n"


def _dimensions(nph, nvp, nv, nri, ndi):
    nt = nvp * nph + nv * nvp * 3 + 2 * nph * nv
    rows = max(1, int(nri * nv * nt / (nph * nvp)))
    columns = max(1, int(ndi * nv * nt / (nv * nv)))
    source = (
        f"      parameter (npm={nph},nvpm={nvp},nvm={nv})\n"
        "      parameter (nd=npm*nvpm,nr=nvm*nvm)\n"
        f"      parameter (nri={rows},ndi={columns})\n"
        "      parameter (nf=nd,nfp=ndi,nrp=nri)\n"
        "      parameter (ni=nvm,nj=nvm)\n"
        "      parameter (nch=1,nsu=1)\n"
    )
    return source, 8 * (rows * nph * nvp + columns * nv * nv)


def _command(args, cwd, timeout):
    try:
        result = subprocess.run(args, cwd=cwd, capture_output=True, text=True,
                                check=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"Tomography command exceeded {timeout} s: {args[0]}") from exc
    except subprocess.CalledProcessError as exc:
        status_path = Path(cwd) / "dop.status"
        if Path(args[0]).name == "dopp" and (Path(cwd) / "emap.par.new").exists():
            return subprocess.CompletedProcess(args, exc.returncode, exc.stdout, exc.stderr)
        if Path(args[0]).name == "dopp" and status_path.exists():
            _status(status_path)  # Prefer the scientific failure to a generic exit code.
        raise RuntimeError(
            f"Tomography command failed: {' '.join(args)}\n"
            f"{exc.stdout or ''}\n{exc.stderr or ''}"
        ) from exc
    return result


def _status(path):
    if not path.exists():
        raise RuntimeError("Missing dop.status: solver failed or workdir contains old Fortran assets. "
                           "Refresh them with copy_fortran_code(workdir, overwrite=True).")
    result = {}
    for line in path.read_text().splitlines():
        tokens = line.replace("=", " ").split()
        if len(tokens) >= 2:
            try:
                result[tokens[0]] = float(tokens[1].replace("D", "e"))
            except ValueError:
                result[tokens[0]] = " ".join(tokens[1:])
    required = ("converged", "optimizer_status", "fit_statistic", "fit_target", "iterations")
    if any(key not in result for key in required):
        raise RuntimeError("Incomplete dop.status; cannot validate solver result.")
    for key in required:
        if not isinstance(result[key], (int, float)) or not np.isfinite(result[key]):
            raise RuntimeError(f"Invalid/nonfinite solver status field: {key}.")
    for key, value in result.items():
        if isinstance(value, (int, float)) and not np.isfinite(value):
            raise RuntimeError(f"Invalid/nonfinite solver status field: {key}.")
    if result["optimizer_status"] == 1:
        raise RuntimeError(
            "Optimizer converged, but the alpha schedule did not reach the RMS target: "
            f"{result['fit_statistic']:.6g} > {result['fit_target']:.6g}. See dop.log."
        )
    if result["converged"] != 1 or result["optimizer_status"] != 0:
        raise RuntimeError(f"Optimizer did not converge: {result}. See dop.log.")
    if result["fit_statistic"] > result["fit_target"]:
        raise RuntimeError(
            f"Converged optimizer did not reach the RMS target: "
            f"{result['fit_statistic']:.6g} > {result['fit_target']:.6g}. See dop.log."
        )
    return result


def run(dop, nri=0.9, ndi=0.7):
    """Compile by content fingerprint, always execute, then validate the result."""
    from .pydoppler import _parse_dopout

    workdir = Path(dop.workdir)
    for name in ("dopin", "dop.f", "makefile", "cclock.f"):
        if not (workdir / name).is_file():
            raise FileNotFoundError(f"Missing {workdir / name}; install sources and run Dopin().")
    for executable in ("make", "gfortran"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"Required executable {executable!r} was not found on PATH.")
    if not np.isfinite([nri, ndi]).all() or min(nri, ndi) <= 0:
        raise ValueError("Projection allocation factors must be positive and finite.")
    with _run_lock(workdir):
        started = time.time()
        dop.fit_diagnostics = None
        dop._output_hash = None
        dop._fit_sigma = None
        dop._fit_covariance = None
        # Invalidate previous provenance before any operation that may fail.
        (workdir / "run.json").write_text(json.dumps(
            {"started_unix": started, "success": False, "state": "preparing"},
            indent=2) + "\n")
        phases, widths, vp, flux, sigma = read_input(workdir / "dopin", dop.iw > 0)
        valid = phase_mask(phases, dop.pb0, dop.pb1)
        if not np.any(valid) or not np.any(flux[:, valid] != 0):
            raise ValueError("There must be nonzero data in at least one unmasked spectrum.")
        dop.valid_phase_mask = valid
        full_covariance = (dop.iw == 2 or (
            dop.iw == 1 and dop.covariance_mode == "full"
            and getattr(dop, "normalised_flux_covariance", None) is not None))
        effective_iw = 2 if full_covariance else int(dop.iw)
        _parameters(dop, effective_iw)  # Validate before allocation or serialization.
        nv = dop.map_size if dop.map_size is not None else max(4, int(dop.overs * len(vp)))
        if not np.isfinite(nv) or int(nv) != nv or nv < 4:
            raise ValueError("map_size must be an integer of at least four pixels.")
        nv = int(nv) + int(nv) % 2
        _, projected_memory = _dimensions(len(phases), len(vp), nv, nri, ndi)
        covariance_memory = 8 * len(vp)**2 * (len(phases) + 6) if full_covariance else 0
        if projected_memory + covariance_memory > dop.max_memory_bytes:
            raise MemoryError(
                "Projection/covariance workspace exceeds max_memory_bytes before preparation; "
                "reduce map_size or the spectral window, or increase the limit deliberately."
            )
        effective_iw = _prepare_precision(dop, flux, sigma)
        dop.effective_iw = effective_iw
        parameters = _parameters(dop, effective_iw)
        _write_changed(workdir / "dop.in", parameters)
        compiler = _command(["gfortran", "--version"], workdir, dop.run_timeout).stdout
        manifest_path = workdir / "run.json"
        manifest = {"started_unix": started, "success": False, "parameters": parameters,
                    "effective_iw": effective_iw, "map_size": nv,
                    "covariance_rcond": dop.covariance_rcond,
                    "input_sha256": _sha(workdir / "dopin"), "compiler": compiler.splitlines()[0]}
        manifest["preprocessing"] = preprocessing_provenance(workdir, manifest["input_sha256"])
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        dop.fit_diagnostics = None
        try:
            for attempt in range(4):
                dimensions, memory = _dimensions(len(phases), len(vp), nv, nri, ndi)
                if effective_iw == 2:
                    memory += 8 * len(phases) * len(vp) ** 2
                if memory > dop.max_memory_bytes:
                    raise MemoryError(
                        f"Projection/covariance allocation needs approximately {memory / 1e9:.2f} GB; "
                        "reduce map_size/window or increase max_memory_bytes deliberately."
                    )
                _write_changed(workdir / "emap.par", dimensions)
                files = ("dop.f", "emap.par", "makefile", "cclock.f")
                source_hashes = {name: _sha(workdir / name) for name in files}
                fingerprint = hashlib.sha256(json.dumps(
                    {"sources": source_hashes, "compiler": compiler,
                     "compiler_path": shutil.which("gfortran")}, sort_keys=True).encode()).hexdigest()
                build_path = workdir / ".build.json"
                try:
                    build = json.loads(build_path.read_text())
                except (OSError, ValueError):
                    build = {}
                binary = workdir / "dopp"
                if (build.get("fingerprint") != fingerprint or not binary.exists()
                        or build.get("binary_sha256") != _sha(binary)):
                    dop._log(logging.INFO, "Compiling tomography solver.")
                    _command(["make", "-B", "dopp"], workdir, dop.run_timeout)
                    build_path.write_text(json.dumps({"fingerprint": fingerprint,
                                                     "binary_sha256": _sha(binary)}) + "\n")
                for name in ("dop.out", "dop.log", "dop.status", "dop.normalization", "emap.par.new"):
                    (workdir / name).unlink(missing_ok=True)
                result = _command([str(binary.resolve())], workdir, dop.run_timeout)
                if result.stderr:
                    dop._log(logging.DEBUG, result.stderr.rstrip())
                if (workdir / "emap.par.new").exists():
                    nri *= 2
                    ndi *= 2
                    continue
                break
            else:
                raise RuntimeError("Projection allocation failed after four bounded attempts.")
            diagnostics = _status(workdir / "dop.status")
            output = _parse_dopout(workdir / "dop.out")
            if output.nph != len(phases) or output.nvp != len(vp) or output.nv != nv:
                raise RuntimeError("Output dimensions do not match this run's inputs.")
            for name in ("im", "dm", "dmr", "vp", "pha", "dpha"):
                if not np.isfinite(getattr(output, name)).all():
                    raise RuntimeError(f"Nonfinite values in solver output {name}.")
            if not np.allclose(output.vp, vp / 1e5, rtol=1e-12, atol=1e-9):
                raise RuntimeError("Output velocity grid does not match this run.")
            if (not np.allclose(output.pha, phases, rtol=1e-12, atol=1e-13)
                    or not np.allclose(output.dpha, widths, rtol=1e-12, atol=1e-13)):
                raise RuntimeError("Output phases/exposure widths do not match this run.")
            if not np.allclose(output.dm, flux, rtol=1e-11, atol=np.max(np.abs(flux)) * 1e-12):
                raise RuntimeError("Output data differ from input; refusing untracked preprocessing.")
            dop.fit_diagnostics = diagnostics
            dop._output_input_hash = manifest["input_sha256"]
            dop._output_hash = _sha(workdir / "dop.out")
            dop._fit_sigma = None if sigma is None else sigma.copy()
            # Only the full-covariance path has validated this preparation
            # against both the serialized data and its supplied uncertainties.
            # In diagonal/manual-input runs, current wrapper covariance may
            # belong to a previous preparation and must not be reused.
            dop._fit_covariance = (
                copy.deepcopy(dop.normalised_flux_covariance) if effective_iw == 2 else None
            )
            dop._fit_rcond = dop.covariance_rcond
            dop._fit_phase_mask = valid.copy()
            manifest.update(success=True, diagnostics=diagnostics, source_sha256=source_hashes,
                            output_sha256=dop._output_hash, binary_sha256=_sha(binary))
            if effective_iw == 2:
                manifest["precision_sha256"] = _sha(workdir / "dop.precision")
            dop._log(logging.INFO, f"Converged tomography; RMS statistic={diagnostics['fit_statistic']:.6g}.")
            return diagnostics
        except Exception as exc:
            manifest["error"] = str(exc)
            raise
        finally:
            manifest["elapsed_seconds"] = time.time() - started
            manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
