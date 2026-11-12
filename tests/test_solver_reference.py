"""The NumPy solver against stored results of the Fortran solver it replaced.

``tests/data/fortran_reference.npz`` holds small synthetic problems together
with the output of the Fortran solver (``dop_core.f90`` and ``dop.f90``) at the
commit, source hashes and compiler recorded in its metadata:

- the final status, map and reconstructed trail of complete inversions;
- for some cases, the building blocks (projection, baseline, likelihood,
  entropy, Newton step) evaluated at stored test vectors.

Building blocks do not depend on the optimizer's path and are compared to
rounding accuracy. Complete inversions agree to rounding on the machine that
generated the data; the tolerances leave room for other numerical libraries.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pydoppler import solver

REFERENCE = np.load(Path(__file__).parent / "data" / "fortran_reference.npz")
CASES = json.loads(str(REFERENCE["metadata"]))["cases"]
VELOCITY = REFERENCE["velocity"]

EXACT = ("converged", "optimizer_status", "alpha_steps", "target_bracketed", "rms_plateau",
         "statistic_rank", "baseline_raised", "likelihood", "weight_mode")
CLOSE = ("fit_statistic", "alpha", "alpha_start", "entropy", "objective", "flux_scale",
         "information", "baseline", "rms_floor_estimate")


def _inputs(name):
    case = CASES[name]
    settings = case["settings"]
    flux = REFERENCE[f"{name}/flux"]
    nph, nvp = flux.shape
    sigma = np.full(flux.shape, case["sigma"]) if settings["iw"] else None
    precision = ranks = None
    if case["has_precision"]:
        precision = np.repeat(REFERENCE[f"{name}/precision"][None], nph, axis=0)
        ranks = np.full(nph, nvp)
    arguments = (REFERENCE[f"{name}/phases"], REFERENCE[f"{name}/widths"], VELOCITY, flux,
                 sigma, precision, ranks)
    return arguments, settings


@pytest.mark.parametrize("name", sorted(CASES))
def test_inversion_matches_fortran(name):
    arguments, settings = _inputs(name)
    expected = CASES[name]["status"]
    result = solver.solve(*arguments, **settings)
    for key in EXACT:
        assert result.status[key] == expected[key], key
    for key in CLOSE:
        assert result.status[key] == pytest.approx(expected[key], rel=1e-6, abs=1e-12), key
    assert abs(result.status["total_iterations"] - expected["total_iterations"]) <= max(
        2, 0.1*expected["total_iterations"])
    if not expected["converged"]:
        assert CASES[name]["fortran_exit"] != 0
        assert result.image is None and result.model is None
        return
    assert result.status["stationarity"] <= settings["ac"]
    image, model = REFERENCE[f"{name}/image"], REFERENCE[f"{name}/model"]
    np.testing.assert_allclose(result.image, image, rtol=0, atol=1e-6*np.max(np.abs(image)))
    np.testing.assert_allclose(result.model, model, rtol=0, atol=1e-6*np.max(np.abs(model)))
    flux = arguments[3]
    np.testing.assert_allclose(result.data, flux, rtol=0, atol=1e-12*np.max(np.abs(flux)))
    assert result.identification == pytest.approx(CASES[name]["identification"], abs=1e-3)


@pytest.mark.parametrize("name", sorted(name for name, case in CASES.items() if case["kernel"]))
def test_building_blocks_match_fortran(name):
    arguments, settings = _inputs(name)
    expected = CASES[name]["kernel"]

    def stored(key):
        return REFERENCE[f"{name}/kernel/{key}"]

    def close(actual, reference, rtol=1e-11):
        reference = np.asarray(reference, dtype=float)
        np.testing.assert_allclose(actual, reference, rtol=rtol,
                                   atol=rtol*np.max(np.abs(reference)))

    problem = solver.Tomography(*arguments, **settings)
    problem.build_projection()
    data = problem.dat.ravel().copy()
    baseline, baseline_map = problem.baseline(data)
    target = data + baseline
    problem.set_weights()

    # Grid, sparsity pattern and statistical scales.
    assert problem.nol == expected["nol"]
    if settings["instrument_sigma"] > 0:
        # Whether a far Gaussian tail (~1e-20) rounds to zero depends on erfc.
        assert problem.matrix.nnz == pytest.approx(expected["nnz"], rel=1e-3)
    else:
        assert problem.matrix.nnz == expected["nnz"]
    assert problem.statistic_rank == expected["statistic_rank"]
    assert problem.ih == expected["ih"]
    np.testing.assert_array_equal(problem.indx, stored("indx"))
    np.testing.assert_array_equal(problem.indy, stored("indy"))
    for value, key in ((problem.va, "va"), (problem.dvmap, "dvmap"),
                       (problem.flux_scale, "flux_scale"), (problem.information, "information"),
                       (problem.baseline_level, "baseline_level")):
        assert value == pytest.approx(expected[key], rel=1e-12), key
    close(target, stored("target").ravel())
    close(baseline, stored("baseline").ravel())
    close(baseline_map, stored("vmb"))

    # Projection and its adjoint, likelihood, entropy, and the Newton step.
    ps, d1, d2 = stored("ps"), stored("d1"), stored("d2")
    forward = problem.fproj(ps)
    close(forward, stored("forward").ravel())
    close(problem.bproj(stored("t").ravel()), stored("back"))
    likelihood, gradient = problem.like(target, forward)
    assert likelihood == pytest.approx(expected["likelihood"], rel=1e-11)
    close(gradient, stored("dha"))
    entropy, entropy_gradient, default, logratio = problem.entr(ps)
    assert entropy == pytest.approx(expected["entropy"], rel=1e-11)
    close(default, stored("chi"))
    close(entropy_gradient, stored("dsa"))
    fla, fmu = problem.lamun(target, ps, d1, d2, forward, default, logratio, 0.7)
    assert (fla, fmu) == pytest.approx((expected["fla"], expected["fmu"]), rel=1e-9)
    assert problem.rrms(target, forward) == pytest.approx(expected["rms"], rel=1e-11)
