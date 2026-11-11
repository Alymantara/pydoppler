from pathlib import Path

import pydoppler


def test_set_workdir_updates_base_and_creates_directory(tmp_path: Path):
    dop = pydoppler.spruit(workdir=tmp_path / "pydoppler-workdir", interactive=False)

    new_base = tmp_path / "runs"
    resolved = dop.set_workdir(new_base)

    assert resolved == new_base.resolve()
    assert dop.workdir == new_base.resolve()
    assert dop.workdir.is_dir()

    run_dir = dop.make_run_dir(prefix="run")
    assert run_dir == (new_base / "run").resolve()
    assert (run_dir / "dop.f90").is_file()


def test_make_run_dir_creates_unique_directories(tmp_path: Path):
    base = tmp_path / "pydoppler-workdir"
    dop = pydoppler.spruit(workdir=base, interactive=False)

    run1 = dop.make_run_dir(prefix="run")
    assert run1 == (base / "run").resolve()
    assert (run1 / "makefile").is_file()

    run2 = dop.make_run_dir(prefix="run")
    assert run2 == (base / "run-1").resolve()
    assert (run2 / "dop.f90").is_file()

