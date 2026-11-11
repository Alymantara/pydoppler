from pathlib import Path

import pydoppler


def test_copy_fortran_code_copies_expected_files(tmp_path: Path):
    copied = pydoppler.copy_fortran_code(tmp_path)

    assert (tmp_path / "dop.f90").is_file()
    assert (tmp_path / "makefile").is_file()
    assert (tmp_path / "dop_core.f90").is_file()
    assert not (tmp_path / "emap.par").exists()  # sizes are set at run time
    assert all(path.is_file() for path in copied)


def test_copy_fortran_code_respects_overwrite(tmp_path: Path):
    dop_f = tmp_path / "dop.f90"
    dop_f.write_text("dummy", encoding="utf-8")

    copied = pydoppler.copy_fortran_code(tmp_path, overwrite=False)
    assert dop_f.read_text(encoding="utf-8") == "dummy"
    assert dop_f not in copied

    pydoppler.copy_fortran_code(tmp_path, overwrite=True)
    assert dop_f.read_text(encoding="utf-8") != "dummy"


def test_fortran_makefile_enables_checked_builds_without_suppressing_errors(tmp_path: Path):
    pydoppler.copy_fortran_code(tmp_path)

    makefile = (tmp_path / "makefile").read_text(encoding="utf-8")

    assert "FC = gfortran" in makefile
    assert "dop_core.f90 dop.f90" in makefile
    assert "-fcheck=all" in makefile
    assert "FFLAGS = -O -w" not in makefile
    assert "-fallow-argument-mismatch" not in makefile
    assert "$(FC) $(FFLAGS)" in makefile


def test_install_sample_script_avoids_overwrite(tmp_path: Path):
    first = pydoppler.install_sample_script(tmp_path)
    assert first.name == "sample_script.py"
    assert first.is_file()

    second = pydoppler.install_sample_script(tmp_path)
    assert second.name.startswith("sample_script-")
    assert second.suffix == ".py"
    assert second.is_file()


def test_copy_test_data_copies_dataset(tmp_path: Path):
    copied = pydoppler.copy_test_data(tmp_path)

    assert (tmp_path / "ugem99" / "ugem0all.fas").is_file()
    assert (tmp_path / "sample_script.py").is_file()
    assert all(path.is_file() for path in copied)
