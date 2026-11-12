from pathlib import Path

import pydoppler


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


def test_copy_test_data_respects_overwrite(tmp_path: Path):
    script = tmp_path / "sample_script.py"
    script.write_text("dummy", encoding="utf-8")

    copied = pydoppler.copy_test_data(tmp_path, overwrite=False)
    assert script.read_text(encoding="utf-8") == "dummy"
    assert script not in copied

    pydoppler.copy_test_data(tmp_path, overwrite=True)
    assert script.read_text(encoding="utf-8") != "dummy"


def test_install_assets_adds_the_sample_script_once(tmp_path: Path):
    dop = pydoppler.spruit(auto_install=True, workdir=tmp_path)
    script = tmp_path / "sample_script.py"
    assert script.is_file()

    script.write_text("edited", encoding="utf-8")
    dop.install_assets(tmp_path)
    assert script.read_text(encoding="utf-8") == "edited"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["sample_script.py"]

    dop.install_assets(tmp_path, force=True)
    assert script.read_text(encoding="utf-8") != "edited"
