import shutil
import subprocess
import sys
from pathlib import Path


def test_import_does_not_import_matplotlib():
    code = """
import sys
sys.modules.pop("matplotlib", None)
import pydoppler
assert "matplotlib" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_import_from_checkout_parent_exposes_api(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    checkout = tmp_path / "pydoppler"
    checkout.mkdir()

    shutil.copyfile(repo_root / "__init__.py", checkout / "__init__.py")
    shutil.copytree(
        repo_root / "pydoppler",
        checkout / "pydoppler",
        ignore=shutil.ignore_patterns("__pycache__", "*.py[cod]"),
    )

    code = """
import pydoppler
assert hasattr(pydoppler, "copy_test_data")
assert hasattr(pydoppler, "spruit")
"""
    subprocess.run([sys.executable, "-c", code], check=True, cwd=tmp_path)
