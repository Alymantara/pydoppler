"""Regenerate the README figures and their provenance record.

Run from the repository root:  python docs/make_readme_figures.py

The sample script runs in a temporary directory; its figures are saved as SVG
with rasterized image/mesh layers (vector axes and labels), and
docs/readme_images/provenance.json records the hashes of the sources and outputs.
"""

import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "docs" / "readme_images"
FIGURES = {"Trail": "Trail.svg", "Doppler Map": "Doppler_Map.svg",
           "Reconstruction": "Reconstruction.svg"}
SOURCES = ("sample_script.py", "pydoppler/pydoppler.py", "pydoppler/preprocessing.py",
           "pydoppler/runtime.py", "pydoppler/fortran_code/dop_core.f90",
           "pydoppler/fortran_code/dop.f90")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    os.environ["MPLBACKEND"] = "Agg"
    sys.path.insert(0, str(ROOT))
    import matplotlib

    matplotlib.use("Agg", force=True)
    matplotlib.rcParams["svg.hashsalt"] = "pydoppler"  # deterministic element ids
    matplotlib.rcParams["svg.fonttype"] = "path"
    import matplotlib.pyplot as plt
    import sample_script

    sample_script.SAVE_PDFS = False
    sample_script.SHOW_PLOTS = True
    plt.show = lambda *args, **kwargs: None  # keep figures open after main()
    with tempfile.TemporaryDirectory() as tmp:
        cwd = Path.cwd()
        os.chdir(tmp)
        try:
            sample_script.main()
            run_output = sha256(Path(tmp) / "pydoppler-workdir" / "dop.out")
        finally:
            os.chdir(cwd)
    outputs = {}
    for label, name in FIGURES.items():
        fig = plt.figure(label)
        for ax in fig.axes:
            for artist in [*ax.images, *ax.collections]:
                if artist.__class__.__name__ in ("AxesImage", "QuadMesh"):
                    artist.set_rasterized(True)
        path = OUTDIR / name
        fig.savefig(path, format="svg", dpi=110, bbox_inches="tight",
                    metadata={"Date": None})
        text = "\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n"
        path.write_text(text)
        outputs[f"docs/readme_images/{name}"] = sha256(path)
    provenance = {
        "generator": "docs/make_readme_figures.py (sample_script.main(), Agg backend)",
        "run_output_sha256": run_output,
        "sources": {name: sha256(ROOT / name) for name in SOURCES},
        "outputs": outputs,
    }
    (OUTDIR / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
