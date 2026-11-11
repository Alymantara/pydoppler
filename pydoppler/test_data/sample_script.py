"""PyDoppler tutorial: U Gem H-alpha Doppler tomography (Echevarria et al. 2007).

Run from a terminal with ``python sample_script.py`` (or ``%run sample_script.py``
in IPython). Everything is written to ./pydoppler-workdir.
"""

import os
from pathlib import Path

import pydoppler

# Save PDF figures; SHOW_PLOTS below controls interactive display.
SAVE_PDFS = True
# Set to True to display figures interactively at the end.
SHOW_PLOTS = False


def main() -> None:
    # Fortran execution + outputs are isolated in this directory.
    workdir = Path.cwd() / "pydoppler-workdir"
    if SAVE_PDFS and not SHOW_PLOTS:
        os.environ.setdefault("MPLBACKEND", "Agg")
        os.environ.setdefault("MPLCONFIGDIR", str(workdir / ".mplconfig"))

    import matplotlib.pyplot as plt

    # Copy the bundled U Gem spectra and the solver sources.
    # <<< COMMENT OUT IF YOU DONT NEED THE TEST DATASET >>>
    pydoppler.copy_test_data(workdir, overwrite=True)
    pydoppler.copy_fortran_code(workdir, overwrite=True)

    dop = pydoppler.spruit(workdir=workdir)

    # Basic data for the tomography to work
    dop.object = 'U Gem'
    dop.base_dir = workdir / 'ugem99'  # Directory with spectra and phase file
    dop.list = 'ugem0all.fas'  # Phase file: name, mid-exposure phase, exposure width
    dop.lam0 = 6562.8  # Rest wavelength in Angstrom
    dop.delw = 35  # Rest-frame line-window half-width in Angstrom
    dop.overs = 0.3  # Map pixels per axis / spectral velocity-channel count
    dop.gama = 36.0  # Systemic redshift in km/s, removed before inversion
    dop.nbins = 28  # Diagnostic phase bins; the inversion uses all exposures

    # Read in the individual spectra and orbital phase information
    dop.Foldspec()

    # Subtract the continuum and resample onto a common velocity grid.
    # The phase file's third column holds each exposure width (600 s =
    # 0.0393 cycles, Echevarria et al. 2007), so use it instead of delta_phase.
    # The bundled spectra have no uncertainties, so this fit is unweighted
    # (iw=0); see the README for a noise-weighted fit.
    dop.Dopin(
        continuum_band=[6500, 6537, 6591, 6620],
        plot_median=False,
        poly_degree=2,
        use_list_dpha=True,
        show=False,
    )

    # Perform tomography
    dop.Syncdop()

    # The map is shown as a fraction of its peak; zero flux is zero.
    _cbar, _data = dop.Dopmap(
        colorbar=True,
        cmaps=plt.cm.magma,
        smooth=False,
        remove_mean=False,
        show=False,
    )

    # Overlay projected donor, Keplerian and ballistic velocities.
    # K1 fixes the velocity scale; M1/P/q/i provide a consistency check.
    qm = 0.35  # M2/M1
    k1 = 107  # Observed primary semi-amplitude in km/s
    inc = 70  # Degrees
    m1 = 1.2  # Solar masses
    porb = 0.1769061911  # Days

    pydoppler.stream(qm, k1, porb, m1, inc)

    # Compare input and reconstructed trails on one common flux scale.
    # Visual agreement alone is not a significance or uncertainty estimate;
    # use residual diagnostics and synthetic recovery tests for interpretation.
    _cbar2, _cbar3, dmr, dm = dop.Reco(
        colorbar=True,
        limits=None,  # Shared scale includes the data and model ranges
        cmaps=plt.cm.magma,
        show=False,
    )

    if SAVE_PDFS:
        outdir = workdir / "output_images"
        outdir.mkdir(parents=True, exist_ok=True)
        figures = {
            "Trail": "Trail.pdf",
            "Doppler Map": "Doppler_Map.pdf",
            "Reconstruction": "Reconstruction.pdf",
        }
        existing = set(plt.get_figlabels())
        for label, filename in figures.items():
            if label not in existing:
                continue
            fig = plt.figure(label)
            fig.savefig(outdir / filename, dpi=200, bbox_inches="tight")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
