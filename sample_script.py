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

    # Import sample data
    # <<< COMMENT OUT IF YOU DONT NEED THE TEST DATASET >>>
    pydoppler.copy_test_data(workdir, overwrite=True)

    pydoppler.copy_fortran_code(workdir, overwrite=True)

    # Load base object for tomography
    dop = pydoppler.spruit(workdir=workdir)

    # Basic data for the tomography to work
    dop.object = 'U Gem'
    dop.base_dir = workdir / 'ugem99' # Base directory for input spectra
    dop.list = 'ugem0all.fas'		# Name of the input file
    dop.lam0 = 6562.8  # Rest wavelength in Angstrom
    dop.delta_phase = 0.003  # Full exposure width in orbital cycles
    dop.delw = 35  # Rest-frame line-window half-width in Angstrom
    dop.overs = 0.3  # Map pixels per axis / spectral velocity-channel count
    dop.gama = 36.0  # Systemic redshift in km/s, removed before inversion
    dop.nbins = 28  # Diagnostic phase bins; the inversion uses all exposures

    # Read in the individual spectra and orbital phase information
    dop.Foldspec()

    # Subtract the continuum and resample onto a common velocity grid.
    # This bundled dataset has no uncertainties: its fit is illustrative.
    # Supply measured flux errors for weighted fits and calibrated diagnostics.
    dop.Dopin(
        continuum_band=[6500, 6537, 6591, 6620],
        plot_median=False,
        poly_degree=2,
        show=False,
    )

    # Perform tomography
    dop.Syncdop()

    # This routine will display the outcome of the Doppler tomography.
    # You can overplot contours and streams.
    _cbar, _data = dop.Dopmap(
        limits=[0.05, 0.99],
        colorbar=True,
        cmaps=plt.cm.magma,
        smooth=False,
        remove_mean=False,
        show=False,
    )

    # Overlay projected donor, Keplerian and ballistic velocities.
    # K1 fixes the velocity scale; M1/P/q/i provide a consistency check.
    qm=0.35  # M2/M1
    k1 = 107  # Observed primary semi-amplitude in km/s
    inc=70  # Degrees
    m1=1.2  # Solar masses
    porb=0.1769061911  # Days

    pydoppler.stream(qm,k1,porb,m1,inc)

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
        for stale_name in (
            "Average_Spec.pdf",
            "Residuals.pdf",
            "Average_Spec.png",
            "Residuals.png",
            "Trail.png",
            "Doppler_Map.png",
            "Reconstruction.png",
        ):
            (outdir / stale_name).unlink(missing_ok=True)
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
