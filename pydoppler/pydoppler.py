"""
pydoppler.py  —  Native Python Doppler tomography (Spruit 1998).

No Fortran compiler required.
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import savgol_filter

from .mem_engine import doppler_mem, MEMResult

plt.rcParams.update({'font.size': 12})

# ─────────────────────────────────────────────────────────────────────────────
# Terminal colour helpers
# ─────────────────────────────────────────────────────────────────────────────

_BOLD  = "\033[1m"
_DIM   = "\033[2m"
_GREEN = "\033[32m"
_CYAN  = "\033[36m"
_YELLOW= "\033[33m"
_RED   = "\033[31m"
_RESET = "\033[0m"
_CHECK = "✓"
_CROSS = "✗"
_ARROW = "→"
_WAVE  = "~"

def _header(title: str):
    w = 60
    print(f"\n{_BOLD}{_CYAN}{'─'*w}{_RESET}")
    print(f"{_BOLD}{_CYAN}  {title}{_RESET}")
    print(f"{_BOLD}{_CYAN}{'─'*w}{_RESET}")

def _info(msg: str):
    print(f"  {_DIM}{_ARROW}{_RESET}  {msg}")

def _ok(msg: str):
    print(f"  {_GREEN}{_CHECK}{_RESET}  {msg}")

def _warn(msg: str):
    print(f"  {_YELLOW}!{_RESET}  {msg}")

def _progress_bar(current: int, total: int, label: str = "", width: int = 30) -> str:
    frac   = current / max(total, 1)
    filled = int(frac * width)
    bar    = "█" * filled + "░" * (width - filled)
    pct    = int(frac * 100)
    return f"  [{bar}] {pct:3d}%  {label}"


# ─────────────────────────────────────────────────────────────────────────────
# rebin_trail  (unchanged from original, but cleaned up)
# ─────────────────────────────────────────────────────────────────────────────

def rebin_trail(waver: np.ndarray, flux: np.ndarray,
                input_phase: np.ndarray,
                nbins: int, delp: float,
                rebin_wave=None) -> tuple[np.ndarray, np.ndarray]:
    """Phase-bin a trail spectrum onto a uniform phase grid."""
    phase     = np.linspace(0, 2, nbins * 2 + 1, endpoint=True) - 1.0 / nbins / 2.0
    phase     = np.concatenate((phase, [2.0 + 1.0 / nbins / 2.0]))
    phase_dec = phase - np.floor(phase)

    trail = np.zeros((waver.size, phase.size))
    tots  = trail.copy()

    for i in range(input_phase.size):
        dist = phase_dec - (input_phase[i] + delp / 2.0)
        dist[np.abs(dist) > 1.0 / nbins] = 0.0
        dist[dist > 0] = 0.0
        weights = np.abs(dist) / (1.0 / nbins)

        dist2 = phase_dec - (input_phase[i] - delp / 2.0)
        dist2[np.abs(dist2) > 1.0 / nbins] = 0.0
        dist2[dist2 > 0] = 0.0
        dist2[np.abs(dist2) > 0] = 1.0 - np.abs(dist2[np.abs(dist2) > 0]) / (1.0 / nbins)
        weights += dist2

        temp = trail.copy().T
        for j in range(phase.size):
            if rebin_wave is None:
                temp[j] = flux[i] * weights[j]
            else:
                temp[j] = np.interp(waver, waver, flux[i]) * weights[j]
        trail += temp.T
        tots  += weights

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        trail = np.where(tots > 0, trail / tots, 0.0)
    return trail, phase


# ─────────────────────────────────────────────────────────────────────────────
# Main class
# ─────────────────────────────────────────────────────────────────────────────

class spruit:
    """
    PyDoppler — Doppler tomography via Maximum Entropy Method.

    Pure Python / NumPy implementation: no Fortran compiler required.

    Typical usage
    -------------
    >>> dop = spruit()
    >>> dop.base_dir = 'ugem99'
    >>> dop.list     = 'ugem0all.fas'
    >>> dop.lam0     = 6562.8
    >>> dop.Foldspec()
    >>> dop.Dopin(continnum_band=[6500,6537,6591,6620])
    >>> dop.Syncdop()          # runs MEM in pure Python
    >>> cbar, data = dop.Dopmap()
    """

    def __init__(self):
        # ── Target / data ──────────────────────────────────────────────────
        self.object      = 'disc'
        self.base_dir    = '.'
        self.list        = 'phases.txt'
        self.lam0        = 6562.83      # rest wavelength (Å)
        self.delw        = 80           # half-window in wavelength (Å)
        self.gama        = 0.0          # systemic velocity (km/s)
        self.delta_phase = 0.001        # exposure half-width in phase
        self.nbins       = 20

        # ── MEM parameters ─────────────────────────────────────────────────
        self.ih     = 0       # 0=log-likelihood, 1=chi²
        self.iw     = 0       # 0=no error bars, 1=use error bars
        self.pb0    = 0.95    # blank phase range start
        self.pb1    = 1.05    # blank phase range end
        self.ns     = 7       # Gaussian smoothing width for default map
        self.ac     = 8e-4    # convergence criterion
        self.nim    = 150     # max iterations per alpha
        self.al0    = 0.002   # starting alpha
        self.alf    = 1.7     # alpha reduction factor
        self.nal    = 0       # number of alphas (0 → single)
        self.clim   = 1.6     # target reduced chi²
        self.norm   = 1       # normalise spectra
        self.wid    = 1e10    # central absorption width (cm/s)
        self.af     = 0.0     # central absorption amplitude
        self.overs  = 0.3     # velocity map oversampling
        self.ipri   = 2

        # ── Plot options ───────────────────────────────────────────────────
        self.psname  = 'doppler'
        self.output  = 'pdf'
        self.plotlim = 1.3
        self.lobcol  = 'white'
        self.verbose = True

        # ── Internal state ─────────────────────────────────────────────────
        self.wave           = []
        self.flux           = []
        self.pha            = None
        self.input_files    = None
        self.input_phase    = None
        self.trsp           = None
        self.normalised_flux = None
        self.normalised_wave = None
        self.vell            = None

        # Result from Syncdop
        self._mem_result: Optional[MEMResult] = None


    # ──────────────────────────────────────────────────────────────────────────
    # Foldspec
    # ──────────────────────────────────────────────────────────────────────────

    def Foldspec(self):
        """Read and phase-bin input spectra."""
        _header("Foldspec — reading & phase-binning spectra")

        list_path = os.path.join(self.base_dir, self.list)
        try:
            # Support 2-column (file phase) or 3-column (file phase exposure)
            raw = np.genfromtxt(list_path, dtype=None, encoding='utf-8')
            if raw.ndim == 0:
                raw = raw.reshape(1)
            files  = np.array([r[0] for r in raw])
            phases = np.array([float(r[1]) for r in raw])
            inputs = {'files': files, 'phase': phases}
        except IOError:
            raise IOError(f"Phase file '{list_path}' not found. "
                          "Check base_dir and list.")

        n_spec = len(inputs['files'])
        _info(f"Found {n_spec} spectra in {list_path}")

        # Read first spectrum to get wavelength grid
        w1st = np.loadtxt(
            os.path.join(self.base_dir, inputs['files'][0]),
            unpack=True)
        wo = w1st[0]

        if self.nbins is None:
            ph_arr = inputs['phase']
            self.nbins = int(1.5 / np.abs(ph_arr[2] - ph_arr[1]))
        _info(f"Phase bins: {self.nbins}")

        wave, flux = [], []
        print()
        for z in range(n_spec):
            fname = inputs['files'][z]
            ph    = inputs['phase'][z]
            w, f  = np.loadtxt(os.path.join(self.base_dir, fname), unpack=True)
            if z == 0:
                wave.append(w); flux.append(f)
            else:
                wave.append(wo); flux.append(np.interp(wo, w, f))
            bar = _progress_bar(z + 1, n_spec,
                                 label=f"{fname}  φ={ph:.4f}  npix={w.size}")
            print(f"\r{bar}", end='', flush=True)
        print()

        _ok(f"Loaded {n_spec} spectra, {wo.size} pixels each")

        # The MEM works on the individual spectra at their own phases; nbins
        # only sets the phase grid of the trail plot (rebin_trail in Dopin).
        # The old per-bin loop here indexed spectra by bin number (IndexError
        # when nbins > number of spectra) and its result was discarded anyway.
        self.wave        = wave
        self.flux        = flux
        self.input_files = np.array(inputs['files'])
        self.input_phase = np.array(inputs['phase'], dtype=float)
        self.pha         = self.input_phase
        self.trsp        = None   # set by Dopin (continuum-normalised spectra)


    # ──────────────────────────────────────────────────────────────────────────
    # Dopin
    # ──────────────────────────────────────────────────────────────────────────

    def Dopin(self, poly_degree: int = 2,
              continnum_band=False,
              plot_median: bool = False,
              xlim=None,
              two_orbits: bool = True,
              vel_space: bool = True,
              verbose: bool = False):
        """
        Continuum-normalise spectra and compute trail spectrum.

        Parameters
        ----------
        poly_degree     : degree of polynomial fit to continuum
        continnum_band  : [x1,x2,x3,x4] wavelength limits for two continuum
                          bands. If False, opens an interactive plot.
        plot_median     : show median spectrum above trail
        xlim            : wavelength limits for trail display
        two_orbits      : show two orbital cycles in trail
        vel_space       : x-axis in velocity (else wavelength)
        """
        _header("Dopin — continuum normalisation")

        lam = self.lam0
        cl  = 2.997e5

        if lam < min(self.wave[0]) or lam > max(self.wave[0]):
            raise ValueError(
                f"Input wavelength {lam:.2f} Å is outside data range "
                f"[{min(self.wave[0]):.2f}, {max(self.wave[0]):.2f}] Å")

        # ── Average spectrum + continuum fit ───────────────────────────────
        fig = plt.figure(num="Average Spectrum", figsize=(8, 6))
        plt.clf()
        ax  = fig.add_subplot(211)
        avgspec = np.sum(self.flux, axis=0)
        ax.plot(self.wave[0], avgspec / len(self.flux), color='steelblue',
                lw=1.2, label='Average spectrum')
        ax.set_xlabel(r'Wavelength / Å')
        ax.set_ylabel('Flux')
        ax.set_title(f"{self.object}  —  average spectrum")
        plt.tight_layout()
        plt.draw()

        if not continnum_band:
            _info("Interactive: click 4 points to define continuum bands")
            xor = []
            for i in range(4):
                xx = plt.ginput(1, timeout=-1)
                xor.append(xx[0][0])
                ax.axvline(x=xx[0][0], ls='--', color='tomato', alpha=0.7)
                plt.draw()
        else:
            xor = list(continnum_band)
            for i, x in enumerate(xor):
                lbl = 'Cont. bands' if i == 0 else ''
                ax.axvline(x=x, ls='--', color='tomato', alpha=0.7, label=lbl)
            plt.draw()

        lop  = ((self.wave[0] > xor[0]) & (self.wave[0] < xor[1])) | \
               ((self.wave[0] > xor[2]) & (self.wave[0] < xor[3]))
        yor  = avgspec[lop] / len(self.flux)
        z    = np.polyfit(self.wave[0][lop], yor, poly_degree)
        pz   = np.poly1d(z)
        cont = pz(self.wave[0])

        ax.plot(self.wave[0], cont, 'r--', lw=1.5, label='Continuum fit')
        ax.set_xlim(xor[0] - 10, xor[3] + 10)
        ax.set_ylim(yor.min() * 0.8, avgspec.max() / len(self.flux) * 1.1)
        ax.legend(fontsize=11)

        vell = ((self.wave[0] / self.lam0)**2 - 1) * cl / \
               (1 + (self.wave[0] / self.lam0)**2)

        ax2 = fig.add_subplot(212)
        ax2.plot(vell, avgspec / len(self.flux) - cont, 'k', lw=1.1,
                 label='Bkg-subtracted')
        ax2.axhline(y=0, ls='--', color='grey', alpha=0.5)
        hlim = self.delw / self.lam0 * cl
        ax2.axvline(x=-hlim, color='darkorange', lw=1.5, label='Map limits')
        ax2.axvline(x= hlim, color='darkorange', lw=1.5)
        ax2.set_xlim(-hlim * 1.5, hlim * 1.5)
        qq = np.abs(vell) < hlim * 1.5
        ymax = np.max(avgspec[qq] / len(self.flux) - cont[qq])
        ax2.set_ylim(-0.05 * ymax, ymax * 1.1)
        ax2.set_xlabel('Velocity / km s$^{-1}$')
        ax2.set_ylabel('Flux')
        ax2.legend(fontsize=11)
        plt.tight_layout()
        plt.savefig('Average_Spec.png', dpi=150, bbox_inches='tight')
        plt.draw()

        _ok("Continuum plot saved → Average_Spec.png")

        # ── Individual normalisation ───────────────────────────────────────
        nufac = (1.0 + self.gama / 2.998e5) * np.sqrt(1.0 - (self.gama / 2.998e5)**2)
        lop2  = (self.wave[0] / nufac > self.lam0 - self.delw) & \
                (self.wave[0] / nufac < self.lam0 + self.delw)
        self.normalised_wave = self.wave[0][lop2] / nufac

        vell_t = ((self.normalised_wave / self.lam0)**2 - 1.0) * cl / \
                 (1.0 + (self.normalised_wave / self.lam0)**2)
        self.vell = np.linspace(vell_t[0], vell_t[-1], vell_t.size)
        self.normalised_flux = np.zeros((len(self.flux), lop2.sum()))

        print()
        for ct, flu in enumerate(self.flux):
            polmask = ((self.wave[0] / nufac > xor[0]) & (self.wave[0] / nufac < xor[1])) | \
                      ((self.wave[0] / nufac > xor[2]) & (self.wave[0] / nufac < xor[3]))
            z2 = np.polyfit(self.wave[0][polmask] / nufac, flu[polmask], 3)
            pz2 = np.poly1d(z2)
            lf  = pz2(self.normalised_wave)
            nf  = flu[lop2] - lf
            self.normalised_flux[ct] = np.interp(self.vell, vell_t, nf)
            print(f"\r{_progress_bar(ct+1, len(self.flux), label=f'spectrum {ct+1}/{len(self.flux)}')}",
                  end='', flush=True)
        print()

        _ok(f"Velocity range: {self.vell.min():.1f} → {self.vell.max():.1f} km/s  "
            f"({self.vell.size} pixels)")

        self.pha  = self.input_phase
        self.trsp = self.normalised_flux

        # ── Trail spectrum plot ────────────────────────────────────────────
        trail, temp_phase = rebin_trail(self.vell, self.normalised_flux,
                                        self.input_phase, self.nbins,
                                        self.delta_phase)

        if xlim is None:
            rr = np.ones(self.normalised_wave.size, dtype=bool)
        else:
            rr = (self.normalised_wave > xlim[0]) & (self.normalised_wave < xlim[1])
        waver = self.vell if vel_space else self.normalised_wave[rr]

        fig2 = plt.figure('Trail Spectrum', figsize=(7, 9))
        plt.clf()

        rows = 6
        si   = 0 if plot_median else 2
        lo   = 2 if plot_median else 0

        if plot_median:
            ax1 = plt.subplot2grid((rows, 1), (0, 0), rowspan=2)
            ax1.plot(waver, np.nanmedian(self.normalised_flux, axis=0),
                     color='#8e44ad', lw=1.3, label='Median')
            ax1.axhline(y=0, ls='--', color='tomato', alpha=0.6)
            ax1.set_xticklabels([])
            ax1.legend(fontsize=11)
            ax1.minorticks_on()

        ax2 = plt.subplot2grid((rows, 1), (lo, 0), rowspan=4 + si)
        ax2.minorticks_on()

        if vel_space:
            x1 = (min(waver) - self.lam0) / self.lam0 * 2.998e5 if not vel_space else min(waver)
            x2 = (max(waver) - self.lam0) / self.lam0 * 2.998e5 if not vel_space else max(waver)
        else:
            x1, x2 = min(waver), max(waver)

        lim_two = 2 if two_orbits else 1
        ax2.imshow(trail.T, interpolation='nearest', cmap='binary',
                   aspect='auto', origin='lower',
                   extent=(min(waver), max(waver),
                           temp_phase[0], temp_phase[-1] + 1.0 / self.nbins))
        ax2.set_xlim(min(waver), max(waver))
        ax2.set_ylim(temp_phase[0], lim_two + 1.0 / self.nbins / 2.0)
        ax2.axvline(x=0 if vel_space else self.lam0, ls='--', color='darkorange', alpha=0.8)
        ax2.set_xlabel('Velocity / km s$^{-1}$' if vel_space else r'Wavelength / Å')
        ax2.set_ylabel('Orbital phase')
        ax2.set_title(f"{self.object}  —  trail spectrum")

        plt.tight_layout(h_pad=0)
        plt.savefig('Trail.png', dpi=150, bbox_inches='tight')
        plt.draw()
        _ok("Trail spectrum saved → Trail.png")


    # ──────────────────────────────────────────────────────────────────────────
    # Syncdop  —  pure Python MEM
    # ──────────────────────────────────────────────────────────────────────────

    def Syncdop(self, nv: Optional[int] = None):
        """
        Run Maximum Entropy Method Doppler tomography (pure Python).

        Parameters
        ----------
        nv : int, optional
            Velocity map grid size (nv × nv).  If None, computed from
            oversampling parameter self.overs.
        """
        _header("Syncdop — MEM Doppler tomography  (native Python)")

        pi2 = 2.0 * np.pi

        nph = self.pha.size
        nvp = self.vell.size

        # --- determine map size ---
        if nv is None:
            nv0 = int(self.overs * nvp)
            nv  = max(nv0, int(min(1.5 * nv0, nph / 3.0)))
            if nv % 2 == 1:
                nv += 1
        _info(f"Map grid:    {nv} × {nv} pixels")
        _info(f"Input data:  {nph} phases × {nvp} velocity channels")

        # Velocity axis in cm/s (MEM engine expects cm/s)
        vp_cgs = self.vell * 1e5

        # Phase bins in radians
        pha_rad  = self.pha * pi2
        dpha_rad = np.ones(nph) * self.delta_phase * pi2

        # Input data: (nph, nvp)
        dat = self.trsp.copy()

        # Progress output
        _alpha_count = [0]
        def progress(ial, nit, del_, rr):
            _alpha_count[0] = ial

        # Estimate sparse projection-matrix memory: each map pixel touches
        # ~(map pixel / data channel width + 2) channels per phase;
        # CSR stores 8-byte values + 4-byte column indices.
        nol_est = int(np.pi * (nv / 2)**2)
        dvi     = abs(self.vell[-1] - self.vell[0]) / max(nvp - 1, 1)
        dvm     = 2.0 * np.max(np.abs(self.vell)) / max(nv - 1, 1)
        nnz_est = nol_est * nph * (dvm / max(dvi, 1e-30) + 2.0)
        mem_mb  = max(1, int(12 * nnz_est / 1e6))
        _info(f"Projection matrix: sparse, ~{mem_mb} MB")

        t0 = time.time()
        print()
        print(f"  {'Alpha step':>10}  {'α':>10}  {'nit':>6}  {'χ':>8}  {'target':>8}")
        print(f"  {'─'*10}  {'─'*10}  {'─'*6}  {'─'*8}  {'─'*8}")

        # Custom per-step callback for pretty output
        _prev_ial = [0]
        def _alpha_cb(ial, nit, del_, rr):
            if ial != _prev_ial[0]:
                _prev_ial[0] = ial

        self._mem_result = doppler_mem(
            pha   = pha_rad,
            dpha  = dpha_rad,
            vp    = vp_cgs,
            dat   = dat,
            nv    = nv,
            ns    = self.ns,
            ac    = self.ac,
            nim   = self.nim,
            al0   = self.al0,
            alf   = self.alf,
            nal   = self.nal,
            clim  = self.clim,
            ih    = self.ih,
            iw    = self.iw,
            pb0   = self.pb0,
            pb1   = self.pb1,
            norm  = self.norm,
            wid   = self.wid,
            af    = self.af,
            verbose = self.verbose,
            progress_fn = _alpha_cb,
        )

        t_elapsed = time.time() - t0
        print()

        r = self._mem_result
        if r.converged:
            _ok(f"Converged!  χ = {r.chi2:.4f}  (target ≤ {self.clim:.2f})  "
                f"[{t_elapsed:.1f} s]")
        else:
            _warn(f"Not fully converged.  χ = {r.chi2:.4f}  (target ≤ {self.clim:.2f})  "
                  f"[{t_elapsed:.1f} s]")

        _info(f"Final α = {r.alpha:.2e}   iterations = {r.n_iter}")
        _info(f"Entropy S = {r.entropy:.4e}")


    # ──────────────────────────────────────────────────────────────────────────
    # Dopmap
    # ──────────────────────────────────────────────────────────────────────────

    def Dopmap(self, cmaps=None, limits=None,
               colorbar: bool = True, negative: bool = False,
               remove_mean: bool = False, corrx: int = 0, corry: int = 0,
               smooth: bool = False,
               stretch: str = 'linear',
               percentile_lo: float = 2.0,
               percentile_hi: float = 99.5):
        """
        Display the Doppler velocity map.

        Parameters
        ----------
        cmaps : colormap, optional
            Default is ``plt.cm.magma_r`` (pale background, dark emission).
            Pass any matplotlib colormap, e.g. ``plt.cm.afmhot`` for the
            dark-background look.
        limits : [vmin, vmax], optional
            Colour-scale limits on the *stretched* image (0–1 range after
            stretch).  If None, set from ``percentile_lo``/``percentile_hi``.
        stretch : str
            Intensity stretch before display: ``'linear'`` (default),
            ``'sqrt'`` or ``'log'``.  ``'sqrt'``/``'log'`` compress the bright
            disc core so fainter features (donor, stream hot-spot) stand out.
        percentile_lo, percentile_hi : float
            Percentile clip applied to the stretched image to set the colour
            limits (ignored when ``limits`` is supplied).  Default 2 / 99.5.
        smooth : bool
            Apply a mild Gaussian interpolation when displaying.
        colorbar : bool
            Show a colour bar.

        Returns
        -------
        (figure, image_2d_array)   —  the matplotlib figure and the
            normalised/stretched 2-D array that was plotted.
        """
        if self._mem_result is None:
            raise RuntimeError("Run Syncdop() first.")

        _header("Dopmap — velocity map")

        # Default colormap: magma_r (pale background, dark emission)
        if cmaps is None:
            cmaps = plt.cm.magma_r

        r    = self._mem_result
        data = r.image.copy()
        va   = r.va           # km/s

        # ── mask outside the |v| < va circle ─────────────────────────────
        # (not data == 0: the baseline map is subtracted inside |v| < va
        #  only, so the thin rim va <= |v| < va + dv/10 would show as a ring)
        nvm  = data.shape[0]
        dvm  = 2.0 * va / (nvm - 1)
        vgrid = dvm * (np.arange(1, nvm + 1) - nvm / 2.0 - 0.5)
        outside = np.hypot(vgrid[:, None], vgrid[None, :]) >= va
        data[outside] = np.nan

        # ── normalise: map [low-emission … high-emission] → [0 … 1] ─────────
        # We clip at percentile_lo to suppress the very lowest baseline noise
        # and at percentile_hi to avoid a single bright pixel dominating.
        dmin = np.nanpercentile(data, percentile_lo)
        dmax = np.nanpercentile(data, percentile_hi)
        data_norm = (data - dmin) / max(dmax - dmin, 1e-30)
        data_norm = np.clip(data_norm, 0.0, 1.0)

        # ── intensity stretch ─────────────────────────────────────────────
        # sqrt pulls up faint features (disc ring, donor) relative to the
        # bright disc core.  log is more aggressive.
        if stretch == 'sqrt':
            stretched = np.sqrt(data_norm)
        elif stretch == 'log':
            stretched = np.log1p(data_norm * 99.0) / np.log(100.0)
        else:  # linear
            stretched = data_norm

        # ── colour limits on the stretched image ──────────────────────────
        finite = stretched[np.isfinite(stretched)]
        if limits is None:
            # Keep full range after stretch
            vmin = float(np.nanmin(finite))
            vmax = float(np.nanmax(finite))
        else:
            vmin, vmax = limits

        _info(f"Velocity range:  {-va:.0f} → {va:.0f} km/s  "
              f"({data.shape[0]}×{data.shape[1]} pixels)")
        _info(f"Stretch:         {stretch}")
        _info(f"Colour limits:   [{vmin:.3f}, {vmax:.3f}]  "
              f"(data clipped at p{percentile_lo:.0f}–p{percentile_hi:.1f})")

        # ── axis extents ──────────────────────────────────────────────────
        dv    = 2.0 * va / (data.shape[0] - 1)
        vpmin = -va - 0.5 * dv
        vpmax =  va + 0.5 * dv

        interp = 'gaussian' if smooth else 'nearest'

        # ── figure ────────────────────────────────────────────────────────
        fig = plt.figure(num='Doppler Map', figsize=(7, 7))
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_facecolor(cmaps(1.0 if negative else 0.0))
        ax.minorticks_on()

        plot_data = -stretched if negative else stretched
        _vmin = -vmax if negative else vmin
        _vmax = -vmin if negative else vmax

        img = ax.imshow(plot_data, interpolation=interp,
                        cmap=cmaps, aspect='equal', origin='lower',
                        extent=(vpmin, vpmax, vpmin, vpmax),
                        vmin=_vmin, vmax=_vmax)

        ax.set_xlim(-va, va)
        ax.set_ylim(-va, va)
        ax.set_xlabel('$V_x$ / km s$^{-1}$', fontsize=13)
        ax.set_ylabel('$V_y$ / km s$^{-1}$', fontsize=13)
        ax.set_title(f"{self.object}  —  Doppler map", fontsize=13, pad=10)

        if colorbar:
            cbar = fig.colorbar(img, fraction=0.046, pad=0.04)
            cbar.set_label(f'Normalised flux  ({stretch} stretch)', fontsize=11)
            # Annotate the colour bar with physically meaningful percentiles
            cbar.ax.axhline(y=(np.nanpercentile(finite, 50) - _vmin) /
                              max(_vmax - _vmin, 1e-9),
                            color='grey', lw=0.8, ls='--', alpha=0.8)
        else:
            cbar = None

        plt.tight_layout()
        plt.savefig('Doppler_Map.png', dpi=150, bbox_inches='tight')
        plt.draw()
        _ok("Doppler map saved → Doppler_Map.png")

        return fig, stretched


    # ──────────────────────────────────────────────────────────────────────────
    # Reco
    # ──────────────────────────────────────────────────────────────────────────

    def Reco(self, cmaps=plt.cm.binary, limits=None, colorbar: bool = True):
        """
        Plot original and reconstructed trail spectra side by side.

        Returns
        -------
        (fig, dm_input, dm_reconstructed)
        """
        if self._mem_result is None:
            raise RuntimeError("Run Syncdop() first.")

        _header("Reco — reconstructed trail spectra")

        r  = self._mem_result
        vp = r.vp     # km/s

        trail_dm,  phase = rebin_trail(vp, r.dm,      r.phases, self.nbins,
                                        self.delta_phase)
        trail_dmr, _     = rebin_trail(vp, r.dm_reco, r.phases, self.nbins,
                                        self.delta_phase)

        mx_dm  = np.nanmax(trail_dm)  or 1.0
        mx_dmr = np.nanmax(trail_dmr) or 1.0

        if limits is None:
            obs = trail_dm / mx_dm
            limits = [np.nanpercentile(obs, 2.0), np.nanpercentile(obs, 99.5)]

        # pixel-centre velocities -> image edges
        dvel   = (vp[-1] - vp[0]) / (vp.size - 1)
        x1, x2 = vp.min() - 0.5 * dvel, vp.max() + 0.5 * dvel

        # constrained layout leaves room for the shared colour bar
        # (tight_layout ignores colour-bar axes, which made it overlap)
        fig = plt.figure(num='Reconstruction', figsize=(11, 7),
                         layout='constrained')
        fig.clf()
        ax1, ax2 = fig.subplots(1, 2, sharey=True)

        # origin='lower' so phase increases upwards, matching the extent
        # and the Dopin trail plot (origin='upper' drew it flipped)
        common = dict(interpolation='nearest', cmap=cmaps, aspect='auto',
                      origin='lower',
                      extent=(x1, x2, phase[0], phase[-1] + 1.0 / self.nbins),
                      vmin=limits[0], vmax=limits[1])

        ax1.imshow(trail_dm.T / mx_dm, **common)
        ax1.set_xlabel('Velocity / km s$^{-1}$', fontsize=12)
        ax1.set_ylabel('Orbital phase',           fontsize=12)
        ax1.set_title('Observed trail',           fontsize=12)
        ax1.minorticks_on()

        ax2.imshow(trail_dmr.T / mx_dmr, **common)
        ax2.set_xlabel('Velocity / km s$^{-1}$', fontsize=12)
        ax2.set_title('Reconstructed trail',      fontsize=12)
        ax2.minorticks_on()
        ax2.yaxis.set_tick_params(labelleft=False)

        if colorbar:
            fig.colorbar(ax2.get_images()[0], ax=[ax1, ax2],
                         shrink=0.9, pad=0.01, label='Normalised flux')

        fig.suptitle(f"{self.object}  —  trail comparison", fontsize=13)
        fig.savefig('Reconstruction.png', dpi=150, bbox_inches='tight')
        plt.draw()
        _ok("Reconstruction plot saved → Reconstruction.png")

        return fig, r.dm, r.dm_reco


# ─────────────────────────────────────────────────────────────────────────────
# stream  (pure Python — ported from original, unchanged physics)
# ─────────────────────────────────────────────────────────────────────────────

def stream(q, k1, porb, m1, inc,
           colors='black', both_lobes: bool = False,
           title: bool = True, label=None,
           tidal_radius: bool = False,
           resonance_32: bool = False,
           resonance_21: bool = False,
           lw: float = 1.2,
           alpha: float = 0.85):
    """
    Overplot binary system geometry on a Doppler map.

    Draws the ballistic stream from L1, the Keplerian stream, the Roche
    lobe of the donor, and optionally the tidal truncation radius and
    the 3:2 / 2:1 orbital resonance circles.

    Parameters
    ----------
    q    : float
        Mass ratio M2/M1.
    k1   : float
        Radial-velocity semi-amplitude of the primary (km/s).
    porb : float
        Orbital period (days).
    m1   : float
        Primary mass (solar masses).
    inc  : float
        Orbital inclination (degrees).
    colors : str or color
        Colour for stream lines and markers.  Default ``'black'`` suits
        the pale ``magma_r`` default Doppler map; use ``'white'`` on dark
        colormaps such as ``afmhot``.
    both_lobes : bool
        If True, plot the full donor Roche lobe contour (both hemispheres).
    title : bool
        Add system parameters as the axes title.
    label : str, optional
        Text annotation placed at the lower-left of the axes.
    tidal_radius : bool
        Overplot the tidal truncation radius of the accretion disc
        (Paczynski 1977 approximation: r_tid ≈ 0.60 / (1+q) * a).
        Plotted as a dashed circle.
    resonance_32 : bool
        Overplot the 3:2 mean-motion resonance circle (Whitehurst & King
        1991; relevant for SU UMa-type superhumps).
    resonance_21 : bool
        Overplot the 2:1 resonance circle.
    lw : float
        Line width for all curves.
    alpha : float
        Line opacity.

    Notes
    -----
    All radii are converted to km/s via the Keplerian velocity at that
    radius:  v_Kep = sqrt(G M1 / r) * sin(i).
    """
    xl, yl, xi, yi, wout, wkout = stream_calculate(q)

    azms = -70.0
    az   = np.arctan2(yi, xi)
    for i in range(len(az)):
        if xi[i] < 0.0:
            az[i] += np.pi
    az = np.degrees(az)

    cut = 0
    for j in range(az.size):
        cut = j
        if az[j] < azms:
            break

    vxi  = np.real(wout[:cut])
    vyi  = np.imag(wout[:cut])
    vkxi = np.real(wkout[:cut])
    vkyi = np.imag(wkout[:cut])

    # ── physical scale ────────────────────────────────────────────────────
    porb_s = porb * 24.0 * 3600.0
    omega  = 2.0 * np.pi / porb_s
    gg     = 6.667e-8
    msun   = 1.989e33
    cm_    = q / (1.0 + q)
    a      = (gg * m1 * msun * (1.0 + q))**(1.0 / 3.0) / omega**(2.0 / 3.0)
    vfs    = 1e5           # cm/s per km/s
    vs     = omega * a / vfs   # orbital velocity scale (km/s)

    inc_r = np.radians(inc)
    si    = np.sin(inc_r)

    xl_c = xl - cm_
    vx   = vxi  * si * vs
    vy   = vyi  * si * vs
    vkx  = vkxi * si * vs
    vky  = vkyi * si * vs
    xl_v = xl_c * vs * si
    yl_v = yl   * vs * si

    ax = plt.gca()

    # ── stream curves ─────────────────────────────────────────────────────
    ax.plot(vx,  vy,  color=colors, lw=lw, alpha=alpha, label='Ballistic stream')
    ax.plot(vkx, vky, color=colors, lw=lw, alpha=alpha, ls='--',
            label='Keplerian stream')

    # ── Roche lobe of donor ───────────────────────────────────────────────
    n4 = yl_v.size // 4
    ax.plot(yl_v[n4:3*n4], xl_v[n4:3*n4], color=colors, lw=lw, alpha=alpha,
            label='Roche lobe')
    if both_lobes:
        yl_top = np.concatenate((yl_v[3*n4:], yl_v[:n4]))
        xl_top = np.concatenate((xl_v[3*n4:], xl_v[:n4]))
        ax.plot(yl_top, xl_top, color=colors, lw=lw, alpha=alpha * 0.5, ls=':')

    # ── system markers ────────────────────────────────────────────────────
    vy1 = cm_ * vs * si                   # K1 (WD velocity)
    ax.plot(0.0, 0.0,   'x', ms=8, c=colors, alpha=0.4, mew=1.5)   # CoM
    ax.plot(0.0, -vy1,  '+', ms=10, c=colors, alpha=0.8, mew=1.5,
            label=f'WD  ($K_1$ = {vy1:.0f} km/s)')
    ax.plot(0.0, (1.0 - cm_) * vs * si,
            '+', ms=10, c=colors, alpha=0.8, mew=1.5, label='Donor CoM')

    # ── tidal truncation radius ────────────────────────────────────────────
    # Paczynski 1977 / Papaloizou & Pringle 1977:
    #   r_tid / a  ≈  0.60 / (1 + q)        (commonly used approximation)
    # Keplerian velocity at r_tid:
    #   v_tid = sqrt(G M1 / r_tid) * sin(i) = vs * sqrt(a / r_tid) * sin(i)
    if tidal_radius:
        r_tid_a = 0.60 / (1.0 + q)        # r_tid in units of a
        v_tid   = vs * si / np.sqrt(r_tid_a)   # Keplerian v at r_tid (km/s)
        phi     = np.linspace(0, 2 * np.pi, 300)
        ax.plot(v_tid * np.cos(phi), v_tid * np.sin(phi) - vy1,
                color=colors, lw=lw, ls=':', alpha=alpha,
                label=f'Tidal radius  ({v_tid:.0f} km/s)')
        _info(f"Tidal radius:    r/a = {r_tid_a:.3f}   "
              f"v_Kep = {v_tid:.1f} km/s  (at WD)")

    # ── 3:2 resonance ─────────────────────────────────────────────────────
    # r_{j:k} / a = [(j-k)/j]^{2/3} / (1+q)^{1/3}
    # (Warner 1995, eq. 2.61)
    if resonance_32:
        j, k = 3, 2
        r_a  = ((j - k) / j)**(2.0 / 3.0) / (1.0 + q)**(1.0 / 3.0)
        v_r  = vs * si / np.sqrt(r_a)
        phi  = np.linspace(0, 2 * np.pi, 300)
        ax.plot(v_r * np.cos(phi), v_r * np.sin(phi) - vy1,
                color=colors, lw=lw * 1.3, ls='-.', alpha=alpha,
                label=f'3:2 resonance  ({v_r:.0f} km/s)')
        _info(f"3:2 resonance:   r/a = {r_a:.3f}   v_Kep = {v_r:.1f} km/s")

    # ── 2:1 resonance ─────────────────────────────────────────────────────
    if resonance_21:
        j, k = 2, 1
        r_a  = ((j - k) / j)**(2.0 / 3.0) / (1.0 + q)**(1.0 / 3.0)
        v_r  = vs * si / np.sqrt(r_a)
        phi  = np.linspace(0, 2 * np.pi, 300)
        ax.plot(v_r * np.cos(phi), v_r * np.sin(phi) - vy1,
                color=colors, lw=lw * 1.3, ls=(0, (5, 2, 1, 2)), alpha=alpha,
                label=f'2:1 resonance  ({v_r:.0f} km/s)')
        _info(f"2:1 resonance:   r/a = {r_a:.3f}   v_Kep = {v_r:.1f} km/s")

    # ── title & annotation ────────────────────────────────────────────────
    if title:
        ax.set_title(
            rf"$i$={inc:.1f}°, $M_1$={m1} M$_\odot$, "
            rf"$q$={q}, $P_{{orb}}$={porb*24:.2f} hr",
            fontsize=11
        )

    if label is not None:
        ax.text(0.12, 0.1, label, ha='center', va='center',
                transform=ax.transAxes, color=colors, fontsize=11)

    plt.tight_layout()


def stream_calculate(qm, ni: int = 100, nj: int = 100):
    """Calculate Roche lobe and ballistic stream trajectory from L1."""
    nmax  = 10000
    xout  = np.zeros(nmax)
    yout  = np.zeros(nmax)
    rout  = np.zeros(nmax)
    wout  = np.zeros(nmax, dtype=complex)
    wkout = np.zeros(nmax, dtype=complex)

    if abs(qm - 1.0) < 1e-4:
        qm = 1e-4
    if qm <= 0.0:
        raise ValueError("Mass ratio q must be > 0.")

    rl1 = rlq1(qm)
    x, y = lobes(qm, rl1, ni, nj)

    cm_  = qm / (1.0 + qm)
    z1   = complex(-cm_, 0.0)
    z2   = complex(1.0 - cm_, 0.0)
    wm1  = complex(0.0, -cm_).conjugate()

    eps  = 1e-3
    z    = complex(rl1 - cm_ - eps, 0.0)
    w    = 0.0 + 0.0j

    t = 0.0; dt = 1e-4; it = 0; ph = 0.0; phmax = 6.0
    r = 1.0; ist = 0; isa = 0

    while it < nmax and ph < phmax:
        dz, dw = _intrk(z, w, dt, z1, z2, qm)
        z += dz; w += dw; t += dt

        if abs(dz) / (abs(z) + 1e-30) > 0.02:
            dt /= 2.0
        if abs(dz) / (abs(z) + 1e-30) < 0.005:
            dt *= 2.0

        dph = -np.imag(z * np.conj(z - dz)) / (abs(z) * abs(z - dz) + 1e-30)
        ph += dph

        wi   = w + 1j * z
        rold = r
        r    = abs(z - z1)

        if ist == 0 and rold < r:
            ist = 1

        vk = 1.0 / np.sqrt(r * (1.0 + qm) + 1e-30)
        no = np.conj(z - z1) / (r + 1e-30)
        wk = -vk * no * 1j + wm1

        xout[it]  = np.real(z) + cm_
        yout[it]  = -np.imag(z)
        rout[it]  = np.sqrt(xout[it]**2 + yout[it]**2)
        wout[it]  = wi
        wkout[it] = np.conj(wk)
        it += 1

    return x, y, xout[:it], yout[:it], wout[:it], wkout[:it]


def rlq1(q: float) -> float:
    """L1 Lagrangian point radius (Spruit)."""
    if abs(1.0 - q) < 1e-4:
        return 0.5
    rl = 0.0
    rn = 1.0 - q
    for _ in range(100):
        rl = rn
        f  = q / (1.0 - rl)**2 - 1.0 / rl**2 + (1.0 + q) * rl - q
        fa = 2.0 * q / (1.0 - rl)**3 + 2.0 / rl**3 + (1.0 + q)
        rn = rl - f / fa
        if abs(rl / rn - 1.0) < 1e-6:
            break
    return rn


def lobes(q: float, rs: float, ni: int, nj: int):
    """Compute Roche lobe contour."""
    r  = np.zeros((ni, nj))
    ch = np.zeros(ni)
    ps = np.zeros(nj)

    r, ch, ps = _surface(q, rs, ni, nj, r, ch, ps)
    x  = np.array([1.0 - r[i, 0] * np.cos(ch[i]) for i in range(ni)])
    y  = np.array([-r[i, 0] * np.sin(ch[i])        for i in range(ni)])

    r, ch, ps = _surface(1.0 / q, 1.0 - rs, ni, nj, r, ch, ps)
    x2 = np.array([r[i, 0] * np.cos(ch[i]) for i in range(ni)])
    y2 = np.array([r[i, 0] * np.sin(ch[i]) for i in range(ni)])

    xt = np.concatenate((x2[::-1], x, x[::-1], x2))
    yt = np.concatenate((y2[::-1], -y, y[::-1], -y2))
    return xt, yt


def _pot(q, x, y, z):
    r  = np.sqrt(x*x + y*y + z*z) + 1e-30
    rh = np.sqrt(x*x + y*y)
    st = rh / r
    cf = (x / rh) if rh > 0 else 1.0
    r2 = 1.0 / (1.0 + q)
    r1 = np.sqrt(1.0 + r**2 - 2.0 * r * cf * st)
    p  = -1.0/r - 1.0/q/r1 - 0.5*(1.0/q+1.0)*(r2**2+(r*st)**2 - 2.0*r2*r*cf*st)
    pr = (1.0/r**2 + 1.0/q/(r1**3)*(r - cf*st)
          - 0.5*(1.0/q+1)*2.0*(r*st*st - r2*cf*st))
    return p, pr


def _surface(q, rs, nc, nop, r, ch, ps):
    r  = np.zeros((100, 100))
    ch = np.zeros(nc)
    ps = np.zeros(nop)
    ch[0] = 0.0
    for i in range(1, nc):
        ch[i] = float(i - 1) * np.pi / (nc - 1.0)
    for j in range(1, nop):
        ps[j] = float(j - 1) * 2.0 * np.pi / nop
    rs1 = 1.0 - rs
    fs, _ = _pot(q, rs1, 0.0, 0.0)
    im = 20
    for i in range(nop):
        cp = np.cos(ps[i]); sp = np.sin(ps[i])
        rx = (1.0 - np.pi / nc) * rs1
        r[0, i] = rs1
        for k in range(1, nc):
            xc = np.cos(ch[k]); sc = np.sin(ch[k])
            yc = sc * cp;       zc = sc * sp
            j = 0; f = 1.0
            while j < im and abs(f - fs) > 1e-4 or j == 0:
                r1 = rx
                f, pr = _pot(q, r1 * xc, r1 * yc, r1 * zc)
                rx = r1 - (f - fs) / (pr + 1e-30)
                if rx > rs1:
                    rx = rs1
                j += 1
            r[k, i] = rx
    return r, ch, ps


def _eqmot(z, w, z1, z2, qm):
    zr1 = z - z1; zr2 = z - z2
    wp  = (-(qm * zr2 / (abs(zr2)**3 + 1e-30) +
              zr1 / (abs(zr1)**3 + 1e-30)) / (1.0 + qm)
            - 2j * w + z)
    return w, wp


def _intrk(z, w, dt, z1, z2, qm):
    zp, wp = _eqmot(z, w, z1, z2, qm)
    hz0, hw0 = zp * dt, wp * dt
    zp, wp = _eqmot(z + hz0/2, w + hw0/2, z1, z2, qm)
    hz1, hw1 = zp * dt, wp * dt
    zp, wp = _eqmot(z + hz1/2, w + hw1/2, z1, z2, qm)
    hz2, hw2 = zp * dt, wp * dt
    zp, wp = _eqmot(z + hz2, w + hw2, z1, z2, qm)
    hz3, hw3 = zp * dt, wp * dt
    return (hz0 + 2*hz1 + 2*hz2 + hz3) / 6.0, (hw0 + 2*hw1 + 2*hw2 + hw3) / 6.0
