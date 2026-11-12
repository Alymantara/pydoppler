"""
mem_engine.py  —  Pure NumPy/SciPy port of Spruit's dop.f MEM tomography engine.

Spruit H.C. (1998) "Doppler Imaging of Accretion Discs"

Algorithm overview
------------------
Given a time-series of spectra sampled at different orbital phases, we seek
the velocity map I(vx, vy) on a 2-D grid that, when projected through the
phase-dependent line-of-sight operator, best reproduces the observed spectra
while maximising image entropy (flattest consistent map).

The projection:
    d[phase, v_los] = sum_{vx,vy} P[phase,v_los; vx,vy] * I[vx,vy]

P is built once then held as a scipy sparse (CSR) matrix.
The MEM iteration follows Lucy (1994) RL with 2-D acceleration.
"""

from __future__ import annotations

import numpy as np
import scipy.ndimage as ndi
from scipy.sparse import csr_matrix
from dataclasses import dataclass, field
from typing import Optional
import warnings


# ---------------------------------------------------------------------------
# Public result container
# ---------------------------------------------------------------------------

@dataclass
class MEMResult:
    """Container for Doppler MEM output."""
    image: np.ndarray          # 2-D velocity map [nv, nv]
    va: float                  # half-width of velocity grid (km/s)
    vp: np.ndarray             # input velocity axis (km/s)
    dm: np.ndarray             # input data (phase, vp) after pre-processing
    dm_reco: np.ndarray        # reconstructed data
    phases: np.ndarray         # input phases (radians)
    alpha: float               # final regularisation alpha
    chi2: float                # final reduced chi²
    converged: bool
    n_iter: int
    entropy: float


# ---------------------------------------------------------------------------
# Utility: piecewise-quadratic weight function (Spruit's dist())
# ---------------------------------------------------------------------------

def _dist(x: np.ndarray) -> np.ndarray:
    """Spruit's piecewise-quadratic weight (vectorised)."""
    out = np.where(x < 0.0, 0.0,
          np.where(x < 0.5, 2.0 * x**2,
          np.where(x < 1.0, 1.0 - 2.0 * (1.0 - x)**2,
                             1.0)))
    return out


# ---------------------------------------------------------------------------
# 1.  Build projection matrix
# ---------------------------------------------------------------------------

def build_projection(pha: np.ndarray, dpha: np.ndarray,
                     vp: np.ndarray, nv: int) -> tuple[csr_matrix, np.ndarray, np.ndarray, float]:
    """
    Build the sparse projection matrix P  (shape: nol × nil).

    Convention (matching Spruit's Fortran):
        forward projection (map → data):  d = P.T @ image_vec   (nol→nil)
        back    projection (data → map):  b = P   @ data_vec    (nil→nol)

    P[ir, inp] = weight of map pixel ir on data point inp.
    P.T[inp, ir] is what projects image into data.

    Parameters
    ----------
    pha   : (nph,)  mid-exposure phases in radians
    dpha  : (nph,)  exposure lengths in radians
    vp    : (nvp,)  velocity axis of input spectra (cm/s)
    nv    : int     number of pixels on each side of velocity map

    Returns
    -------
    P     : CSR sparse (nol, nil)  — normalised projection matrix
    indx  : (nol,) x-pixel index for each map point inside the circle
    indy  : (nol,) y-pixel index
    va    : float   half-width of velocity grid (cm/s)
    """
    nph  = len(pha)
    nvp  = len(vp)
    nil  = nph * nvp

    va   = max(abs(vp[0]), abs(vp[-1]))
    dv   = 2.0 * va / (nv - 1)
    dvi  = (vp[-1] - vp[0]) / (nvp - 1)

    # Pre-compute cos/sin of phase bin edges
    phl  = pha - dpha / 2.0
    phu  = pha + dpha / 2.0
    cfl, sfl = np.cos(phl), np.sin(phl)
    cfu, sfu = np.cos(phu), np.sin(phu)

    # Velocity bin edges, 0-based: data column k spans [edges[k], edges[k+1]].
    # (Fortran prom uses 1-based vpl(iv) for bin iv -> data index
    #  ip+(iv-1)*nph; edges[k] == vpl(k+1).  An earlier port mixed the
    #  Fortran edge indices with a 0-based data index, which shifted every
    #  projected spectrum by +1 velocity channel.)
    edges = np.empty(nvp + 1)
    edges[1:nvp] = 0.5 * (vp[:-1] + vp[1:])
    edges[0]     = vp[0]  - (vp[1]  - vp[0])  / 2.0
    edges[nvp]   = vp[-1] + (vp[-1] - vp[-2]) / 2.0

    # Collect map pixels inside velocity circle (ix-major order, 1-based)
    fv2 = nv / 2.0
    ixg, iyg = np.meshgrid(np.arange(1, nv + 1), np.arange(1, nv + 1),
                           indexing='ij')
    vxg = dv * (ixg - fv2 - 0.5)
    vyg = dv * (iyg - fv2 - 0.5)
    inside = vxg**2 + vyg**2 < (va + dv / 10.0)**2

    indx = ixg[inside].astype(np.int32)
    indy = iyg[inside].astype(np.int32)
    nol  = len(indx)

    # ---- Vectorised over all (map pixel, phase) pairs: shape (nol, nph) ----
    vx = (dv * (indx - fv2 - 0.5))[:, None]
    vy = (dv * (indy - fv2 - 0.5))[:, None]

    # Line-of-sight projection: v_los = vx*sin(phi) - vy*cos(phi)
    # (Spruit 1998, eq. for the doppler mapping operator)
    vp0 = vx * sfl[None, :] - vy * cfl[None, :]
    vp1 = vx * sfu[None, :] - vy * cfu[None, :]

    swap = vp0 > vp1
    vpma = np.where(swap, vp0, vp1) + dv / 2.0
    vpmi = np.where(swap, vp1, vp0) - dv / 2.0

    # Data columns overlapping [vpmi, vpma]: first column whose upper edge
    # exceeds vpmi, last column whose lower edge is below vpma.  Columns
    # off the data grid are dropped (as Fortran's ki>=1, ka<=nvpm clamps do).
    ki = np.maximum(0, np.searchsorted(edges[1:], vpmi, side='right'))
    ka = np.minimum(nvp - 1,
                    np.searchsorted(edges[:-1], vpma, side='left') - 1)

    # Loop only over the (small) number of bins each pair can touch
    width = vpma - vpmi + 1e-30
    nspan = int(max(0, (ka - ki).max() + 1)) if ki.size else 0
    ip_idx = np.broadcast_to(np.arange(nph)[None, :], (nol, nph))
    ir_idx = np.broadcast_to(np.arange(nol)[:, None], (nol, nph))

    rows, cols, vals = [], [], []
    for j in range(nspan):
        iv = ki + j
        valid = iv <= ka
        if not valid.any():
            continue
        ivc = np.where(valid, iv, 0)          # safe index for masked entries
        x0 = np.clip((edges[ivc]     - vpmi) / width, 0.0, 1.0)
        x1 = np.clip((edges[ivc + 1] - vpmi) / width, 0.0, 1.0)
        w  = _dist(x1) - _dist(x0)
        keep = valid & (w > 0)
        rows.append(ir_idx[keep])
        cols.append(ip_idx[keep] + ivc[keep] * nph)   # phase-major data index
        vals.append(w[keep])

    rows = np.concatenate(rows) if rows else np.empty(0, dtype=np.int64)
    cols = np.concatenate(cols) if cols else np.empty(0, dtype=np.int64)
    vals = np.concatenate(vals) if vals else np.empty(0)

    # Normalise each map pixel's row to unit total weight
    s_ir = np.bincount(rows, weights=vals, minlength=nol)
    inv  = np.zeros(nol)
    pos  = s_ir > 0
    inv[pos] = 1.0 / s_ir[pos]
    vals = vals * inv[rows]

    # Each (pixel, data point) pair occurs once, so no duplicates to sum
    P = csr_matrix((vals, (rows, cols)), shape=(nol, nil))

    return P, indx, indy, va


# ---------------------------------------------------------------------------
# 2.  Gaussian smoothing of default image (chims/chimg)
# ---------------------------------------------------------------------------

_KERNEL_CACHE: dict = {}


def _spruit_kernel(ns: int) -> np.ndarray:
    """Spruit's chimg kernel: exp(-i^2/ns^2) for |i| <= int(1.5*ns), normalised."""
    k = _KERNEL_CACHE.get(ns)
    if k is None:
        ns1 = int(1.5 * ns)
        i = np.arange(-ns1, ns1 + 1, dtype=np.float64)
        k = np.exp(-i**2 / float(ns)**2)
        k /= k.sum()
        _KERNEL_CACHE[ns] = k
    return k


def _smooth_image(x_flat: np.ndarray, indx: np.ndarray, indy: np.ndarray,
                  nv: int, ns: int) -> np.ndarray:
    """
    Default-map smearing (Fortran chims/chimg): separable Gaussian
    exp(-i^2/ns^2) truncated at 1.5*ns, with fold-back boundaries
    (scipy 'reflect' == Fortran x2(i, 1-j)).
    """
    img = np.zeros((nv, nv))
    img[indx - 1, indy - 1] = x_flat
    if ns > 0:
        k = _spruit_kernel(ns)
        img = ndi.correlate1d(img, k, axis=1, mode='reflect')   # j-direction
        img = ndi.correlate1d(img, k, axis=0, mode='reflect')   # i-direction
    return img[indx - 1, indy - 1]


# ---------------------------------------------------------------------------
# 3.  Likelihood and entropy
# ---------------------------------------------------------------------------

def _likelihood(d_obs: np.ndarray, d_pred: np.ndarray,
                weights: np.ndarray, ih: int):
    """
    Likelihood H and gradient dH/dpsi (back-projected later).

    ih=0 : log-likelihood (Poisson / Cash)
    ih=1 : weighted χ²
    """
    eps = 1e-30
    if ih == 0:
        h = np.sum(weights * (d_obs * np.log(d_pred + eps) - d_pred))
        dph = weights * (d_obs / (d_pred + eps) - 1.0)
    else:
        diff = d_obs - d_pred
        h    = -np.sum(weights * diff**2) * len(d_obs)
        dph  =  2.0 * len(d_obs) * weights * diff
    return h, dph


def _entropy(ps: np.ndarray, chi: np.ndarray,
             indx: np.ndarray, indy: np.ndarray, nv: int,
             ns: int) -> tuple[float, np.ndarray]:
    """
    Entropy with floating default chi = smooth(ps)  (Fortran entr, nc=1):

        S       = -sum[ ps (ln(ps/chi) - 1) + chi ]
        dS/dps  = -ln(ps/chi) + smooth((ps - chi)/chi)

    The second gradient term comes from chi's own dependence on ps.
    """
    eps  = 1e-30
    tem  = np.log(ps + eps) - np.log(chi + eps)
    s    = -np.sum(ps * (tem - 1.0)) - np.sum(chi)
    dsa  = -tem + _smooth_image((ps - chi) / (chi + eps), indx, indy, nv, ns)
    return s, dsa


# ---------------------------------------------------------------------------
# 4.  Lucy 2-D acceleration  (lamun)
# ---------------------------------------------------------------------------

def _lucy_step(d_obs: np.ndarray, d_pred: np.ndarray,
               ps: np.ndarray, dps1: np.ndarray, dps2: np.ndarray,
               P: csr_matrix, chi: np.ndarray,
               dchi1: np.ndarray, dchi2: np.ndarray,
               weights: np.ndarray, al: float, ih: int):
    """
    Compute optimal λ, μ for Lucy's 2-D acceleration.
    Returns fla, fmu.
    """
    dph1 = P.T.dot(dps1)   # forward: image → data
    dph2 = P.T.dot(dps2)
    eps  = 1e-30

    if ih == 0:
        phw  = d_obs * weights
        dpp  = d_obs / (d_pred + eps) - 1.0
        xl   = dph1 / (d_pred + eps)
        xm   = dph2 / (d_pred + eps)
        hl   = np.sum(weights * dph1 * dpp)
        hm   = np.sum(weights * dph2 * dpp)
        hll  = -np.sum(xl**2 * phw)
        hmm  = -np.sum(xm**2 * phw)
        hlm  = -np.sum(xl * xm * phw)
    else:
        diff = d_obs - d_pred
        nil  = len(d_obs)
        hl   = 2.0 * nil * np.sum(weights * diff * dph1)
        hm   = 2.0 * nil * np.sum(weights * diff * dph2)
        hll  = -2.0 * nil * np.sum(weights * dph1**2)
        hmm  = -2.0 * nil * np.sum(weights * dph2**2)
        hlm  = -2.0 * nil * np.sum(weights * dph1 * dph2)

    if al > 0.0:
        sup  = np.sum(ps)
        pl   = dps1 / (ps + eps)
        pm   = dps2 / (ps + eps)
        tl   = dchi1 / (chi + eps)
        tm   = dchi2 / (chi + eps)
        pz   = ps - chi
        sl   = np.sum(-dps1 * np.log(ps / (chi + eps)) + pz * tl)
        sm   = np.sum(-dps2 * np.log(ps / (chi + eps)) + pz * tm)
        # Fortran lamun adds pz*(tl**2 - tll) etc.; with a single default
        # component (nc=1, epc=1) tll == tl**2, so those terms vanish.
        sll  = np.sum(-ps * (pl - tl)**2)
        smm  = np.sum(-ps * (pm - tm)**2)
        slm  = np.sum(-ps * (pl - tl) * (pm - tm))
    else:
        sl = sm = sll = smm = slm = 0.0

    ql  = hl  + al * sl
    qm  = hm  + al * sm
    qll = hll + al * sll
    qmm = hmm + al * smm
    qlm = hlm + al * slm

    dd = qmm * qll - qlm**2
    dr = abs(dd) / (qll**2 + qmm**2 + eps)
    if dr < 1e-3:
        fla, fmu = 0.0, -qm / (qmm + eps)
    else:
        fmu = (qlm * ql - qll * qm) / dd
        fla = (qlm * qm - qmm * ql) / dd
    return fla, fmu


# ---------------------------------------------------------------------------
# 5.  Inner RL-MEM iteration  (rlfldn)
# ---------------------------------------------------------------------------

def _rlfldn(d_obs: np.ndarray, P: csr_matrix, ps: np.ndarray,
            indx: np.ndarray, indy: np.ndarray, nv: int,
            al: float, ac: float, nim: int, ns: int, ih: int,
            weights: np.ndarray, flim: float = 1e5,
            callback=None) -> tuple[np.ndarray, int, float, float, int]:
    """
    Single RL-MEM iteration sequence at fixed alpha.

    Returns
    -------
    ps_out  : updated map
    nit     : number of iterations
    hs      : final Q value
    entropy : final S
    ierr    : 0=OK, 1=not converged, 2=not monotonic, 4=stuck
    """
    nil  = len(d_obs)
    nol  = len(ps)
    spt  = np.sum(d_obs)
    sup  = np.sum(ps)
    eps  = 1e-30

    del_  = 1.0
    nit   = 0
    hs    = 0.0
    hs0   = 0.0
    ierr  = 0

    dps0  = None
    s     = 0.0
    am_hist = []

    while del_ > ac and nit < nim:
        nit += 1

        # Forward projection: image (nol) → data (nil)
        d_pred = P.T.dot(ps)
        d_pred = np.maximum(d_pred, eps)

        # Measure convergence (no-alpha mode)
        if al <= 0.0:
            rm  = np.sqrt(np.mean((d_obs - d_pred)**2))
            del_ = rm / (np.mean(d_obs) + eps)

        # Likelihood
        h, dph = _likelihood(d_obs, d_pred, weights, ih)
        dha    = P.dot(dph)   # back-projection: data (nil) → image (nol)

        # Adjust delta-H
        su  = np.sum(ps * dha)
        sup = np.sum(ps)
        dha = dha + (-su / sup + spt / sup - 1.0)

        # Entropy
        if al > 0.0:
            chi = _smooth_image(ps, indx, indy, nv, ns)
            chi = np.maximum(chi, eps)
            s, dsa = _entropy(ps, chi, indx, indy, nv, ns)
            su2 = np.sum(ps * dsa)
            dsa = dsa - su2 / sup
        else:
            chi = ps.copy()
            dsa = np.zeros(nol)
            s   = 0.0

        # Q value
        hs0 = hs
        hs  = h + al * s
        if nit > 1 and hs < hs0:
            ierr = 2  # not monotonic (warning, continue)

        # Update direction
        dps = ps * (dha + al * dsa)
        if dps0 is None:
            dps0 = dps.copy()

        # Lucy 2-D acceleration
        if al > 0.0:
            dchi1 = _smooth_image(dps0, indx, indy, nv, ns)
            dchi2 = _smooth_image(dps,  indx, indy, nv, ns)
        else:
            dchi1 = dchi2 = np.zeros(nol)

        fla, fmu = _lucy_step(d_obs, d_pred, ps, dps0, dps, P, chi,
                               dchi1, dchi2, weights, al, ih)

        # Normalise direction
        amu = np.sqrt(fla**2 + fmu**2) + eps
        fla /= amu
        fmu /= amu

        # Positivity limit
        combined = fla * dps0 + fmu * dps
        neg = combined < 0
        if neg.any():
            flm_arr = -combined[neg] / (ps[neg] + eps)
            flm = min(flim, 0.9 / flm_arr.max())
        else:
            flm = flim
        amu = min(amu, flm)

        # Stuck check
        am_hist.append(amu)
        if len(am_hist) >= 20 and max(am_hist[-20:]) < 0.01:
            ierr = 4
            return ps, nit, hs, s, ierr

        # Update map
        up  = fmu * dps + fla * dps0
        ps  = ps + amu * up
        ps  = np.maximum(ps, eps)   # positivity
        dps0 = up.copy()

        # Convergence measure
        if al > 0.0:
            dj   = np.abs(dha + al * dsa) / (np.abs(dha) + al * np.abs(dsa) + eps)
            del_ = np.sqrt(np.mean(dj**2))

        # Renormalise
        su = np.sum(ps)
        if su > 0:
            ps = ps * spt / su

        if callback:
            callback(nit, del_, hs, s)

    if nit >= nim:
        ierr = ierr + 1

    return ps, nit, hs, s, ierr


# ---------------------------------------------------------------------------
# 6.  Phase normalisation, blanking, weights, baseline
# ---------------------------------------------------------------------------

def _blanks(pha: np.ndarray, pb0: float, pb1: float) -> np.ndarray:
    """Return boolean mask: True where phase should be blanked."""
    pi2   = 2.0 * np.pi
    nphmi = int(np.floor(pha.min() / pi2))
    nphma = int(np.floor(pha.max() / pi2))
    npb1  = int(np.floor(pb1 / pi2))
    pb0   = pb0 - npb1 * pi2
    pb1   = pb1 - npb1 * pi2
    bl    = np.zeros(len(pha), dtype=bool)
    for ic in range(nphmi, nphma + 2):
        bl |= (pha > pb0 + ic * pi2) & (pha < pb1 + ic * pi2)
    return bl


def _norms(dat: np.ndarray, bl: np.ndarray) -> np.ndarray:
    """Normalise spectra to constant integrated flux per phase."""
    dat = dat.copy()
    nph, nvp = dat.shape
    wi = np.zeros(nph)
    for i in range(nph):
        if not bl[i]:
            wi[i] = dat[i].sum()
    ws = wi[~bl].mean() if (~bl).any() else 1.0
    wi = wi / ws
    wim = wi.min()
    wlm = 0.3
    wl  = max(0.0, wlm - wim)
    for i in range(nph):
        if wi[i] != 0.0:
            dat[i] *= (1.0 + wl) / (wi[i] + wl)
    return dat


def _build_weights(dat: np.ndarray, dat_base: np.ndarray,
                   bl: np.ndarray, nph: int, nvp: int,
                   eb: Optional[np.ndarray], iw: int) -> np.ndarray:
    """Build 1-D weight array (length nil = nph*nvp)."""
    nil = nph * nvp
    w   = np.ones(nil)
    for i in range(nph):
        for j in range(nvp):
            inp = i + j * nph
            if dat_base[i, j] == 0.0:
                w[inp] = 0.0
            if bl[i]:
                w[inp] = 0.0
    if iw == 1 and eb is not None:
        for i in range(nil):
            if eb.flat[i] != 0:
                w[i] *= dat.flat[i] / eb.flat[i]
    total = w.sum()
    if total > 0:
        w = w / total * nph
    return w


def _baseline(P: csr_matrix, indx, indy, nv, nph, nvp, va,
              dat_flat: np.ndarray, bfa: float = 0.1):
    """Compute baseline flux from a flat circular map.

    P has shape (nil, nol) where nil = nph*nvp.
    Forward projection: d = P.T @ vmb  (shape nil)
    """
    dv  = 2.0 * va / (nv - 1)
    fv2 = nv / 2.0

    vx  = dv * (indx - fv2 - 0.5)
    vy  = dv * (indy - fv2 - 0.5)
    vmb = (np.sqrt(vx**2 + vy**2) < va).astype(np.float64)

    # P is (nol, nil); forward projection: d = P.T @ image_vec (nol→nil)
    datb_flat = P.T.dot(vmb)   # shape (nil,)

    # Value near v=0: Fortran datb(nvp/2*nph), 1-based -> 0-based index
    d0 = datb_flat[(nvp // 2) * nph - 1]
    if d0 == 0:
        d0 = datb_flat.mean() or 1.0
    avg       = dat_flat.mean()
    bf        = bfa * avg
    datb_flat = bf * datb_flat / d0
    vmb       = vmb * bf / d0
    return datb_flat, vmb, bf


# ---------------------------------------------------------------------------
# 7.  Alpha sequence outer loop  (alfs)
# ---------------------------------------------------------------------------

def _next_alpha(al0: float, ial: int, alf: float,
                rr: float, clim: float,
                rr1: float, rr2: float) -> float:
    """Spruit's adaptive alpha schedule."""
    if rr > np.sqrt(alf) * clim or ial < 3:
        return al0 / alf**(ial - 1)
    dl2 = np.log(rr1 / (rr2 + 1e-30) + 1e-30)
    dl  = np.log(rr / (clim + 1e-30) + 1e-30)
    # Guard (not in Fortran alfa): if chi did not decrease over the last
    # step (dl2 <= 0, e.g. at the data's noise floor) the log-linear
    # interpolation extrapolates backwards and alpha explodes (seen: 3e59).
    # Fall back to the plain geometric step in that case.
    if dl2 <= 0.0:
        return al0 / alf**(ial - 1)
    x   = dl / dl2
    if x > 1.0 or x <= 0.0:
        return al0 / alf**(ial - 1)
    alfr = np.exp(x * np.log(alf))
    return al0 / alf**(ial - 2) / alfr


# ---------------------------------------------------------------------------
# 8.  Top-level entry point
# ---------------------------------------------------------------------------

def doppler_mem(pha: np.ndarray, dpha: np.ndarray, vp: np.ndarray,
                dat: np.ndarray,
                # Algorithm parameters
                nv: int         = 60,
                ns: int         = 7,
                ac: float       = 8e-4,
                nim: int        = 150,
                al0: float      = 2e-3,
                alf: float      = 1.7,
                nal: int        = 0,
                clim: float     = 1.6,
                # Data options
                ih: int         = 0,
                iw: int         = 0,
                pb0: float      = 0.95,
                pb1: float      = 1.05,
                norm: int       = 1,
                wid: float      = 1e10,
                af: float       = 0.0,
                eb: Optional[np.ndarray] = None,
                # Progress
                verbose: bool   = True,
                progress_fn     = None,
                ) -> MEMResult:
    """
    Run MEM Doppler tomography in pure Python.

    Parameters
    ----------
    pha   : (nph,)      Mid-exposure phases in radians
    dpha  : (nph,)      Phase bin widths in radians
    vp    : (nvp,)      Velocity axis of input spectra in cm/s
    dat   : (nph, nvp)  Normalised input spectra (emission line, bgd-subtracted)
    nv    : int         Velocity map grid size (nv × nv pixels)
    ns    : int         Gaussian smoothing half-width for default map
    ac    : float       Convergence criterion
    nim   : int         Max iterations per alpha
    al0   : float       Starting regularisation alpha
    alf   : float       Alpha reduction factor between steps
    nal   : int         Max number of alphas (0 ⟹ use single alpha; <0 ⟹ abs(nal) alphas regardless)
    clim  : float       Target reduced chi² (or chi equivalent)
    ih    : int         Likelihood: 0=log, 1=chi²
    iw    : int         Use supplied error bars: 0=no, 1=yes
    pb0,pb1 : float     Phase range to blank (e.g. for eclipse)
    norm  : int         1=normalise spectra to constant integrated flux
    eb    : optional    Error bar array (nph, nvp)
    verbose : bool      Print progress
    progress_fn : callable(alpha_idx, nit, del_, rr)  — optional callback

    Returns
    -------
    MEMResult
    """
    pi2  = 2.0 * np.pi
    eps  = 1e-30

    nph, nvp = dat.shape
    nil      = nph * nvp

    # --- convert phases to radians if needed --
    pha  = np.asarray(pha, dtype=np.float64)
    dpha = np.asarray(dpha, dtype=np.float64)
    vp   = np.asarray(vp,  dtype=np.float64)
    dat  = np.asarray(dat, dtype=np.float64).copy()

    # --- default alpha if 0 nal ---
    if nal == 0:
        nal = -1   # use single alpha (Spruit convention)

    # --- blanking array ---
    # convert pb0,pb1 from normalised phase to radians
    bl = _blanks(pha, pb0 * pi2, pb1 * pi2)

    # --- auto alpha ---
    if al0 <= 0.0:
        al0 = 0.003 * (nv / nvp)**2 * min(nph, nvp) / nvp

    # --- build projection matrix ---
    P, indx, indy, va = build_projection(pha, dpha, vp, nv)
    nol = len(indx)

    # --- flatten data in Spruit's order: index = ip + iv*nph ---
    # (Fortran dat(nph,nvp) is column-major, i.e. C-order of dat.T)
    dat_flat = dat.T.ravel().copy()

    # --- baseline (Fortran basel), computed from the raw spectra ---
    datb_flat, vmb, bf = _baseline(P, indx, indy, nv, nph, nvp, va, dat_flat)

    # Add baseline; floor points below bf/3 (Fortran main program)
    dat_work = dat_flat + datb_flat
    low      = dat_work < bf / 3.0
    dat_work = np.where(low, 1e-5 * bf, dat_work)

    # --- phase normalisation, applied after the baseline as in dop.f ---
    if norm == 1:
        dat_work = _norms(dat_work.reshape(nvp, nph).T, bl).T.ravel()

    # --- build weights ---
    dat_base_2d = datb_flat.reshape(nvp, nph).T   # (nph, nvp)
    w = _build_weights(dat_work.reshape(nvp, nph).T, dat_base_2d, bl, nph, nvp,
                       eb.T.ravel() if eb is not None else None, iw)

    # --- eb relative ---
    if eb is not None:
        eb_flat = (eb.T).ravel()
        with np.errstate(divide='ignore', invalid='ignore'):
            eb_flat = np.where(dat_work > 0, eb_flat / dat_work, 0.0)
    else:
        eb_flat = np.ones(nil)

    # --- initialise map (flat) ---
    sdat = dat_work.sum()
    ps   = np.full(nol, sdat / nol)

    # --- convergence metric helpers ---
    def rrms(d1, d0):
        return np.sqrt(np.sum((d1 - d0)**2) / max(len(d1) - 1, 1))

    def rmsn(d1, d0):
        d1r = d1.reshape(nvp, nph)
        d0r = d0.reshape(nvp, nph)
        s = 0.0
        c = 0
        for iv in range(1, nvp - 1):
            for ip in range(1, nph - 1):
                s0 = (d0r[iv-1, ip-1] + d0r[iv-1, ip] + d0r[iv-1, ip+1] +
                      d0r[iv,   ip-1]                  + d0r[iv,   ip+1] +
                      d0r[iv+1, ip-1] + d0r[iv+1, ip] + d0r[iv+1, ip+1])
                s1 = (d1r[iv-1, ip-1] + d1r[iv-1, ip] + d1r[iv-1, ip+1] +
                      d1r[iv,   ip-1]                  + d1r[iv,   ip+1] +
                      d1r[iv+1, ip-1] + d1r[iv+1, ip] + d1r[iv+1, ip+1])
                s += ((s1 - s0) / 8.0 - (d1r[iv, ip] - d0r[iv, ip]))**2
                c += 1
        return np.sqrt(s / max(c, 1))

    # --- outer alpha loop ---
    nal1    = abs(nal)
    rr      = 999.0   # large sentinel; updated after first alpha step
    rr1     = rr
    rr2     = rr
    best_ps   = ps.copy()
    best_dat1 = None
    alx       = al0
    sx        = 0.0
    nco       = 0
    ier       = 0
    converged = False
    last_rr   = rr   # track the actual last chi

    for ial in range(1, nal1 + 1):
        if nal < 0:
            al = al0 * alf**(1 - ial)
        else:
            al = _next_alpha(al0, ial, alf, rr, clim, rr1, rr2)

        def _cb(nit_, del__, hs_, s_):
            if progress_fn:
                progress_fn(ial, nit_, del__, rr)

        ps_new, nit, hs, s, ierr_inner = _rlfldn(
            dat_work, P, ps.copy(), indx, indy, nv,
            al, ac, nim, ns, ih, w,
            callback=_cb if (progress_fn or verbose) else None
        )

        dat1 = P.T.dot(ps_new)   # forward projection → reconstructed data

        # Convergence
        rr0  = rr
        rr_abs = rrms(dat_work, dat1)
        if iw == 0:
            rn = rmsn(dat_work, dat1)
            rr = rr_abs / max(rn, eps)
        else:
            rr = rr_abs

        rr1 = rr2
        rr2 = rr

        nco  += 1
        best_ps   = ps_new.copy()
        best_dat1 = dat1.copy()
        alx  = al
        sx   = s
        last_rr = rr

        if verbose:
            _emit_alpha_line(ial, al, nit, rr, clim, ierr_inner)

        cont = (nco == 0) or (ierr_inner == 0) or (ierr_inner == 2)

        if rr <= clim and nal >= 0:
            converged = True
            break

        if not cont:
            ier = 1
            rr  = rr0
            break

    if nco == 0:
        ier = 2

    # --- build output image, baseline map removed (Fortran: vmo - vmb) ---
    # (kept in data units; Fortran's extra /dv^2 * dnu scaling is omitted)
    image = np.zeros((nv, nv))
    image[indx - 1, indy - 1] = best_ps - vmb

    # Input and reconstruction both with the baseline removed, on the same
    # (normalised) scale -- as written to dop.out by the Fortran.
    if best_dat1 is None:
        best_dat1 = datb_flat.copy()
    dm_reco = (best_dat1 - datb_flat).reshape(nvp, nph).T
    dm_in   = (dat_work  - datb_flat).reshape(nvp, nph).T

    return MEMResult(
        image     = image,
        va        = va / 1e5,          # → km/s
        vp        = vp / 1e5,          # → km/s
        dm        = dm_in,
        dm_reco   = dm_reco,
        phases    = pha / pi2,         # → normalised phase
        alpha     = alx,
        chi2      = last_rr,
        converged = converged or (last_rr <= clim),
        n_iter    = nit,
        entropy   = sx,
    )


# ---------------------------------------------------------------------------
# 9.  Pretty-printing helper (used by doppler_mem and spruit)
# ---------------------------------------------------------------------------

def _emit_alpha_line(ial, al, nit, rr, clim, ierr):
    """Internal: print a single alpha-step summary line."""
    status = "✓" if rr <= clim else ("!" if ierr else "→")
    print(f"  {status}  α={al:.2e}  nit={nit:4d}  χ={rr:.4f}  (target ≤ {clim:.2f})")
