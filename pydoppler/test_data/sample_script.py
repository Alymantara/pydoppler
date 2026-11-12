"""
PyDoppler native — U Gem sample script.
Run from any directory (pydoppler.test_data() copies the data there):
    python sample_script.py
"""

import pydoppler
import matplotlib.pyplot as plt

pydoppler.test_data()

dop = pydoppler.spruit()
dop.object      = 'U Gem'
dop.base_dir    = 'ugem99'
dop.list        = 'ugem0all.fas'
dop.lam0        = 6562.8
dop.delta_phase = 0.003
dop.delw        = 35
dop.overs       = 0.3
dop.gama        = 36.0
dop.nbins       = 28

# MEM parameters — increase nal/nim for better convergence
dop.nim  = 150
dop.nal  = 10
dop.clim = 1.6
dop.al0  = 2e-3
dop.ns   = 7

dop.Foldspec()
dop.Dopin(continnum_band=[6500, 6537, 6591, 6620],
          plot_median=False, poly_degree=2)
dop.Syncdop()   # nv auto-sized; pass nv=100 for higher resolution

fig, data = dop.Dopmap(colorbar=True)   # defaults: linear stretch, magma_r

pydoppler.stream(0.35, 107, 0.1769061911, 1.2, 70,
                 tidal_radius=True,
                 resonance_32=True,
                 resonance_21=True)
plt.savefig('Doppler_Map.pdf', bbox_inches='tight')

fig2, dm, dmr = dop.Reco(colorbar=True, limits=[0.05, 0.95])
plt.show()
