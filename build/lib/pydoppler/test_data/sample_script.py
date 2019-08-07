import pydoppler
import matplotlib.pyplot as plt
# Import sample data
# <<< COMMENT OUT IF YOU DONT NEED THE TEST DATASET >>>
pydoppler.test_data()

# Load base object for tomography
dop = pydoppler.spruit()

# Basic data for the tomography to work
dop.object = 'U Gem'
dop.base_dir = 'ugem99' # Base directory for input spectra
dop.list = 'ugem0all.fas'		# Name of the input file
dop.lam0 = 6562.8 # Wavelength zero in units of the original spectra
dop.delta_phase = 0.003
dop.delw = 35	# size of Doppler map in wavelenght
dop.overs = 0.3 # between 0-1, Undersampling of the spectra. 1= Full resolution
dop.gama = 36.0  # km /s
dop.nbins = 28

# Read in the individual spectra and orbital phase information
dop.Foldspec()

# Normalise the spectra
dop.Dopin(continnum_band=[6500,6537,6591,6620],
        plot_median=False,poly_degree=2)

# Perform tomography
dop.Syncdop(ndi=0.7,nri=0.9)

# This routine will display the outcome of the Doppler tomography.
# You can overplot contours and streams.
cb,data = dop.Dopmap(limits=[0.05,0.99],colorbar=False,cmaps=plt.cm.magma_r,
                     smooth=False,remove_mean=False)

# Overplot the donor contours, keplerian and ballistic streams
qm=0.35
k1 = 107
inc=70
m1=1.2
porb=0.1769061911

pydoppler.stream(qm,k1,porb,m1,inc)

# Always check that reconstructed spectra looks like the original one. A good
# rule of thumb "If a feature on the Doppler tomogram isn not in the trail,
# most likely its not real!"
cb2,cb3,dmr,dm = dop.Reco(colorbar=False,limits=[.05,0.95],cmaps=plt.cm.magma_r)
