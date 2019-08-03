# PyDoppler

    This is the repository for a python wrapper for Henk Spruit's doppler tomography software.
    This code can

    The original code and IDL interface can be found at:
    *  https://wwwmpa.mpa-garching.mpg.de/~henk/

    AT the moment, there are many features that haven't been implemented. However, the code will be updated
    continuously. If you have any queries/comments/bug report please send an e-mail to

    * jvhs @ st-andrews.ac.uk

  ## Requirements & Installation

  The doppler tomography code is written in fortran. Please ensure that you have a fortran compiler installed, in particular Gfortran.
  At the moment, only gfortran is supported.


  Python >3.0 version is required (I am using 3.7, No tests have been done for backwards compatability with python >2.0).

  YOu can download and install PyDoppler via pip. In a command line, just type:

  ```
  pip install pydoppler
  ```

  If you need to upgrade the package to the latest version, you can do this with
  ```
  pip install pydoppler -upgrade
  ```


  #  Section 1: Load data

  ##  Section 1.1: Load test data

  You can start to test pydoppler with a test dataset kindly provided by J. Echevarria and published
  in Echevarria et al. 2007, AJ, 134, 262 (https://ui.adsabs.harvard.edu/abs/2007AJ....134..262E/abstract).

  It will create a subdirectory from your working directory (called ugem99) which will contain text files
  for each spectra (txtugem40*). The format of each spectrum file is two columns: Wavelength and flux.
  Wavelength is assumed to be in Angstrom.

  In addition, a phase file (ugem0all.fas) is included which contains the name of the spectrum file and the corresponding 
  orbital phase. This is a two column file:

  txtugem4004 0.7150
  txtugem4005 0.7794
         ...

  ```python
  import pydoppler

  pydoppler.test_data()

  ```
  ##  Section 1.1: Load test data

  #  Section 2:  Doppler tomography


  ```python
  import pydoppler

  # Load base object for tomography
  dop = pydoppler.spruit()

  # Basic data for the tomography to work
  dop.object = 'U Gem'
  dop.base_dir = 'ugem99' # Base directory for input spectra
  dop.list = 'ugem0all.fas'		# Name of the input file
  dop.lam0 = 6562.8 # Wavelength zero in units of the original spectra
  dop.dellx = 220  # Delta in pixels
  dop.delta_phase = 0.003
  dop.delw = 35	# size of Doppler map in wavelenght
  dop.overs = 0.3 # between 0-1
  dop.gama = 36.0  # km /s
  dop.nbins = 28


  ```

  ```python
  # Read in maps
  do.Foldspec()

  # Normalise the spectra
  do.Dopin(continnum_band=[6500,6537,6591,6620],
  		 plot_median=False,poly_degree=2)

  # Perform tomography
  do.Syncdop()

  # Read and plot map
  cb,data = do.Dopmap(limits=[0.05,0.99],colorbar=True,cmaps=cm.gist_stern_r,
  					smooth=False,remove_mean=False)

  # Overplot the donor contours, keplerian and ballistic streams
  qm=0.35
  k1 = 107
  inc=70
  m1=1.2
  porb=0.1769061911
  pydoppler.stream(qm,k1,porb,m1,inc)

  # plot trail spectra
  cb2,cb3,dmr = do.Reco(colorbar=True,limits=[.05,0.95],cmaps=cm.gist_stern_r)
  ```
