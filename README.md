# PyDoppler

  This is the repository for a python wrapper for Henk Spruit's doppler tomography software.
  This code can will produce a trail spectra of a dataset, and perform
  Doppler tomography maps. It is intended to be a light-weight for single lien datasets.
  The code will be able to:


  The original code and IDL interface can be found at:
   *  https://wwwmpa.mpa-garching.mpg.de/~henk/

  At the moment, there are many features that haven't been implemented. However, the code will be updated
  continuously. If you have any queries/comments/bug report please send an e-mail to:
   * jvhs1 (at) st-andrews.ac.uk

  If you make use of this software, please acknowledge the original code and this repository:
   * Spruit 1998, arXiv, astro-ph/9806141 (https://ui.adsabs.harvard.edu/abs/1998astro.ph..6141S/abstract)

  ## Requirements & Installation

  The doppler tomography code is written in fortran. Please ensure that you have a fortran compiler installed.
  At the moment, only gfortran is supported.


  Python >3.0 version is required (No tests have been done for backwards compatibility with python 2.X).

  You can download and install PyDoppler via pip. In a command line, just type:

  ```
  pip install pydoppler
  ```

  If you need to upgrade the package to the latest version, you can do this with
  ```
  pip install pydoppler -upgrade
  ```


  ##  Section 1: Load data

  ###  Section 1.1: Test case - accreting white dwarf U Gem

  You can start to test PyDoppler with a test dataset kindly provided by J. Echevarria and published
  in Echevarria et al. 2007, AJ, 134, 262 (https://ui.adsabs.harvard.edu/abs/2007AJ....134..262E/abstract).

  It will create a subdirectory from your working directory (called ugem99) which will contain text files
  for each spectra (txtugem40*). The format of each spectrum file is two columns: Wavelength and flux.
  Wavelength is assumed to be in Angstrom.

  In addition, a phase file (ugem0all.fas) is included which contains the name of the spectrum file and the corresponding
  orbital phase. This is a two column file:
```
  txtugem4004 0.7150
  txtugem4005 0.7794
         .
         .
         .
```

  ```python
  import pydoppler

  pydoppler.test_data()

  ```
  ###  Section 1.1: Load test data
  The python wrapper is still in an initial release. I recommend to stick to the file formats
  and directory tree in order for PyDoppler to work properly.

  ```
  wrk_dir
  └── data_dir (ugem99)
  │   ├── individual_spectra (N spectra)
  │   ├── phases_file
  └── fortran_code
  ```
  ##  Section 2:  Doppler tomography
  Before running any routines, verify that you have added all the relevant
  parameters into the PyDoppler object

  ```python
  import pydoppler

  # Load base object for tomography
  dop = pydoppler.spruit()

  # Basic data for the tomography to work
  dop.object = 'U Gem'
  dop.base_dir = 'ugem99' # Base directory for input spectra
  dop.list = 'ugem0all.fas'		# Name of the input file
  dop.lam0 = 6562.8 # Wavelength zero in units of the original spectra
  dop.delta_phase = 0.003  # Exposure time in terms of orbital phase
  dop.delw = 35	# size of Doppler map in wavelength
  dop.overs = 0.3 # between 0-1, Undersampling of the spectra. 1= Full resolution
  dop.gama = 36.0  # km /s
  dop.nbins = 28  # Number of bins. Only supported the number of spectra at the moment
  ```

  ### Section 2.1: Trail spectra of original data
  This routine reads in the raw data and prepares the files for further
  processing.
  ```python
  # Read in the individual spectra and orbital phase information
  dop.Foldspec()
  ```
  ```
  001 txhugem4004  0.715 2048
  002 txhugem4005  0.7794 2048
  003 txhugem4006  0.8348 2048
  004 txhugem4007  0.8942 2048
  005 txhugem4008  0.9518 2048
  006 txhugem4009  0.0072 2048
  007 txhugem4010  0.0632 2048
  008 txhugem4011  0.1186 2048
  009 txhugem4012  0.1745 2048
  010 txhugem4013  0.2344 2048
  011 txhugem4014  0.2904 2048
  012 txhugem4015  0.3724 2048
  013 txhugem4016  0.4283 2048
  014 txhugem4017  0.4866 2048
  015 txhugem4018  0.5425 2048
  016 txhugem4019  0.5979 2048
  017 txhugem4020  0.6544 2048
  018 txhugem4021  0.7098 2048
  019 txhugem4022  0.7652 2048
  020 txhugem4023  0.8195 2048
  021 txhugem4024  0.8772 2048
  022 txhugem4025  0.9269 2048
  023 txhugem4026  0.9614 2048
  024 txhugem4027  0.9959 2048
  025 txhugem4028  0.0304 2048
  026 txhugem4029  0.0648 2048
  027 txhugem4030  0.1027 2048
  028 txhugem4031  0.1372 2048
```

  ### Section 2.2: Normalise the data and set doppler files
  You will need to define a continnum band - one at each side of the emission line -
  to fit and later subtract the continuum. This normalised spectra will be put in
  in a dopin file to be read by the fortran code.
  ```python  
  # Normalise the spectra
      dop.Dopin(continnum_band=[6500,6537,6591,6620],
      		 plot_median=False,poly_degree=2)
  ```


 <div class="row">
   <div class="column">
     <img src="pydoppler/test_data/output_images/Average_Spec.png" width="333" height="450" />
   </div>
   <div class="column">
     <img src="pydoppler/test_data/output_images/Trail.png" width="333" height="450" />
   </div>
 </div>
  ### Section 2.3: Run the fortran code
  Now, let's run the tomography software!
  ```python
  # Perform tomography
  dop.Syncdop()
  ```
  ### Section 2.4: Plot the tomography map
  This routine will display the outcome of the Doppler tomography. You can overplot
  contours and streams.
  ```python
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
  ```
  ### Section 2.5: Spectra reconstrunction
  Always check that reconstructed spectra looks like the original one. A good
  rule of thumb "If a feature on the Doppler tomogram isn not in the trail, most likely
  its not real!"

  ```python
  # plot trail spectra
  cb2,cb3,dmr = do.Reco(colorbar=True,limits=[.05,0.95],cmaps=cm.gist_stern_r)
  ```
  ## Section 3: Troubleshoot
  This is an early version of the wrapper. Things will go wrong. If you find a
  bug or need a feature, I will try my best to work it out. If you think you can
  add to this wrapper, I encourage you push new changes and contribute to it.
