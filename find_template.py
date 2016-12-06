# ==================================================================
# Analyse reduced VIMOS data using pPFX
# ==================================================================
# warrenj 20150216 Process to analyse the reduced VIMOS data.
# warrenj 20160913 Ported to python

import numpy as np
import glob
from astropy.io import fits
from scipy import ndimage # for gaussian blur

from ppxf import ppxf
import ppxf_util as util
from errors2 import remove_anomalies, determine_goodpixels
from checkcomp import checkcomp
cc= checkcomp()

quiet = True

def find_template (galaxy, z=0.01, discard=2, set_range=[4200,10000]):
	print '     Finding templates to use'
# ----------===============================================---------
# ----------============= Input parameters  ===============---------
# ----------===============================================---------
	dir = '%s/Data/vimos/'  %(cc.base_dir)
	templatesDirectory = '%s/models/miles_library' % (cc.home_dir)

	c = 299792.458
	vel = 0.0 # Initial estimate of the galaxy velocity and
	sig = 200.0 # velocity dispersion in km/s in the rest frame
	FWHM_gal = 2.5/(1+z) # VIMOS documentation (and fits header)
	moments = 4 # number of componants to calc with ppxf (see 
				# keyword moments in ppxf.pro for more details)
	degree = 4 # order of addative Legendre polynomial used to 
		   # correct the template continuum shape during the fit 
	# File for output: an array containing the calculated dynamics of the
	# galaxy. 
	output_temp_weighting = '%s/analysis/%s/templates.txt' % (dir,
		galaxy)

	# Tessellation input
	tessellation_File = '%s/analysis/%s/voronoi_2d_binning_output.txt' % (
		dir, galaxy)


# ----------===============================================---------
# ----------=============== Run analysis  =================---------
# ----------===============================================---------

## ----------=============== Miles library =================---------
	# Finding the template files
	templateFiles = glob.glob(templatesDirectory + '/m0[0-9][0-9][0-9]V') 

	# v1 is wavelength, v2 is spectrum
	v1, v2 = np.loadtxt(templateFiles[0], unpack='True')

	# Using same keywords as fits headers
	CRVAL_temp = v1[0]		# starting wavelength
	NAXIS_temp = np.shape(v2)[0]   # Number of entries
	# wavelength increments (resolution?)
	CDELT_temp = (v1[NAXIS_temp-1]-v1[0])/(NAXIS_temp-1)

	lamRange_template = CRVAL_temp + [0, CDELT_temp*(NAXIS_temp-1)]

	log_temp_template, logLam_template, velscale = \
		util.log_rebin(lamRange_template, v1)

	FWHM_tem = 2.5 # Miles library has FWHM of 2.5A.
	FWHM_dif = np.sqrt(abs(FWHM_gal**2 - FWHM_tem**2))

	nfiles = len(templateFiles)
	templates = np.zeros((len(log_temp_template), nfiles))

	## Reading the contents of the files into the array templates. 
	## Including rebinning them.
	for i in range(nfiles):
		v1, v2 = np.loadtxt(templateFiles[i], unpack='True')
		if FWHM_tem < FWHM_gal:
			sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels
			v2 = ndimage.gaussian_filter1d(v2,sigma)
		## Rebinning templates logarthmically
		log_temp_template, logLam_template, _ = util.log_rebin(lamRange_template, 
			v2, velscale=velscale)
	## ****************************************************************
	## ^^^ this has changed from the last time we called this: we called v1 before...
		templates[:,i] = log_temp_template

	#templates /= np.median(log_temp_template)
## ----------========= Reading Tessellation  ===============---------

	## Reads the txt file containing the output of the binning_spaxels
	## routine. 
	x,y,bin_num = np.loadtxt(tessellation_File, usecols=(0,1,2), \
		unpack=True, skiprows=1)

	n_bins = max(bin_num) + 1
	## Contains the order of the bin numbers in terms of index number.
	order = np.sort(bin_num)
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = glob.glob(dir+"cubes/%s.cube.combined.fits" % (galaxy)) 
		
	galaxy_data, header = fits.getdata(dataCubeDirectory[0], 0, header=True)
	galaxy_badpix = fits.getdata(dataCubeDirectory[0], 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	s = galaxy_data.shape

	rows_to_remove = range(discard)
	rows_to_remove.extend([s[1]-1-i for i in range(discard)])
	cols_to_remove = range(discard)
	cols_to_remove.extend([s[2]-1-i for i in range(discard)])

	galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
	galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
	galaxy_badpix = np.delete(galaxy_badpix, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_badpix, cols_to_remove, axis=2)

	# Check for nan is data set.
	galaxy_badpix[np.isnan(galaxy_data)] = 1
	galaxy_data[galaxy_badpix==1] = 0

	s = galaxy_data.shape

	# Collapse to single spectrum
	gal_spec = np.nansum(np.nansum(galaxy_data, axis=2),axis=1)

	n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	gal_spec, lam = remove_anomalies(gal_spec, window=201, repeats=3, lam=lam, 
		set_range=set_range)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)

	## smooth spectrum to fit with templates resolution
	if FWHM_gal < FWHM_tem:
		sigma = FWHM_dif/2.355/CDELT_spec # Sigma difference in pixels
		gal_spec = ndimage.gaussian_filter1d(gal_spec, sigma)	

	## rebin spectrum logarthmically
	bin_log, logLam_bin, _ = util.log_rebin(lamRange, gal_spec, velscale=velscale)

	noise = bin_log*0+1

	dv = (logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodPixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z) 
	lambdaq = np.exp(logLam_bin)
	start = [vel, sig] # starting guess


## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------

	pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet)


	with open('%s/analysis/%s/templates.txt' % (dir, galaxy), 'w') as f:
		for i in range(nfiles):
			if pp.weights[i] != 0.0:
				f.write(str(i) + '   ' + str(pp.weights[i]) + '\n')
