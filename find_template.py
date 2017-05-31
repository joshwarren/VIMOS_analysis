# ==================================================================
#  Find templates to use in fit
# ==================================================================
# warrenj 20150216 Process to analyse the reduced VIMOS data.
# warrenj 20160913 Ported to python

import numpy as np
import glob
from astropy.io import fits
from scipy import ndimage # for gaussian blur
from ppxf import ppxf
import ppxf_util as util
from errors2 import remove_anomalies, use_templates, determine_goodpixels, \
	get_stellar_templates, get_dataCubeDirectory
from checkcomp import checkcomp
cc= checkcomp()

quiet = True
c = 299792.458


def setup(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000], 
	use_all_temp=False):
# ----------===============================================---------
# ----------============= Input parameters  ===============---------
# ----------===============================================---------
	dir = '%s/Data/vimos'  %(cc.base_dir)
	templatesDirectory = '%s/models/miles_library' % (cc.home_dir)

	FWHM_gal = 2.5/(1+z) # VIMOS documentation (and fits header)
	moments = 4 # number of componants to calc with ppxf (see 
				# keyword moments in ppxf.pro for more details)
	degree = 4 # order of addative Legendre polynomial used to 
		   # correct the template continuum shape during the fit 

## ----------=============== Miles library =================---------
	stellar_templates = get_stellar_templates(galaxy, FWHM_gal, use_all_temp=True)
	templates = stellar_templates.templates
	velscale = stellar_templates.velscale
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = get_dataCubeDirectory(galaxy)
	
	f = fits.open(dataCubeDirectory)
	galaxy_data, header = f[0].data, f[0].header
	galaxy_noise = f[1].data
	galaxy_badpix = f[3].data

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
	galaxy_noise = np.delete(galaxy_noise, rows_to_remove, axis=1)
	galaxy_noise = np.delete(galaxy_noise, cols_to_remove, axis=2)
	galaxy_badpix = np.delete(galaxy_badpix, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_badpix, cols_to_remove, axis=2)

	s = galaxy_data.shape

	# Check for nan is data set.
	galaxy_badpix[np.isnan(galaxy_data)] = 1
	galaxy_data[galaxy_badpix==1] = np.nan
	galaxy_noise[galaxy_badpix==1] = np.nan

	# Collapse to single spectrum
	gal_spec = np.nansum(galaxy_data, axis=(1,2))
	gal_noise = np.sqrt(np.nansum(galaxy_noise**2, axis=(1,2)))
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	gal_spec, lam, cut = remove_anomalies(gal_spec, window=201, repeats=3, lam=lam, 
		set_range=set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
	gal_noise = gal_noise[cut]


	## smooth spectrum to fit with templates resolution
	if FWHM_gal < stellar_templates.FWHM_tem:
		sigma = stellar_templates.FWHM_dif/2.355/CDELT_spec # Sigma difference in pixels
		gal_spec = ndimage.gaussian_filter1d(gal_spec, sigma)
		gal_noise = np.sqrt(ndimage.gaussian_filter1d(gal_noise**2, sigma))


	## rebin spectrum logarthmically
	bin_log, logLam_bin, _ = util.log_rebin(lamRange, gal_spec, velscale=velscale)
	gal_noise, logLam_bin, _ = util.log_rebin(lamRange, gal_noise**2, velscale=velscale)
	gal_noise = np.sqrt(gal_noise)

	noise = gal_noise + 0.0000000000001

	dv = (stellar_templates.logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodpixels = determine_goodpixels(logLam_bin,stellar_templates.lamRange_template,vel, z) 
	lambdaq = np.exp(logLam_bin)
	start = [vel, sig] # starting guess

	return templates, bin_log, noise, velscale, start, goodpixels, moments, \
		degree, dv, lambdaq, not quiet, quiet



## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
def find_template(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000]):
	print '     Finding templates to use'

	templates, bin_log, noise, velscale, start, goodpixels,	moments, degree, dv, \
		lambdaq, plot, quiet = setup(galaxy, z=z, vel=vel, sig=sig, discard=discard, 
		set_range=set_range, use_all_temp=True)


	pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodpixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet)

	with open('%s/Data/vimos/analysis/%s/templates.txt' % (cc.base_dir, galaxy), 'w') as f:
		for i in range(templates.shape[1]):
			if pp.weights[i] != 0.0:
				f.write(str(i) + '   ' + str(pp.weights[i]) + '\n')