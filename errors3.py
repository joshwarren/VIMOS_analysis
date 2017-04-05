## ==================================================================
## Propergate uncertainty
## ==================================================================
## warrenj 20150216 Process to progerate the uncertainty using Monty
## Carlo methods to get uncertainty in velocity space.
## warrenj 20160118 Python version of errors.pro
## warrenj 20160422 Will futher separate ionised gases.
## warrenj 20160916 Addapted for stellar population work (No addative 
## polynomials, only multaplicative)

import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
from scipy import ndimage # for gaussian blur
import math
import os
import sys
from checkcomp import checkcomp
cc = checkcomp()
if cc.device == -1:
	cc = checkcomp(override='glamdring')
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
# import of matplotlib.pyplot is within errors routine
from errors2 import set_lines, use_templates, determine_goodpixels, \
	remove_anomalies, get_stellar_templates, get_emission_templates, saveAll
from ppxf import ppxf
import ppxf_util as util


#-----------------------------------------------------------------------------
def set_params():
	quiet = True
	gas = 3 # 0   No gas emission lines
			# 1   Probe ionised gas
			# 2   Seperate gases heated by shocks (OIII and NI) and by SF gas
			#     (Hb and Hd)
			# 3   All gas seperate.
	reps = 0 ## number of monte carlo reps per bin.
	discard = 0
	set_range = np.array([4200,10000])
	FWHM_gal = 2.5 # VIMOS documentation (and fits header)
	stellar_moments = 4 # number of componants to calc with ppxf (see 
						# keyword moments in ppxf.pro for more details)
	gas_moments = 2
	degree = -1  # order of addative Legendre polynomial used to 
				#; correct the template continuum shape during the fit
	mdegree = 10 # order of multaplicative Legendre polynomials used. 
	return quiet, gas, reps, discard, set_range, FWHM_gal, stellar_moments, \
		gas_moments, degree, mdegree
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def errors3(i_gal=None, bin=None):
	if i_gal is None: i_gal=int(sys.argv[1])
	if bin is None: bin=int(sys.argv[2])
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
	quiet, gas, reps, discard, set_range, FWHM_gal, stellar_moments, gas_moments, \
		degree, mdegree = set_params()
	
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[i_gal]

	c = 299792.458

	if cc.device == 'glamdring':
		dir = cc.base_dir
	else:
		dir = '%s/Data/vimos' % (cc.base_dir)


	data_file = dir + "/analysis/galaxies.txt"
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2,3))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]



	tessellation_File = "%s/analysis/%s/" % (dir, galaxy) +\
		"voronoi_2d_binning_output_pop.txt"

	FWHM_gal = FWHM_gal/(1+z) # Adjust resolution in Angstrom

## ----------===============================================---------
## ----------=============== Run analysis  =================---------
## ----------===============================================---------
	stellar_templates = get_stellar_templates(galaxy, FWHM_gal)
	velscale = stellar_templates.velscale
## ----------========= Reading Tessellation  ===============---------

	## Reads the txt file containing the output of the binning_spaxels
	## routine. 
	x,y,bin_num = np.loadtxt(tessellation_File, usecols=(0,1,2), \
		unpack=True, skiprows=1).astype(int)#, dtype='int,int,int')

	n_bins = max(bin_num) + 1
	## Contains the order of the bin numbers in terms of index number.
	order = np.sort(bin_num)
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = "%s/cubes/%s.cube.combined.corr.fits" % (dir,galaxy)
		
	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	# Normalise each spaxel for population pipeline
	norm = np.nanmedian(galaxy_data, axis=0)
	galaxy_data /= norm
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_noise /= norm
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

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
	# galaxy_badpix[np.isnan(galaxy_data)] = 1
	# galaxy_data[galaxy_badpix==1] = 0
	# galaxy_noise[galaxy_badpix==1] = 0.000000001
## ----------============= Spatially Binning ===============---------
	spaxels_in_bin = np.where(bin_num == bin)[0]
	n_spaxels_in_bin = len(spaxels_in_bin)

	bin_lin = np.nansum(galaxy_data[:,x[spaxels_in_bin], y[spaxels_in_bin]], axis=1)
	bin_lin_noise = np.nansum(galaxy_noise[:,x[spaxels_in_bin],
		y[spaxels_in_bin]]**2, axis=1)
	bin_lin_noise = np.sqrt(bin_lin_noise)
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = remove_anomalies(bin_lin, window=201, repeats=3, 
		lam=lam, set_range=set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
	bin_lin_noise = bin_lin_noise[cut]

	## smooth spectrum to fit with templates resolution
	if FWHM_gal < stellar_templates.FWHM_tem:
		sigma = stellar_templates.FWHM_dif/2.355/CDELT_spec # Sigma difference in pixels
		bin_lin = ndimage.gaussian_filter1d(bin_lin, sigma)
		bin_lin_noise = np.sqrt(ndimage.gaussian_filter1d(bin_lin_noise**2, sigma))
	
	## rebin spectrum logarthmically
	bin_log, logLam_bin, _ = util.log_rebin(lamRange, bin_lin, velscale=velscale)
	bin_log_noise, logLam_bin, _ = util.log_rebin(lamRange, bin_lin_noise**2, 
		velscale=velscale)
	bin_log_noise = np.sqrt(bin_log_noise)

	noise = bin_log_noise+0.0000000000001

	dv = (stellar_templates.logLam_template[0]-logLam_bin[0])*c # km/s
	lambdaq = np.exp(logLam_bin)

## ----------===============================================---------
## ----------=============== Emission lines ================---------
## ----------===============================================---------
	e_templates = get_emission_templates(gas, lamRange, 
		stellar_templates.logLam_template, FWHM_gal)

	if gas:
		templates = np.column_stack((stellar_templates.templates, e_templates.templates))
	else:
		templates = stellar_templates.templates
	component = [0]*len(stellar_templates.templatesToUse) + e_templates.component
	templatesToUse = np.append(stellar_templates.templatesToUse, 
		e_templates.templatesToUse)
	element = ['stellar'] + e_templates.element

	start = [[vel, sig]] * (max(component) + 1)
	moments = [stellar_moments] + [gas_moments] * max(component)


	goodPixels = determine_goodpixels(logLam_bin,stellar_templates.lamRange_template,
		vel, z, gas=gas!=0)
## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
	noise = np.abs(noise)
	bin_log_sav = bin_log
	noise_sav = noise
	saveTo="%s/analysis/%s/pop_MC/bestfit/plots/%s.png" % (dir, galaxy, str(bin))

	pp = ppxf(templates, bin_log, noise, velscale, start, 
			  goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, 
			  component=component, lam=lambdaq, plot=not quiet, 
			  quiet=quiet, save=saveTo, mdegree=mdegree)

## ----------===============================================---------
## ----------================= The MC part =================---------
## ----------===============================================---------
	stellar_output = np.zeros((reps, stellar_moments))
	stellar_errors = np.zeros((reps, stellar_moments))
	if gas:
		gas_output = np.zeros((gas, reps, gas_moments))
		gas_errors = np.zeros((gas, reps, gas_moments))

	for rep in range(reps):
		# print rep
		random = np.random.randn(len(noise))
		# gaussian = 1/(np.sqrt(2*math.pi)*noise)*np.exp(-0.5*((random)/noise)**2)
		# add_noise = (random/abs(random))* \
		#     np.sqrt((-2*np.power(noise,2))*np.log(gaussian*noise))
		add_noise = random*np.abs(noise)
		bin_log = pp.bestfit + add_noise
	
		ppMC = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet, bias=0.1, 
			component=component, mdegree=mdegree)

		stellar_output[rep,:] = ppMC.sol[0][0:stellar_moments]
		stellar_errors[rep,:] = ppMC.error[0][0:stellar_moments]
		for g in range(len(element)-1):
			gas_output[g,rep,:] = ppMC.sol[g+1][0:gas_moments]
			gas_errors[g,rep,:] = ppMC.error[g+1][0:gas_moments]

## ----------============ Write ouputs to file =============---------
	saveAll(galaxy, bin, pp, lambdaq, stellar_output, stellar_errors, bin_log_sav, 
		noise_sav, element, templatesToUse, gas_output=gas_output, 
		gas_errors=gas_errors, opt='pop')


##############################################################################


# Use of plot_results.py


if __name__ == '__main__':
	errors3(5,29) if len(sys.argv)<3 else errors3()




