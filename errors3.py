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
# import of matplotlib.pyplot is within errors routine
from errors2 import set_lines, use_templates, determine_goodpixels, \
	remove_anomalies, get_stellar_templates, get_emission_templates


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
		templatesDirectory = '%s/ppxf/MILES_library/' % (cc.base_dir)	
	else:
		dir = '%s/Data/vimos' % (cc.base_dir)
		templatesDirectory = '%s/models/miles_library/' % (cc.home_dir)
	
	if cc.remote:
		import matplotlib # 20160202 JP to stop lack-of X-windows error
		matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
		import matplotlib.pyplot as plt # used for plotting
	else:
		import matplotlib.pyplot as plt # used for plotting
	from ppxf import ppxf
	import ppxf_util as util



	data_file = dir + "/analysis/galaxies.txt"
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]



	tessellation_File = "%s/analysis/%s/" % (dir, galaxy) +\
		"voronoi_2d_binning_output_pop.txt"
	tessellation_File2 = "%s/analysis/%s/" % (dir, galaxy) +\
		"voronoi_2d_binning_output2_pop.txt"


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

	dataCubeDirectory = glob.glob(dir+"/cubes/%s.cube.combined.corr.fits" % (galaxy)) 
		
	galaxy_data, header = fits.getdata(dataCubeDirectory[0], 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory[0], 1)
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

	bin_lin = np.nansum(galaxy_data[:,x[spaxels_in_bin],
		y[spaxels_in_bin]], axis=1)
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
	templatesToUse = np.append(stellar_templates.templatesToUse, e_templates.templatesToUse)
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

		stellar_output[rep,:] = ppMC.sol[0:stellar_moments][0]
		stellar_errors[rep,:] = ppMC.error[0:stellar_moments][0]
		for g in range(gas):
			gas_output[g,rep,:] = ppMC.sol[0:gas_moments][g]
			gas_errors[g,rep,:] = ppMC.error[0:gas_moments][g]

## ----------============ Write ouputs to file =============---------
	# stellar MC results
	if not os.path.exists("%s/analysis/%s/pop_MC/stellar/errors" % (dir,galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/stellar/errors" % (dir, galaxy))
	bin_file = "%s/analysis/%s/pop_MC/stellar/%s.dat" % (dir, galaxy, str(bin))
	errors_file = "%s/analysis/%s/pop_MC/stellar/errors/%s.dat" % (dir, galaxy, 
		str(bin))
	f = open(bin_file, 'w')
	e = open(errors_file, 'w')
	for i in range(reps):
		f.write(str(stellar_output[i,0]) + "   " + \
			str(stellar_output[i,1]) + "   " + str(stellar_output[i,2]) + \
			"   " + str(stellar_output[i,3]) + '\n')
		e.write(str(stellar_errors[i,0]) + "   " + str(stellar_errors[i,1]) + \
			"   " + str(stellar_errors[i,2]) + "   " + \
			str(stellar_errors[i,3]) + '\n')

	# gas MC results
	gas_dir=[]
	if gas != 0: gas_dir = e_templates.element
	for d in range(len(gas_dir)):
		if not os.path.exists("%s/analysis/%s/pop_MC/gas/%s/errors" % (dir,
			galaxy, gas_dir[d])):
			os.makedirs("%s/analysis/%s/pop_MC/gas/%s/errors" % (dir,
				galaxy, gas_dir[d]))

		gas_file = "%s/analysis/%s/pop_MC/gas/%s/%s.dat" % (dir, galaxy,
			gas_dir[d], str(bin))
		gas_errors_file = "%s/analysis/%s/pop_MC/gas/%s/errors/%s.dat" % (
			dir, galaxy, gas_dir[d], str(bin))

		g = open(gas_file, 'w')
		ger = open(gas_errors_file, 'w')
		
		for i in range(reps):

			g.write(str(gas_output[d,i,0]) + "   " + str(gas_output[d,i,1]) + \
				"   " + str(gas_output[d,i,2]) + "   " + \
				str(gas_output[d,i,3]) + '\n')
			ger.write(str(gas_errors[d,i,0]) + "   " + \
				str(gas_errors[d,i,1]) + "   " + str(gas_errors[d,i,2]) + \
				"   " + str(gas_errors[d,i,3]) + '\n')

			
	## save bestfit spectrum
	if not os.path.exists("%s/analysis/%s/pop_MC/bestfit" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/bestfit" % (dir, galaxy)) 
	bestfit_file = "%s/analysis/%s/pop_MC/bestfit/%s.dat" % (dir, galaxy, 
		str(bin))
   
	s = open(bestfit_file, 'w')
	for i in range(len(pp.bestfit)):
		s.write(str(pp.bestfit[i]) + '\n')

	## save input
	if not os.path.exists("%s/analysis/%s/pop_MC/input" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/input" % (dir, galaxy)) 
	input_file = "%s/analysis/%s/pop_MC/input/%s.dat" % (dir, galaxy, str(bin))
   
	inp = open(input_file, 'w')
	for i in range(len(bin_log_sav)):
		inp.write(str(bin_log_sav[i]) + '\n')
	## save input noise
	if not os.path.exists("%s/analysis/%s/pop_MC/noise_input" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/noise_input" % (dir, galaxy)) 
	input_file = "%s/analysis/%s/pop_MC/noise_input/%s.dat" % (dir, galaxy, str(bin))
   
	inp = open(input_file, 'w')
	for i in range(len(noise_sav)):
		inp.write(str(noise_sav[i]) + '\n')

	## save bestfit LOSVD output
	bestfit_file = "%s/analysis/%s/pop_MC/%s.dat" % (dir, galaxy, str(bin))
   
	b = open(bestfit_file, 'w')
	b.write(element[0] + "   " + str(pp.sol[0][0]) + "   " + str(pp.sol[0][1]) + 
		"   " + str(pp.sol[0][2]) + "   " +  str(pp.sol[0][3]) + '\n')
	if gas:
		for i in range(1, np.shape(pp.sol)[0]):
			b.write(element[i] + "   " + str(pp.sol[i][0]) + "   " + 
				str(pp.sol[i][1]) + '   0.0   0.0 \n')
	else: b.write("stellar " + str(pp.sol[0]) + "   " + str(pp.sol[1]) + 
		"   " + str(pp.sol[2]) + "   " + str(pp.sol[3]) + '\n')

	## save chi2
	if not os.path.exists("%s/analysis/%s/pop_MC/chi2" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/chi2" % (dir, galaxy)) 
	chi2_file = "%s/analysis/%s/pop_MC/chi2/%s.dat" % (dir, galaxy, str(bin))
   
	c2 = open(chi2_file, 'w')
	c2.write(str(pp.chi2) + '\n')

	## save weights
	if not os.path.exists("%s/analysis/%s/pop_MC/temp_weights" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/temp_weights" % (dir, galaxy)) 
	weights_file = "%s/analysis/%s/pop_MC/temp_weights/%s.dat" % (dir, galaxy,
		str(bin))

	w = open(weights_file, 'w')
	for i in range(len(pp.weights)):
		w.write(str(templatesToUse[i]) + "   " + str(pp.weights[i]) + '\n') 

	## save indervidual template bestfits
	if not os.path.exists("%s/analysis/%s/pop_MC/bestfit/matrix" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/bestfit/matrix" % (dir, galaxy)) 
	matrix_file = "%s/analysis/%s/pop_MC/bestfit/matrix/%s.dat" % (dir, galaxy,
		str(bin))
	## save addative polyweights
	if hasattr(pp, 'polyweights'):
		if not os.path.exists("%s/analysis/%s/pop_MC/apweights" % (dir, galaxy)):
			os.makedirs("%s/analysis/%s/pop_MC/apweights" % (dir, galaxy)) 
		polyweights_file = "%s/analysis/%s/pop_MC/apweights/%s.dat" % (dir, galaxy,
			str(bin))

		apw = open(polyweights_file, 'w')
		for i in range(len(pp.polyweights)):
			apw.write(str(pp.polyweights[i]) + '\n')

		l = open(matrix_file, 'w')
		for i in range(len(pp.polyweights),pp.matrix.shape[1]):
			l.write(str(templatesToUse[i]) + "   ")
			for j in range(pp.matrix.shape[0]):
				l.write(str(pp.matrix[j,i]) + "   ")
			l.write('\n')
	else:
		l = open(matrix_file, 'w')
		for i in range(pp.matrix.shape[1]):
			l.write(str(templatesToUse[i]) + "   ")
			for j in range(pp.matrix.shape[0]):
				l.write(str(pp.matrix[j,i]) + "   ")
			l.write('\n')

	## save multiplicative polyweights
	if hasattr(pp, 'mpolyweights'):
		if not os.path.exists("%s/analysis/%s/pop_MC/mpweights" % (dir, galaxy)):
			os.makedirs("%s/analysis/%s/pop_MC/mpweights" % (dir, galaxy)) 
		mpolyweights_file = "%s/analysis/%s/pop_MC/mpweights/%s.dat" % (dir, galaxy,
			str(bin))

		mpw = open(mpolyweights_file, 'w')
		for i in range(len(pp.mpolyweights)):
			mpw.write(str(pp.mpolyweights[i]) + '\n')

	## save lambda input
	if not os.path.exists("%s/analysis/%s/pop_MC/lambda" % (dir, galaxy)):
		os.makedirs("%s/analysis/%s/pop_MC/lambda" % (dir, galaxy)) 
	lambda_file = "%s/analysis/%s/pop_MC/lambda/%s.dat" % (dir, galaxy,
		str(bin))

	l = open(lambda_file, 'w')
	for i in range(len(lambdaq)):
		l.write(str(lambdaq[i]) + '\n')


##############################################################################


# Use of plot_results.py


if __name__ == '__main__':
	errors3(5,29) if len(sys.argv)<3 else errors3()




