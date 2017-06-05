## ==================================================================
## Propergate uncertainty
## ==================================================================
## warrenj 20150216 Process to progerate the uncertainty using Monty
## Carlo methods to get uncertainty in velocity space.
## warrenj 20160118 Python version of errors.pro
## warrenj 20160422 Will futher separate ionised gases.

import numpy as np # for array handling
from astropy.io import fits # reads fits files (is from astropy)
from scipy import ndimage # for gaussian blur
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
import ppxf_util as util
from rolling_stats import *

#-----------------------------------------------------------------------------
class set_params(object):
	def __init__(opt='kin'):
		self.quiet = True
		self.gas = 1 # 0   No gas emission lines
				# 1   Probe ionised gas
				# 2   Seperate gases heated by shocks (OIII and NI) and by SF gas
				#     (Hb and Hd)
				# 3   All gas seperate.
		self.reps = 0 ## number of monte carlo reps per bin.
		self.discard = 0
		self.set_range = np.array([4200,10000])
		self.FWHM_gal = 2.5 # VIMOS documentation (and fits header)
		self.stellar_moments = 4 # number of componants to calc with ppxf (see 
							# keyword moments in ppxf.pro for more details)
		self.gas_moments = 2
		if 'kin' in opt:
			self.degree = 4  # order of addative Legendre polynomial used to 
							# correct the template continuum shape during the fit
			self.mdegree = 0
		elif 'pop' in opt:
			self.degree = -1
			self.mdegree = 10
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def set_lines (lines, logLam_temp, FWHM_gal):
	# In this routine all lines are free to have independent intensities.
	# One can fix the intensity ratio of different lines (e.g. the [OIII] doublet)
	# by placing them in the same emission template
	lam = np.exp(logLam_temp)
	#lines = lines[where((lines gt min(lam)) and (lines lt max(lam)))]
	sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom
	emission_lines = np.zeros((len(logLam_temp),len(lines)))
	for j in range(len(lines)):
		emission_lines[:,j] = np.exp(-0.5*np.power((lam - lines[j])/sigma,2))
	return emission_lines
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def use_templates(galaxy, glamdring=False):
	if glamdring:
		template_weighting = '%s/analysis/%s/templates.txt' % (cc.base_dir, 
			galaxy)
	else:
		template_weighting = '%s/Data/vimos/analysis/%s/templates.txt' % (
			cc.base_dir, galaxy)

	templatesToUse = np.loadtxt(template_weighting, usecols=(0,), dtype='i')
	return templatesToUse
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
class get_stellar_templates(object):

	def __init__(self, galaxy, FWHM_gal, use_all_temp=False):
		import glob # for searching for files
		# Finding the template files
		# There is some issue with the memory structure of the university macs (HFS+),
		# meaning these templates can only be loaded once if located on the home directory,
		# but more if on the Data partition...
		if cc.device != 'uni':
			templateFiles = glob.glob('%s/models/miles_library/m0[0-9][0-9][0-9]V' % (
				cc.home_dir))
		else:
			templateFiles = glob.glob('%s/Data/idl_libraries/ppxf/' % (cc.base_dir) +
				'MILES_library/m0[0-9][0-9][0-9]V')

		# self.wav is wavelength, v2 is spectrum
		self.wav, v2 = np.loadtxt(templateFiles[0], unpack='True')

		# Using same keywords as fits headers
		CRVAL_temp = self.wav[0]		# starting wavelength
		NAXIS_temp = np.shape(v2)[0]   # Number of entries
		# wavelength increments (resolution?)
		CDELT_temp = (self.wav[NAXIS_temp-1]-self.wav[0])/(NAXIS_temp-1)

		self.lamRange_template = CRVAL_temp + [0, CDELT_temp*(NAXIS_temp-1)]

		log_temp_template, self.logLam_template, self.velscale = \
			util.log_rebin(self.lamRange_template, self.wav)

		self.FWHM_tem = 2.5 # Miles library has FWHM of 2.5A.
		self.FWHM_dif = np.sqrt(abs(FWHM_gal**2 - self.FWHM_tem**2))
		
		if use_all_temp:
			self.ntemp = len(templateFiles)
			self.templatesToUse = np.arange(1, self.ntemp + 1)
		else:
			self.templatesToUse = use_templates(galaxy, cc.device=='glamdring')
			self.ntemp = len(self.templatesToUse)
		self.templates = np.zeros((len(log_temp_template), self.ntemp))
		self.lin_templates = np.zeros((len(log_temp_template), self.ntemp))

		## Reading the contents of the files into the array templates. 
		## Including rebinning them.
		for i in range(self.ntemp):
			if use_all_temp:
				_, self.lin_templates[:,i] = np.loadtxt(templateFiles[i], unpack='True')
			else:
				_, self.lin_templates[:,i] = np.loadtxt(
					templateFiles[self.templatesToUse[i]], unpack='True')
			if self.FWHM_tem < FWHM_gal:
				sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels
				conv_temp = ndimage.gaussian_filter1d(self.lin_templates[:,i],sigma)
			else:
				conv_temp = self.lin_templates[:,i]
			## Rebinning templates logarthmically
			log_temp_template, self.logLam_template, _ = util.log_rebin(
				self.lamRange_template, conv_temp, velscale=self.velscale)
			self.templates[:,i] = log_temp_template
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
class get_emission_templates(object):
	def __init__(self, gas, lamRange, logLam_template, FWHM_gal, quiet=True):
		self.element = []
		self.component = []
		self.templatesToUse = []
		self.templates = []

		## ----------============ All lines together ===============---------
		if gas == 1:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			for i in range(len(line_name)):
				self.templatesToUse = np.append(self.templatesToUse, line_name[i])
				self.templates.append(emission_lines[:,i])
			self.component = self.component + [1]*len(line_name)
			self.element.append('gas')
		## ----------=============== SF and shocks lines ==============---------
		if gas == 2:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			for i in range(len(line_name)):
				

				if 'H' in line_name[i]:
					self.templatesToUse = np.append(self.templatesToUse, line_name[i])
					self.templates.append(emission_lines[:,i])
					self.component = self.component + [1]
					self.element.append['SF']
				else:
					self.templatesToUse = np.append(self.templatesToUse, line_name[i])
					self.templates.append(emission_lines[:,i])
					self.component = self.component + [2] 
					self.element.append['Shocks']
		## ----------=========== All lines inderpendantly ==============---------
		if gas == 3:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			aph_lin = np.sort(line_name)

			for i in range(len(line_name)):
				self.templatesToUse = np.append(self.templatesToUse, line_name[i])

				# line listed alphabetically
				self.component = self.component + \
					[np.where(line_name[i] == aph_lin)[0][0]+1]
				self.templates.append(emission_lines[:,i])
				self.element.append(aph_lin[i])

		self.templates = np.array(self.templates).T
		if gas:	self.ntemp = self.templates.shape[1]
		else: self.ntemp = 0
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def determine_goodpixels(logLam, lamRangeTemp, vel, z, gas=False):
	"""
	warrenj 20150905 Copied from ppxf_determine_goodPixels.pro
	
	PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of
		goodPixels to be used as input keyword for the routine
		PPXF. This is useful to mask gas emission lines or atmospheric
		absorptions. It can be trivially adapted to mask different
		lines. 
	
	INPUT PARAMETERS:
	- LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
		of each pixel of the log rebinned *galaxy* spectrum.
	- LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the
		minimum and maximum wavelength in Angstrom in the stellar
		*template* used in PPXF. 
	- VEL: Estimate of the galaxy velocity in km/s.
	
	V1.0: Michele Cappellari, Leiden, 9 September 2005
	V1.01: Made a separate routine and included additional common
	  emission lines. MC, Oxford 12 January 2012
	V1.02: Included more lines. MC, Oxford, 7 Januray 2014
	20150617 warrenj Added Telluric lines (tell) at 5199 (is a blended sky
	line)
	20160816 warrenj Added NI doublet 5199.36 and 5201.86
	"""

	c = 299792.458 # speed of light in km/s
 
	# dv = lines*0+800d # width/2 of masked gas emission region in km/s
	dv = 800 # width/2 of masked gas emission region in km/s
	# flag = bytearray([0]*len(logLam))
	flag = logLam < 0

	# Marks telluric line
	tell = 5199
	flag |= (logLam > np.log(tell) - z - dv/c) \
		& (logLam < np.log(tell) - z + dv/c) 

	if not gas:
		#                 -----[OII]-----   Hdelta    Hgamma   Hbeta;
		lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, \
		#    -----[OIII]----- -----[NI]----- ----??----   [OI]   
			4958.92, 5006.84, 5199.36, 5201.86, 5528, 5535, 6300.30, \
		#    -----[NII]------  Halpha   -----[SII]-----   
			6548.03, 6583.41,  6562.80, 6716.47, 6730.85])

		for j in range(len(lines)):
			flag |= (logLam > np.log(lines[j]) + (vel- dv)/c) \
				& (logLam < np.log(lines[j]) + (vel+ dv)/c)


	flag |= logLam < np.log(lamRangeTemp[0]) + (vel + 900)/c # Mask edges of
	flag |= logLam > np.log(lamRangeTemp[1]) + (vel - 900)/c # stellar library


	flag[0:3] = 1 # Mask edge of data
	flag[-4:]= 1 # to remove edge effects
	return np.where(flag == 0)[0]
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def saveAll(galaxy, bin, pp, lambdaq, stellar_output, stellar_errors, bin_log_sav, 
	noise_sav, element, templatesToUse, gas_output=None, gas_errors=None, opt='kin'):
	# stellar MC results
	if cc.device == 'glamdring':
		dir = '%s/analysis/%s/%s/MC' % (cc.base_dir, galaxy, opt)
	else:
		dir = '%s/Data/vimos/analysis/%s/%s/MC' % (cc.base_dir, galaxy, opt)

	reps = stellar_output.shape[0]
	gas = len(element) > 1

	def check_directory(d):
		if not os.path.exists(d):
			try:
				os.makedirs(d)
			except OSError: 
				# Another node has created directory in the time that it took for if 
				# statement above
				pass

	check_directory("%s/stellar/errors" % (dir))
	bin_file = "%s/stellar/%s.dat" % (dir, str(bin))
	errors_file = "%s/stellar/errors/%s.dat" % (dir, str(bin))
	with open(bin_file, 'w') as f,  open(errors_file, 'w') as e:
		for i in range(reps):
			f.write(str(stellar_output[i,0]) + "   " + \
				str(stellar_output[i,1]) + "   " + str(stellar_output[i,2]) + \
				"   " + str(stellar_output[i,3]) + '\n')
			e.write(str(stellar_errors[i,0]) + "   " + str(stellar_errors[i,1]) + \
				"   " + str(stellar_errors[i,2]) + "   " + \
				str(stellar_errors[i,3]) + '\n')
	
	# gas MC results
	if gas: gas_dir = [e for e in element if e != 'stellar']
	else: gas_dir=[] 
	for d in range(len(gas_dir)):
		check_directory("%s/gas/%s/errors" % (dir, gas_dir[d]))
		gas_file = "%s/gas/%s/%s.dat" % (dir, gas_dir[d], str(bin))
		gas_errors_file = "%s/gas/%s/errors/%s.dat" % (dir, gas_dir[d], str(bin))

		with open(gas_file, 'w') as g, open(gas_errors_file, 'w') as ger:
			for i in range(reps):
				if gas_output is not None:
					g.write(str(gas_output[d,i,0]) + "   " + str(gas_output[d,i,1]) + \
						"   " + str(gas_output[d,i,2]) + "   " + \
						str(gas_output[d,i,3]) + '\n')
				if gas_errors is not None:
					ger.write(str(gas_errors[d,i,0]) + "   " + \
						str(gas_errors[d,i,1]) + "   " + str(gas_errors[d,i,2]) + \
						"   " + str(gas_errors[d,i,3]) + '\n')

	## save bestfit spectrum
	check_directory("%s/bestfit" % (dir))
	bestfit_file = "%s/bestfit/%s.dat" % (dir, str(bin))
   
	with open(bestfit_file, 'w') as s:
		for i in range(len(pp.bestfit)):
			s.write(str(pp.bestfit[i]) + '\n')

	## save input
	check_directory("%s/input" % (dir))
	input_file = "%s/input/%s.dat" % (dir, str(bin))
   
	with open(input_file, 'w') as inp:
		for i in range(len(bin_log_sav)):
			inp.write(str(bin_log_sav[i]) + '\n')

	## save input noise
	check_directory("%s/noise_input" % (dir))
	input_file = "%s/noise_input/%s.dat" % (dir, str(bin))
   
	with open(input_file, 'w') as inp:
		for i in range(len(noise_sav)):
			inp.write(str(noise_sav[i]) + '\n')

	## save bestfit LOSVD output
	bestfit_file = "%s/%s.dat" % (dir, str(bin))
	with open(bestfit_file, 'w') as b:
		if gas:
			for i in range(np.shape(pp.sol)[0]):
				b.write(element[i])
				for j in pp.sol[i]:
					b.write("   " + str(j)) 
				b.write('\n')
		else: 
			b.write("stellar")
			for j in range(stellar_moments):
				b.write("   " + str(j))
			b.write('\n')

	## save chi2
	check_directory("%s/chi2" % (dir))
	chi2_file = "%s/chi2/%s.dat" % (dir, str(bin))
   
	with open(chi2_file, 'w') as c2:
		c2.write(str(pp.chi2) + '\n')

	## save weights
	check_directory("%s/temp_weights" % (dir))
	weights_file = "%s/temp_weights/%s.dat" % (dir, str(bin))

	with open(weights_file, 'w') as w:
		for i in range(len(pp.weights)):
			w.write(str(templatesToUse[i]) + "   " + str(pp.weights[i]) + '\n') 

	## save indervidual template bestfits
	check_directory("%s/bestfit/matrix" % (dir))
	matrix_file = "%s/bestfit/matrix/%s.dat" % (dir, str(bin))

	## save addative polyweights
	if hasattr(pp, 'polyweights'):
		check_directory("%s/apweights" % (dir))
		polyweights_file = "%s/apweights/%s.dat" % (dir, str(bin))

		with open(polyweights_file, 'w') as apw:
			for i in range(len(pp.polyweights)):
				apw.write(str(pp.polyweights[i]) + '\n')

		with open(matrix_file, 'w') as l:
			for i, j in enumerate(range(len(pp.polyweights),pp.matrix.shape[1])):
				l.write(str(templatesToUse[i]) + "   ")
				for k in range(pp.matrix.shape[0]):
					l.write(str(pp.matrix[k,j]) + "   ")
				l.write('\n')
	else:
		with open(matrix_file, 'w') as l:
			for i in range(pp.matrix.shape[1]):
				l.write(str(templatesToUse[i]) + "   ")
				for j in range(pp.matrix.shape[0]):
					l.write(str(pp.matrix[j,i]) + "   ")
				l.write('\n')

	## save multiplicative polyweights
	if hasattr(pp, 'mpolyweights'):
		check_directory("%s/mpweights" % (dir))
		mpolyweights_file = "%s/mpweights/%s.dat" % (dir, str(bin))

		with open(mpolyweights_file, 'w') as mpw:
			for i in range(len(pp.mpolyweights)):
				mpw.write(str(pp.mpolyweights[i]) + '\n')

	## save lambda input
	check_directory("%s/lambda" % (dir))
	lambda_file = "%s/lambda/%s.dat" % (dir, str(bin))

	with open(lambda_file, 'w') as l:
		for i in range(len(lambdaq)):
			l.write(str(lambdaq[i]) + '\n')
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def remove_anomalies(spec, window=201, repeats=3, lam=None, set_range=None, 
	return_cuts=False, n_sigma=3):
	if set_range is not None and lam is None:
		raise ValueError('lam keyword must be supplied if set_range keyword'+\
			' is supplied')
	elif set_range is not None and lam is not None:
		lam = np.array(lam)
		r = (lam > set_range[0]) * (lam < set_range[1])
	else:
		r = np.ones(len(spec)).astype(bool)
	spec = np.array(spec[r])
	lam = np.array(lam[r])


	x=np.arange(len(spec))
	mask = np.zeros(len(spec)).astype(bool)
	for rep in range(repeats):
		med = rollmed(spec,window)
		std = rollstd(spec,window)
		for i in range(len(mask)):
			mask[i] += spec[i] > med[i]+n_sigma*std[i]
			mask[i] += spec[i] < med[i]-n_sigma*std[i]

		spec[mask] = np.interp(x[mask],x[~mask],spec[~mask])

	if lam is not None:
		if return_cuts:
			return spec, lam, r
		else:
			return spec, lam
	else:
		if return_cuts:
			return spec, r
		else:
			return spec
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def get_dataCubeDirectory(galaxy):
	class mystring2(str):
		def __init__(self, s):
			str.__init__(s)
			self.RAoffset = 0 # offset in arcsec
			self.decoffset = 0


	class mystring(str):
		def __init__(self, s):
			str.__init__(s)
			self.radio = mystring2('')
			self.CO = mystring2('')
			self.xray = mystring2('')


	if cc.device == 'uni':
		dir = '%s/Data/vimos' % (cc.base_dir)
	elif cc.device == 'glamdring':
		dir = '%s/cubes' % (cc.base_dir)
	elif 'home' in cc.device:
		dir = '%s/Data/vimos' % (cc.base_dir)

	dataCubeDirectory = mystring('%s/cubes/%s.cube.combined.corr.fits' %  (dir, galaxy))
	dataCubeDirectory.CO = mystring2("%s/Data/alma/%s-mom0.fits" % (cc.base_dir, galaxy))

	if galaxy == 'eso443-g024':
		dataCubeDirectory.xray = '%s/Data/Chandra/ESO443_full.fits' % (cc.base_dir)
	elif galaxy == 'ic1459':
		dataCubeDirectory.xray = '%s/Data/Chandra/IC1459_full.fits' % (cc.base_dir)
	elif galaxy == 'ic1531':
		dataCubeDirectory.xray = '%s/Data/Chandra/IC1531_full.fits' % (cc.base_dir)
	elif galaxy == 'ic4296':
		dataCubeDirectory.xray = '%s/Data/Chandra/IC4296_full.fits' % (cc.base_dir)
	elif galaxy == 'ngc0612':
		dataCubeDirectory.xray = '%s/Data/Chandra/N612_full.fits' % (cc.base_dir)
	elif galaxy == 'ngc1399':
		dataCubeDirectory.xray = '%s/Data/Chandra/N1399_full.fits' % (cc.base_dir)
	elif galaxy == 'ngc3100':
		dataCubeDirectory.radio = mystring2('%s/Data/VLA/%s/AD270.fits' % (cc.base_dir, 
			galaxy))
		dataCubeDirectory.radio.RAoffset = 2.6
		dataCubeDirectory.radio.decoffset = -1.2
		dataCubeDirectory.CO.RAoffset = 2.6
		dataCubeDirectory.CO.decoffset = -1.2
	elif galaxy == 'ngc3557':
		pass
	elif galaxy == 'ngc7075':
		pass
	elif galaxy == 'pks0718-34':
		pass
	
	

	return dataCubeDirectory
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def errors2(i_gal=None, opt=None, bin=None):
	if i_gal is None: i_gal=int(sys.argv[1])
	if opt is None: opt=sys.argv[2]
	if bin is None: bin=int(sys.argv[3])
	from ppxf import ppxf
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
	params = set_params()
	
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[i_gal]

	c = 299792.458

	if cc.device == 'glamdring':
		dir = cc.base_dir
	else:
		dir = '%s/Data/vimos' % (cc.base_dir)


	data_file = "%s/analysis/galaxies.txt" % (dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2,3))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]



	tessellation_File = "%s/analysis/%s/%s/setup/" % (dir, galaxy, opt) + \
		"voronoi_2d_binning_output.txt"

	FWHM_gal = params.FWHM_gal/(1+z) # Adjust resolution in Angstrom

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

	dataCubeDirectory = get_dataCubeDirectory(galaxy)

	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	s = galaxy_data.shape

	rows_to_remove = range(params.discard)
	rows_to_remove.extend([s[1]-1-i for i in range(params.discard)])
	cols_to_remove = range(params.discard)
	cols_to_remove.extend([s[2]-1-i for i in range(params.discard)])

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

	bin_lin = np.nansum(galaxy_data[:,x[spaxels_in_bin],y[spaxels_in_bin]], axis=1)
	bin_lin_noise = np.nansum(galaxy_noise[:,x[spaxels_in_bin],
		y[spaxels_in_bin]]**2, axis=1)
	bin_lin_noise = np.sqrt(bin_lin_noise)
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = remove_anomalies(bin_lin, window=201, repeats=3, 
		lam=lam, set_range=params.set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
	bin_lin_noise = bin_lin_noise[cut]

	## smooth spectrum to fit with templates resolution
	if FWHM_gal < stellar_templates.FWHM_tem:
		sigma = stellar_templates.FWHM_dif/2.355/CDELT_spec # Change in px
		bin_lin = ndimage.gaussian_filter1d(bin_lin, sigma)
		bin_lin_noise = np.sqrt(ndimage.gaussian_filter1d(bin_lin_noise**2, 
			sigma))
	
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
	e_templates = get_emission_templates(params.gas, lamRange, 
		stellar_templates.logLam_template, FWHM_gal)

	if params.gas:
		templates = np.column_stack((stellar_templates.templates, e_templates.templates))
	else:
		templates = stellar_templates.templates
	component = [0]*len(stellar_templates.templatesToUse) + e_templates.component
	templatesToUse = np.append(stellar_templates.templatesToUse, 
		e_templates.templatesToUse)
	element = ['stellar'] + e_templates.element

	start = [[vel, sig]] * (max(component) + 1)
	moments = [params.stellar_moments] + [params.gas_moments] * max(component)


	goodPixels = determine_goodpixels(logLam_bin,stellar_templates.lamRange_template, 
		vel, z, gas=params.gas!=0)
## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
	noise = np.abs(noise)
	bin_log_sav = bin_log
	noise_sav = noise
	saveTo="%s/analysis/%s/%s/MC/bestfit/plots/%s.png" % (dir, galaxy, opt, str(bin))

	pp = ppxf(templates, bin_log, noise, velscale, start, 
			  goodpixels=goodPixels, moments=moments, degree=params.degree, vsyst=dv, 
			  component=component, lam=lambdaq, plot=not params.quiet, 
			  quiet=params.quiet, save=saveTo, mdegree=params.mdegree)
## ----------===============================================---------
## ----------================= The MC part =================---------
## ----------===============================================---------
	stellar_output = np.zeros((params.reps, params.stellar_moments))
	stellar_errors = np.zeros((params.reps, params.stellar_moments))
	if params.gas:
		gas_output = np.zeros((params.gas, params.reps, params.gas_moments))
		gas_errors = np.zeros((params.gas, params.reps, params.gas_moments))

	for rep in range(params.reps):
		# print rep
		random = np.random.randn(len(noise))
		add_noise = random*np.abs(noise)
		bin_log = pp.bestfit + add_noise
	
		ppMC = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodPixels, moments=moments, degree=params.degree, vsyst=dv, 
			lam=lambdaq, plot=not params.quiet, quiet=params.quiet, bias=0.1, 
			component=component, mdegree=params.mdegree)

		stellar_output[rep,:] = ppMC.sol[0][0:params.stellar_moments]
		stellar_errors[rep,:] = ppMC.error[0][0:params.stellar_moments]
		for g in range(len(element)-1):
			gas_output[g,rep,:] = ppMC.sol[g+1][0:params.gas_moments]
			gas_errors[g,rep,:] = ppMC.error[g+1][0:params.gas_moments]
## ----------============ Write ouputs to file =============---------
	saveAll(galaxy, bin, pp, lambdaq, stellar_output, stellar_errors, bin_log_sav, 
		noise_sav, element, templatesToUse, gas_output=gas_output, 
		gas_errors=gas_errors, opt=opt)
	
##############################################################################


# Use of plot_results.py


if __name__ == '__main__':
	errors2(5,'kin',29) if len(sys.argv)<4 else errors2()




