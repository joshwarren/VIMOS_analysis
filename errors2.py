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
from ppxf import ppxf, create_plot

c = 299792.458 # speed of light in km/s


#-----------------------------------------------------------------------------
class set_params(object):
	def __init__(self,
		opt 			=	'kin', 
		reps 			= 	10,
		quiet	 		= 	True, 
		gas 			= 	1, 
		set_range 		= 	np.array([4000,10000]), 
		stellar_moments = 	2, 
		gas_moments 	=	2, 
		narrow_broad 	=	False, 
		lines 			= 	'all', 
		produce_plot 	= 	True, 
		temp_mismatch 	= 	True,
		start 			= 	None,
		discard			= 	0,
		library 		=	'Miles',
		use_all_temp 	= 	False
		):
		if reps > 0 and reps < 4: raise ValueError('Reps must be either 0 or > 3')
		self.quiet = quiet # True
		self.gas = gas # 0   No gas emission lines
					# 1   Probe ionised gas
					# 2   Seperate gases heated by shocks (OIII and NI) and by 
					#		SF gas (Balmer lines)
					# 3   All gas seperate.
		self.reps = reps ## number of monte carlo reps per bin.
		self.FWHM_gal = 3 # Py3D reduces FHWM from 2.5 VIMOS documentation
		self.set_range = set_range # [2000,7410]
		self.stellar_moments = stellar_moments # 2
		self.gas_moments = gas_moments # 2
		self.narrow_broad = narrow_broad # False; Find 2 components to each gas 
										 # 	line
		self.lines = lines # 'all'; list of which lines to use
		self.produce_plot = produce_plot # True
		self.temp_mismatch = temp_mismatch # False
		self.start = start # None => vel and sigma in VIMOS galaxies.txt files
		self.discard = discard
		self.opt = opt
		self.library = library
		self.use_all_temp = use_all_temp

	@property
	def opt(self):
		return self._opt
	@opt.setter
	def opt(self, opt):
		self._opt = opt
		if 'kin' in opt:
			self.degree = 4  # order of addative Legendre polynomial used to 
							# correct the template continuum shape during the fit
			self.mdegree = 0
			self.gas = 0
		elif 'pop' in opt:
			self.degree = -1
			self.mdegree = 10

	@property
	def lines(self):
		if self.gas == 0:
			return []
		else:
			return self._lines
	@lines.setter
	def lines(self, value):
		if value == 'all' or value =='All' or value =='ALL':
			self._lines = ['Hdelta', 'Hgamma', 'Hbeta', '[OIII]5007d', '[NI]d']
		elif value == 'none' or value == 'None' or value is None:
			self._lines = []
			self.gas=0
		else:
			self._lines = list(value)
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

	def __init__(self, FWHM_gal, library='Miles'):
		self.library = library
		if self.library == 'Miles':
			self.FWHM_tem = 2.5 # Miles library has FWHM of 2.5A.
		elif self.library == 'vazdekis':
			self.FWHM_tem = 1.8
		self.FWHM_gal = FWHM_gal
		self.FWHM_dif = np.sqrt(abs(FWHM_gal**2 - self.FWHM_tem**2))

	def get_templates(self, galaxy, velscale=None, use_all_temp=False):
		import glob # for searching for files
		if self.library == 'vazdekis':
			templateFiles = glob.glob(
				'%s/libraries/python/ppxf/spectra/Rbi1.30z*.fits' % (cc.home_dir))

			f = fits.open(templateFiles[0])
			v2 = f[0].data
			self.wav = np.arange(f[0].header['NAXIS1'])*f[0].header['CDELT1'] + \
				f[0].header['CRVAL1']

			# Using same keywords as fits headers
			CRVAL_temp = f[0].header['CRVAL1']		# starting wavelength
			NAXIS_temp = f[0].header['NAXIS1']   # Number of entries
			# wavelength increments (resolution?)
			CDELT_temp = f[0].header['CDELT1']
			f.close()

			self.lamRange_template = CRVAL_temp + np.array(
				[0, CDELT_temp*(NAXIS_temp-1)])

			log_temp_template, self.logLam_template, self.velscale = \
				util.log_rebin(self.lamRange_template, self.wav, velscale=velscale)
			
			self.ntemp = len(templateFiles)
			self.templatesToUse = np.arange(self.ntemp)
			
			self.templates = np.zeros((len(log_temp_template), self.ntemp))
			self.lin_templates = np.zeros((f[0].header['NAXIS1'], self.ntemp))

			## Reading the contents of the files into the array templates. 
			## Including rebinning them.
			for i in range(self.ntemp):
				f = fits.open(templateFiles[i])
				self.lin_templates[:,i] = f[0].data
				f.close()
				if self.FWHM_tem < self.FWHM_gal:
					sigma = self.FWHM_dif/2.355/CDELT_temp # Sigma difference (px)
					conv_temp = ndimage.gaussian_filter1d(self.lin_templates[:,i],
						sigma)
				else:
					conv_temp = self.lin_templates[:,i]
				## Rebinning templates logarthmically
				log_temp_template, self.logLam_template, _ = util.log_rebin(
					self.lamRange_template, conv_temp, velscale=self.velscale)
				self.templates[:,i] = log_temp_template

		elif self.library=='Miles': # Miles stars
			# Finding the template files
			# There is some issue with the memory structure of the university macs 
			# (HFS+), meaning these templates can only be loaded once if located on 
			# the home directory, but more if on the Data partition...
			if cc.device != 'uni':
				templateFiles = glob.glob(
					'%s/models/miles_library/m0[0-9][0-9][0-9]V' % (cc.home_dir))
			else:
				templateFiles = glob.glob('%s/Data/' % (cc.base_dir) +
					'idl_libraries/ppxf/MILES_library/m0[0-9][0-9][0-9]V')

			# self.wav is wavelength, v2 is spectrum
			self.wav, v2 = np.loadtxt(templateFiles[0], unpack='True')

			# Using same keywords as fits headers
			CRVAL_temp = self.wav[0]		# starting wavelength
			NAXIS_temp = np.shape(v2)[0]   # Number of entries
			# wavelength increments (resolution?)
			CDELT_temp = (self.wav[NAXIS_temp-1]-self.wav[0])/(NAXIS_temp-1)

			self.lamRange_template = CRVAL_temp + np.array(
				[0, CDELT_temp*(NAXIS_temp-1)])

			log_temp_template, self.logLam_template, self.velscale = \
				util.log_rebin(self.lamRange_template, self.wav, velscale=velscale)
			
			if use_all_temp or galaxy is None:
				self.ntemp = len(templateFiles)
				self.templatesToUse = np.arange(self.ntemp)
			else:
				self.templatesToUse = use_templates(galaxy, cc.device=='glamdring')
				self.ntemp = len(self.templatesToUse)
			self.templates = np.zeros((len(log_temp_template), self.ntemp))
			self.lin_templates = np.zeros((len(v2), self.ntemp))

			## Reading the contents of the files into the array templates. 
			## Including rebinning them.
			for i in range(self.ntemp):
				if use_all_temp:
					_, self.lin_templates[:,i] = np.loadtxt(templateFiles[i], 
						unpack='True')
				else:
					_, self.lin_templates[:,i] = np.loadtxt(
						templateFiles[self.templatesToUse[i]], unpack='True')
				if self.FWHM_tem < self.FWHM_gal:
					sigma = self.FWHM_dif/2.355/CDELT_temp # Sigma difference (px)
					conv_temp = ndimage.gaussian_filter1d(self.lin_templates[:,i],
						sigma)
				else:
					conv_temp = self.lin_templates[:,i]
				## Rebinning templates logarthmically
				log_temp_template, self.logLam_template, _ = util.log_rebin(
					self.lamRange_template, conv_temp, velscale=self.velscale)
				self.templates[:,i] = log_temp_template
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
class get_emission_templates(object):
	def __init__(self, gas, lamRange, logLam_template, FWHM_gal, quiet=True,
		lines='all', narrow_broad=False, goodWav=None):
		if narrow_broad:
			raise ValueError('Two components (narrow and broad line componets' +
				' are not currently set up on the version (See git' +
				' errors2_muse.py)')

		if lines == 'all' or lines =='All' or lines =='ALL':
			lines = ['Hdelta', 'Hgamma', 'Hbeta', 'Halpha', '[OII]3726', 
				'[OII]3729', '[SII]6716', '[SII]6731', '[OIII]5007d', '[NI]d', 
				'[OI]6300d', '[NII]6583d']
		self.element = []
		self.component = []
		self.templatesToUse = []
		self.templates = []

		# This is not particularly robust code
		if goodWav is not None:
			w = np.where(np.diff(goodWav) > goodWav[-1]-goodWav[-2])[0]
			mask = np.array(zip(goodWav[w], goodWav[w+1]))
		else: 
			mask = np.array([[lamRange[0],lamRange[0]]])
		if len(mask) == 0:
			mask = np.array([[lamRange[0],lamRange[0]]])
		## ----------============ All lines together ===============---------
		if gas == 1:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			nTemp = 0
			for i in range(len(line_name)):
				if np.sum(mask[:,0] - line_wav[i] > 0) == np.sum(
					mask[:,1] - line_wav[i] > 0) and line_name[i] in lines:
					self.templatesToUse = np.append(self.templatesToUse, 
						line_name[i])
					self.templates.append(emission_lines[:,i])
					nTemp += 1
			self.component = self.component + [1]*nTemp
			self.element.append('gas')
		## ----------=============== SF and shocks lines ==============---------
		if gas == 2:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			for i in range(len(line_name)):
				if np.sum(mask[:,0] - line_wav[i] > 0) == np.sum(
					mask[:,1] - line_wav[i] > 0) and line_name[i] in lines:
					if 'H' in line_name[i]:
						self.templatesToUse = np.append(self.templatesToUse, 
							line_name[i])
						self.templates.append(emission_lines[:,i])
						self.component = self.component + [1]
						self.element.append['SF']
					else:
						self.templatesToUse = np.append(self.templatesToUse, 
							line_name[i])
						self.templates.append(emission_lines[:,i])
						self.component = self.component + [2] 
						self.element.append['Shocks']
		## ----------=========== All lines inderpendantly ==============---------
		if gas == 3:
			emission_lines, line_name, line_wav = util.emission_lines(
				logLam_template, lamRange, FWHM_gal, quiet=quiet)

			aph_lin = np.sort(line_name)

			for i in range(len(line_name)):
				if np.sum(mask[:,0] - line_wav[i] > 0) == np.sum(
					mask[:,1] - line_wav[i] > 0) and line_name[i] in lines:
					self.templatesToUse = np.append(self.templatesToUse,
						line_name[i])

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
def determine_goodpixels(logLam, lamRangeTemp, vel, z, lines='all',	invert=False):
	"""
	warrenj 20150905 Copied from ppxf_determine_goodPixels.pro
	
	ppxf_DETERMINE_GOODPIXELS: Example routine to generate the vector of
		goodPixels to be used as input keyword for the routine
		ppxf. This is useful to mask gas emission lines or atmospheric
		absorptions. It can be trivially adapted to mask different
		lines. 
	
	INPUT PARAMETERS:
	- LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
		of each pixel of the log rebinned *galaxy* spectrum.
	- LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the
		minimum and maximum wavelength in Angstrom in the stellar
		*template* used in ppxf. 
	- VEL: Estimate of the galaxy velocity in km/s.
	
	V1.0: Michele Cappellari, Leiden, 9 September 2005
	V1.01: Made a separate routine and included additional common
	  emission lines. MC, Oxford 12 January 2012
	V1.02: Included more lines. MC, Oxford, 7 Januray 2014
	20150617 warrenj Added Telluric lines (tell) at 5199 (is a blended sky
	line)
	20160816 warrenj Added NI doublet 5199.36 and 5201.86
	20170907 warrenj commented out lines at 5528A and 5535A as I cannot find ID
	"""

	# dv = lines*0+800d # width/2 of masked gas emission region in km/s
	dv = 400 # width/2 of masked gas emission region in km/s
	# flag = bytearray([0]*len(logLam))
	flag = logLam < 0

	# Marks telluric line
	tell = 5199
	flag |= (logLam > np.log(tell) - z - dv/c) \
		& (logLam < np.log(tell) - z + dv/c) 

	# if mask_sky:
	# 	flag |= (logLam > np.log(5450)) & (logLam < np.log(5550)) 
	# 	flag |= (logLam > np.log(5780)) & (logLam < np.log(5950)) 
	# 	flag |= (logLam > np.log(6120)) & (logLam < np.log(6520)) 
	# 	flag |= (logLam > np.log(6720)) & (logLam < np.log(7000)) 
	# 	flag |= (logLam > np.log(7140)) & (logLam < np.log(lamRangeTemp[1])) 

	if lines != 'all':
		#                 -----[OII]-----   Hdelta    Hgamma   Hbeta;
		wav = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 
		#    -----[OIII]----- -----[NI]-----      -----[OI]-----   
			4958.92, 5006.84, 5199.36, 5201.86,  6300.30, 6363.67,
		#    -----[NII]------  Halpha   -----[SII]-----   
			6548.03, 6583.41,  6562.80, 6716.47, 6730.85]) #5528, 5535

		if len(lines) != 0 and lines is not None and lines != 'none':
			# Mask lines in list called 'lines'
			line_names = ['[OII]3729','[OII]3729','Hdelta', 'Hgamma', 'Hbeta', 
				'[OIII]5007d', '[OIII]5007d', '[NI]d', '[NI]d', '[OI]6300d', 
				'[OI]6300d', '[NII]6583d', '[NII]6583d', 'Halpha', '[SII]6716', 
				'[SII]6731']
			wav = np.array([l for i, l in enumerate(wav) if 
				line_names[i] not in lines])
		else:
			pass # Mask all lines

		for j in range(len(wav)):
			flag |= (logLam > np.log(wav[j]) + (vel- dv)/c) \
				& (logLam < np.log(wav[j]) + (vel+ dv)/c)
	else:
		pass # Fit all lines

	flag |= logLam < np.log(lamRangeTemp[0]) + (vel + 900)/c # Mask edges of
	flag |= logLam > np.log(lamRangeTemp[1]) + (vel - 900)/c # stellar library

	if invert:
		flag[0:3] = 0 # Mask edge of data
		flag[-4:]= 0 # to remove edge effects
		return np.where(flag != 0)[0]
	else:
		flag[0:3] = 1 # Mask edge of data
		flag[-4:]= 1 # to remove edge effects
		return np.where(flag == 0)[0]
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def saveAll(galaxy, bin, pp, opt='kin'):
	# stellar MC results
	if cc.device == 'glamdring':
		dir = '%s/analysis/%s/%s/MC' % (cc.base_dir, galaxy, opt)
	else:
		dir = '%s/Data/vimos/analysis/%s/%s/MC' % (cc.base_dir, galaxy, opt)

	reps = pp.MCstellar_kin.shape[0]
	gas = len(pp.element) > 1

	def check_directory(d):
		if not os.path.exists(d):
			try:
				os.makedirs(d)
			except OSError: 
				# Another node has created directory in the time that it took for if 
				# statement above
				pass

	# stellar MC results
	check_directory("%s/stellar/errors" % (dir))
	bin_file = "%s/stellar/%s.dat" % (dir, str(bin))
	errors_file = "%s/stellar/errors/%s.dat" % (dir, str(bin))
	with open(bin_file, 'w') as f,  open(errors_file, 'w') as e:
		for i in range(reps):
			f.write("   ".join([str(s) for s in pp.MCstellar_kin[i,:]]) + '\n')
			e.write("   ".join([str(s) for s in pp.MCstellar_kin_err[i,:]]) + '\n')

	
	# gas MC kinematics results
	for d, gas_dir in enumerate([e for e in pp.element if e != 'stellar']):
		check_directory("%s/gas/%s/errors" % (dir, gas_dir))

		gas_file = "%s/gas/%s/%s.dat" % (dir, gas_dir, str(bin))
		gas_error_file = "%s/gas/%s/errors/%s.dat" % (dir, gas_dir, str(bin))

		with open(gas_file, 'w') as g, open(gas_error_file, 'w') as ger:
			for i in range(reps):
				if pp.MCgas_kin is not None:
					g.write("   ".join([str(s) for s in pp.MCgas_kin[d,i,:]]) + '\n')
				if pp.MCgas_kin_err is not None:
					ger.write("   ".join([str(s) for s in pp.MCgas_kin_err[d,i,:]]) + '\n')

	# gas MC uncertainty spectrum
	for d, gas_dir in enumerate([e for e in pp.templatesToUse if not e.isdigit()]):
		check_directory("%s/gas_uncert_spectrum/%s" %(dir, gas_dir))
		gas_uncert_file = "%s/gas_uncert_spectrum/%s/%s.dat" % (dir, gas_dir, str(bin))
		with open(gas_uncert_file, 'w') as u:
			for i in pp.MCgas_uncert_spec[d,:]:
				u.write(str(i) + '\n')

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
		for i in range(len(pp.galaxy)):
			inp.write(str(pp.galaxy[i]) + '\n')

	## save input noise
	check_directory("%s/noise_input" % (dir))
	input_file = "%s/noise_input/%s.dat" % (dir, str(bin))
   
	with open(input_file, 'w') as inp:
		for i in range(len(pp.noise)):
			inp.write(str(pp.noise[i]) + '\n')

	## save bestfit LOSVD output
	bestfit_file = "%s/%s.dat" % (dir, str(bin))
	with open(bestfit_file, 'w') as b:
		for i in range(len(pp.element)):
			b.write(pp.element[i])
			if gas:
				for j in pp.sol[i]:
					b.write("   " + str(j))
			else: # gas = 0 
				for j in pp.sol:
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
			w.write(str(pp.templatesToUse[i]) + "   " + str(pp.weights[i]) + '\n') 

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
				l.write(str(pp.templatesToUse[i]) + "   ")
				for k in range(pp.matrix.shape[0]):
					l.write(str(pp.matrix[k,j]) + "   ")
				l.write('\n')
	else:
		with open(matrix_file, 'w') as l:
			for i in range(pp.matrix.shape[1]):
				l.write(str(pp.templatesToUse[i]) + "   ")
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
		for i in range(len(pp.lam)):
			l.write(str(pp.lam[i]) + '\n')

	## save plot
	check_directory("%s/bestfit/plots" % (dir))
	plot_file="%s/bestfit/plots/%i.png" % (dir, bin)
	pp.fig.savefig(plot_file, bbox_inches="tight")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def apply_range(spec, lam=None, set_range=None, return_cuts=False):
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
def get_dataCubeDirectory(galaxy, radio_band=None):
	class mystring2(str):
		def __init__(self, s):
			str.__init__(s)
			self.RAoffset = 0 # offset in arcsec
			self.decoffset = 0
			self.band = ''
			self.default_scale = None


	class mystring(str):
		def __init__(self, s):
			str.__init__(s)
			self.radio = mystring2('')
			self.CO = mystring2('')
			self.xray = mystring2('')


	if cc.device == 'uni':
		dir = '%s/Data/vimos' % (cc.base_dir)
	elif cc.device == 'glamdring':
		dir = '%s' % (cc.base_dir)
	elif 'home' in cc.device:
		dir = '%s/Data/vimos' % (cc.base_dir)

	dataCubeDirectory = mystring('%s/cubes/%s.cube.combined.corr.fits' %  (dir, galaxy))
	dataCubeDirectory.CO = mystring2("%s/Data/alma/%s-mom0.fits" % (cc.base_dir, galaxy))

	# Using offsets file
	if cc.device == 'glamdring':
		offsets_file = '%s/offsets.txt' % (cc.base_dir)
	else:
		offsets_file = '%s/Data/offsets.txt' % (cc.base_dir)
	file_headings = np.genfromtxt(offsets_file, dtype=str, max_rows=1)
	galaxies, CO_RA, CO_dec = np.loadtxt(offsets_file, dtype=str, usecols=(0,9,10), 
		skiprows=2, unpack=True)

	if galaxy == 'eso443-g024':
		# dataCubeDirectory.xray = '%s/Data/Chandra/ESO443_full.fits' % (cc.base_dir)

		if radio_band is None or radio_band == 'L':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/ESO443.CBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'L band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_L')[0][0]

		elif radio_band == 'U':
			# Plot linear not log contours: very small extent
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/ESO443.UBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'U band (14.94 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_U')[0][0]
	elif galaxy == 'ic1459':
		pass
		# dataCubeDirectory.xray = '%s/Data/Chandra/IC1459_full.fits' % (cc.base_dir)
	elif galaxy == 'ic1531':
		# dataCubeDirectory.xray = '%s/Data/Chandra/IC1531_full.fits' % (cc.base_dir)

		if radio_band is None or radio_band == 'L':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/IC1531.LBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'L band (1.45 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_L')[0][0]
		elif radio_band == 'X':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/IC1531.XBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'X band (8.46 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_X')[0][0]
	elif galaxy == 'ic4296':
		# dataCubeDirectory.xray = '%s/Data/Chandra/IC4296_full.fits' % (cc.base_dir)

		if radio_band is None or radio_band == 'L':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/IC4296.LBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'L band (1.45 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_L')[0][0]
		elif radio_band == 'C':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/IC4296.CBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'C band (4.87 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
	elif galaxy == 'ngc0612':
		# dataCubeDirectory.xray = '%s/Data/Chandra/N612_full.fits' % (cc.base_dir)

		dataCubeDirectory.CO = mystring2("%s/Data/alma/ngc612-mom0.fits" % (cc.base_dir))
		if radio_band is None or radio_band == 'C' or radio_band == 'CI':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc0612/NGC612.I.13ARCSEC.CBAND.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band == 'CQ':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc0612/NGC612.Q.13ARCSEC.CBAND.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band == 'CU':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc0612/NGC612.U.13ARCSEC.CBAND.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
	elif galaxy == 'ngc1399':
		# dataCubeDirectory.xray = '%s/Data/Chandra/N1399_full.fits' % (cc.base_dir)
		if radio_band is None or radio_band == 'C' or radio_band == 'CI':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc1399/NGC1399.CBAND.ANB.I.1.25.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band == 'CQ':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc1399/NGC1399.CBAND.ANB.Q.1.25.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band == 'CU':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/ngc1399/NGC1399.CBAND.ANB.U.1.25.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
	elif galaxy == 'ngc3100':
		if radio_band is None or radio_band == 'L':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/Southern_RG_VLA/AD254.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'L band (1.45 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_L')[0][0]

		elif radio_band == 'C':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/Southern_RG_VLA/AD270.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'C band (4.86 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band == 'X':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/Southern_RG_VLA/AL508.fits' % (
				cc.base_dir))
			dataCubeDirectory.radio.band = 'X band (8.46 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_X')[0][0]
	elif galaxy == 'ngc3557':
		if radio_band is None or radio_band =='C':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/NGC3557.ABC.6CM.PBCOR.FITS')
			dataCubeDirectory.radio.band = 'C band (4.87 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band =='C2':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) + 
				'Southern_RG_VLA/NGC3557.CD.6CM.PBCOR.FITS')
			dataCubeDirectory.radio.band = 'C band (4.87 GHz)'
			dataCubeDirectory.radio_band.RAoffset = 0
			dataCubeDirectory.radio_band.decoffset = 0
	elif galaxy == 'ngc7075':
		if radio_band is None or radio_band == 'C':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/NGC7075BC.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'C band (4.87 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_C')[0][0]
		elif radio_band =='X':
		# Plot linear, not log contours: very small extent
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/NGC7075.XBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'X band (8.46 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_X')[0][0]
	elif galaxy == 'pks0718-34':
		if radio_band is None or radio_band == 'X':
			dataCubeDirectory.radio = mystring2('%s/Data/VLA/' % (cc.base_dir) +
				'Southern_RG_VLA/0718.XBAND.IMAGE.FITS')
			dataCubeDirectory.radio.band = 'X band (8.46 GHz)'
			col = np.where(file_headings=='VIMOS-VLA_X')[0][0]

	# Extraction CO offsets
	CO_RA[CO_RA=='-'] = 'nan'
	CO_dec[CO_dec=='-'] = 'nan'
	i_gal = np.where(galaxies == galaxy)[0][0]
	dataCubeDirectory.CO.RAoffset = np.float(CO_RA[i_gal].strip('*'))
	dataCubeDirectory.CO.decoffset = np.float(CO_dec[i_gal].strip('*'))

	# Extracting radio offsets (found by eye)
	if 'Data' in dataCubeDirectory.radio and radio_band != 'C2': # check is set
		col *= 2
		dataCubeDirectory.radio.RAoffset, dataCubeDirectory.radio.decoffset = \
			np.genfromtxt(offsets_file, unpack=True, usecols=(col-1, col), dtype=str,
			skip_header=2, missing_values='-', filling_values='nan')[:, i_gal]
		dataCubeDirectory.radio.RAoffset = float(
			dataCubeDirectory.radio.RAoffset.strip('*'))
		dataCubeDirectory.radio.decoffset = float(
			dataCubeDirectory.radio.decoffset.strip('*'))

	return dataCubeDirectory
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def errors2(i_gal=None, opt=None, bin=None):
	if i_gal is None: i_gal=int(sys.argv[1])
	if opt is None: opt=sys.argv[2]
	if bin is None: bin=int(sys.argv[3])
	## ----------===============================================---------
	## ----------============= Input parameters  ===============---------
	## ----------===============================================---------
	params = set_params(opt=opt)
	
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[i_gal]

	if cc.device == 'glamdring':
		dir = cc.base_dir
	else:
		dir = '%s/Data/vimos' % (cc.base_dir)
	tessellation_File = "%s/analysis/%s/%s/setup/" % (dir, galaxy, opt) + \
		"voronoi_2d_binning_output.txt"
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

	bin_lin = np.nansum(galaxy_data[:,x[spaxels_in_bin],y[spaxels_in_bin]], axis=1)/\
		n_spaxels_in_bin
	bin_lin_noise = np.nansum(galaxy_noise[:,x[spaxels_in_bin],
		y[spaxels_in_bin]]**2, axis=1)
	bin_lin_noise = np.sqrt(bin_lin_noise)/n_spaxels_in_bin
	## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = apply_range(bin_lin, lam=lam, set_range=params.set_range, 
		return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])
	bin_lin_noise = bin_lin_noise[cut]


	pp = run_ppxf(galaxy, bin_lin, bin_lin_noise, lamRange, CDELT_spec, params)

	## ----------============ Write ouputs to file =============---------
	saveAll(galaxy, bin, pp, opt=opt)
	



## ----------===============================================---------
## ----------=============== Run analysis  =================---------
## ----------===============================================---------
class run_ppxf(ppxf):
	def __init__(self, galaxy_name, bin_lin, bin_lin_noise, lamRange, CDELT, params, 
		z=0.):

		self.galaxy_name = galaxy_name
		self.bin_lin = bin_lin
		self.bin_lin_noise = bin_lin_noise
		self.lamRange = lamRange
		self.CDELT = CDELT
		self.params = params


		if galaxy_name is not None:
			if cc.device == 'glamdring':
				data_file = "%s/analysis/galaxies.txt" % (cc.base_dir)
			else:
				data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
			# different data types need to be read separetly
			z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
				usecols=(1,2,3))
			galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
			i_gal = np.where(galaxy_gals==galaxy_name)[0][0]
			self.vel = vel_gals[i_gal]
			self.sig = sig_gals[i_gal]
			self.z = z_gals[i_gal]
		else: 
			self.vel = 150.
			self.sig = 100.
			self.z = z

		self.lamRange = self.lamRange/(1 + self.z)
		self.FWHM_gal = self.params.FWHM_gal/(1 + self.z) # Adjust resolution in Angstrom

		self.rebin()
		self.load_stellar_templates()

		if not self.params.temp_mismatch or self.params.gas==0:
			self.load_emission_templates()
			self.run()
		else:
			from copy import copy
			# Fit stellar component only.
			params_sav = copy(self.params)
			self.params.gas = 0
			self.params.produce_plot = False
			self.load_emission_templates()

			save_bin_log = np.array(self.bin_log)
			save_bin_log_noise = np.array(self.bin_log_noise)
			save_lamRange = np.array(self.lamRange)

			self.run() # stars

			MCstellar_kin = np.array(self.MCstellar_kin)
			MCstellar_kin_err = np.array(self.MCstellar_kin_err)


			# Find [OIII] kinematics and amplitude, enforce stellar kinematics
			self.params.stellar_moments *= -1
			self.params.start = [self.sol[0].tolist(), 
				[None]*self.params.gas_moments]
			self.params.gas = params_sav.gas
			self.params.lines = ['[OIII]5007d']
			self.load_emission_templates()

			self.bin_log = np.copy(save_bin_log)
			self.bin_lin_noise = np.copy(save_bin_log_noise)
			self.lamRange = np.copy(save_lamRange)

			self.run() # stars and [OIII]

			MCgas_kin = np.array(self.MCgas_kin)
			MCgas_kin_err = np.array(self.MCgas_kin_err)
			OIII_uncert_spec = np.array(self.MCgas_uncert_spec)

			# Find all other gas amplitudes enforcing [OIII] and stellar kinematics
			self.params.gas_moments *= -1
			self.params.start = [i.tolist() for i in self.sol]
			self.params.lines = params_sav.lines
			self.load_emission_templates()

			self.bin_log = np.copy(save_bin_log)
			self.bin_lin_noise = np.copy(save_bin_log_noise)
			self.lamRange = np.copy(save_lamRange)

			self.run() # all in

			self.MCstellar_kin = MCstellar_kin
			self.MCstellar_kin_err = MCstellar_kin_err
			self.MCgas_kin = MCgas_kin
			self.MCgas_kin_err = MCgas_kin_err
			self.MCgas_uncert_spec[
				self.templatesToUse[self.component!=0]=='[OIII]5007d', :] = \
				OIII_uncert_spec.flatten()
	
			if params_sav.produce_plot:
				self.fig, self.ax = create_plot(self).produce


	def rebin(self):
		self.stellar_templates = get_stellar_templates(self.FWHM_gal, 
			library=self.params.library)
		
		## smooth spectrum to fit with templates resolution
		if self.FWHM_gal < self.stellar_templates.FWHM_tem:
			sigma = self.stellar_templates.FWHM_dif/2.355/self.CDELT # Change in px
			self.bin_lin = ndimage.gaussian_filter1d(self.bin_lin, sigma)
			self.bin_lin_noise = np.sqrt(ndimage.gaussian_filter1d(
				self.bin_lin_noise**2, sigma))
		
		## rebin spectrum logarthmically
		self.bin_log, self.logLam_bin, self.velscale = util.log_rebin(self.lamRange, 
			self.bin_lin)
		bin_log_noise, logLam_bin, _ = util.log_rebin(self.lamRange, 
			self.bin_lin_noise**2)
		self.bin_log_noise = np.sqrt(bin_log_noise)
		self.lambdaq = np.exp(self.logLam_bin)

	## ----------============= Stellar templates ===============---------
	def load_stellar_templates(self):
		if self.params.use_all_temp is None and self.params.gas == 0:
			raise ValueError('No templates to fit... gas or stellar')
		
		self.stellar_templates.get_templates(self.galaxy_name, self.velscale, 
			use_all_temp=self.params.use_all_temp)
		# self.velscale = stellar_templates.velscale
		self.dv = (self.stellar_templates.logLam_template[0] - 
			self.logLam_bin[0])*c # km/s	
	## ----------=============== Emission lines ================---------
	def set_goodPix(self):
		self.goodPixels = determine_goodpixels(self.logLam_bin,
			self.stellar_templates.lamRange_template, self.vel, self.z, 
			lines=self.params.lines, invert=self.params.use_all_temp is None)
		self.goodPixels = np.array([g for g in self.goodPixels if 
			(~np.isnan(self.bin_log[g]))])
	def load_emission_templates(self):
		self.set_goodPix()
		self.e_templates = get_emission_templates(self.params.gas, 
			self.lamRange, self.stellar_templates.logLam_template, 
			self.FWHM_gal, goodWav=self.lambdaq[self.goodPixels], 
			narrow_broad=self.params.narrow_broad, lines=self.params.lines)

		if self.params.use_all_temp is not None:
			if self.params.gas:
				self.templates = np.column_stack((self.stellar_templates.templates, 
					self.e_templates.templates))
			else:
				self.templates = self.stellar_templates.templates
			self.component = [0]*len(self.stellar_templates.templatesToUse
				) + self.e_templates.component
			self.templatesToUse = np.append(self.stellar_templates.templatesToUse.astype(str), 
				self.e_templates.templatesToUse)
			self.element = ['stellar'] + self.e_templates.element
		else:
			self.templates = self.e_templates.templates
			self.component = [i - 1 for i in self.e_templates.component]
			self.templatesToUse = self.e_templates.templatesToUse
			self.element = self.e_templates.element

	def run(self):
		if self.params.start is None:
			if not self.params.narrow_broad:
				start = [[self.vel, self.sig]] * (max(self.component) + 1)
			else:
				start = [[self.vel, self.sig]] * (max(self.component)/2 + 1)
				# Start narrow line component at half of broad line
				start.extend([[self.vel, self.sig/2]] * (max(self.component)/2))
		else:
			start = self.params.start
			for i, e in enumerate(self.element):
				if start[i][0] is None:
					start[i][0] = self.vel
				if start[i][1] is None:
					if 'n_' not in e:
						start[i][1] = self.sig
					else:
						start[i][1] = self.sig/2.
		
		moments = [self.params.stellar_moments] + [self.params.gas_moments] * \
			max(self.component)
	## ----------============== The bestfit part ===============---------
		ppxf.__init__(self, self.templates, self.bin_log, self.bin_log_noise, 
			self.velscale, start, goodpixels=self.goodPixels, 
			mdegree=self.params.mdegree, moments=moments, 
			degree=self.params.degree, vsyst=self.dv, component=self.component, 
			lam=self.lambdaq, plot=not self.params.quiet, quiet=self.params.quiet, 
			produce_plot=False)
		if self.params.gas == 0: 
			self.sol = [self.sol]
			self.error = [self.error]

	## ----------===============================================---------
	## ----------================= The MC part =================---------
	## ----------===============================================---------
		if self.params.use_all_temp is not None:
			self.MCstellar_kin = np.zeros((self.params.reps, 
				abs(self.params.stellar_moments)))
			self.MCstellar_kin_err = np.zeros((self.params.reps, 
				abs(self.params.stellar_moments)))
		self.MCbestfit_uncert = np.zeros((3, len(self.galaxy)))
		MCbestfit_mean = np.zeros(len(self.galaxy))

		if self.params.gas:
			self.MCgas_kin = np.zeros((max(self.component), self.params.reps,
				abs(self.params.gas_moments)))
			self.MCgas_kin_err = np.zeros((max(self.component), self.params.reps, 
				abs(self.params.gas_moments)))

			n_lines = len(self.e_templates.templatesToUse)
			# self.MCgas_weights = np.zeros((n_lines, self.params.reps))

			self.MCgas_uncert_spec = np.zeros((n_lines, 3, len(self.galaxy)))
			MCgas_mean_spec = np.zeros((n_lines, len(self.galaxy)))

			if self.params.reps == 0:
				self.MCgas_uncert_spec = np.zeros((n_lines, len(self.galaxy)))
		else:
			self.MCgas_kin = None
			self.MCgas_kin_err = None
			# self.MCgas_weights = None



		for rep in range(self.params.reps):
			random = np.random.randn(len(self.bin_log_noise))
			add_noise = random*np.abs(self.bin_log_noise)
			self.bin_log = self.bestfit + add_noise

			ppMC = ppxf(self.templates, self.bin_log, self.bin_log_noise, 
				self.velscale, start, goodpixels=self.goodPixels, moments=moments, 
				degree=self.params.degree, vsyst=self.dv, lam=self.lambdaq, 
				plot=not self.params.quiet, quiet=self.params.quiet, bias=0.1, 
				component=self.component, mdegree=self.params.mdegree)

			if self.params.gas == 0: 
				ppMC.sol = [ppMC.sol]
				ppMC.error = [ppMC.error]

			if self.params.use_all_temp is not None:
				self.MCstellar_kin[rep,:] = ppMC.sol[0][
					0:abs(self.params.stellar_moments)]
				self.MCstellar_kin_err[rep,:] = ppMC.error[0][0:
					abs(self.params.stellar_moments)]

			# Find uncertainty in bestfit
			new_mean_bestfit_spec = ((rep + 1) * MCbestfit_mean + ppMC.bestfit)/(
				rep + 2)
			if rep < 3:
				# Save spec until 3 reps have been completed
				self.MCbestfit_uncert[rep, :] = ppMC.bestfit
			else:
				# Finding sigma_N from x_N, mean_N, mean_(N-1) and sigma_(N-1)
				# NB: rep = N-1 due to being zero-based
				self.MCbestfit_uncert = np.sqrt(((ppMC.bestfit - 
					new_mean_bestfit_spec) * (ppMC.bestfit - MCbestfit_mean) + 
					rep * self.MCbestfit_uncert**2)/(rep + 1))
			MCbestfit_mean = np.array(new_mean_bestfit_spec)
			if rep == 2:
				# Calc std at 3rd rep
				self.MCbestfit_uncert = np.std(self.MCbestfit_uncert, axis=0)

			# Gas kinematics from all reps
			for g in range(len(self.element) - 1):
				self.MCgas_kin[g,rep,:] = ppMC.sol[g+1][0:abs(self.params.gas_moments)]
				self.MCgas_kin_err[g,rep,:] = ppMC.error[g+1][
					0:abs(self.params.gas_moments)]

			# Find uncertainty in fitted emission lines
			for i, n in enumerate(self.e_templates.templatesToUse):
				e_line_spec = ppMC.matrix[:, -n_lines + i] * ppMC.weights[
					-n_lines + i]
				new_mean_gas_spec = ((rep + 1) * MCgas_mean_spec[i, :] + 
					e_line_spec)/(rep + 2)

				if rep < 3 or self.params.reps == 0:
					self.MCgas_uncert_spec[i, rep, :] = e_line_spec
				else:
					self.MCgas_uncert_spec[i,:] = np.sqrt(((e_line_spec - 
						new_mean_gas_spec) * (e_line_spec - MCgas_mean_spec[i,:]) + 
						rep * self.MCgas_uncert_spec[i,:]**2)/(rep + 1))
				MCgas_mean_spec[i, :] = np.array(new_mean_gas_spec)

			if rep == 2 and self.params.gas != 0:
				self.MCgas_uncert_spec = np.std(self.MCgas_uncert_spec, axis=1)

		if self.params.produce_plot:
			self.fig, self.ax = create_plot(self).produce
##############################################################################


# Use of plot_results.py


if __name__ == '__main__':
	errors2(5,'kin',29) if len(sys.argv)<4 else errors2()




