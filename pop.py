## ==================================================================
## 		Stellar population
## ==================================================================
## A method to impliment MCMC routines to find the bestfit model to the 
## self.bin given.
## 
##
#######################################################################
# Keywords:
# ab_lines: 			Dict of absorption line strength. NB: entry can be 
#						array to fit more than one self.bin simulataniously.
# uncerts:				Dict of uncertainty in the absorption line strenghts.
# interp:		None 	Dict of interpolated models. If None, models are read 
#						and interpolated everytime routine is called.
# grid_length	40		Int of the length of a side of the cube of interpolated 
#						models. NB: This is not used if interp is given.
#######################################################################
import numpy as np
import sys
import os
from glob import glob
from astropy.io import fits
from scipy.interpolate import LinearNDInterpolator
import emcee
from checkcomp import checkcomp
cc = checkcomp()
if cc.device == -1:
	cc = checkcomp(override='glamdring')
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
from absorption import absorption
from scipy.ndimage.filters import gaussian_filter1d
from tools import moving_weighted_average

# from prefig import Prefig 
# Prefig()

c = 299792.458

def get_vin_dir(instrument):
	if instrument=='vimos':
		if cc.device == 'glamdring': 
			vin_dir = '%s/analysis' % (cc.base_dir)
		else: 
			vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	elif instrument=='muse':
		if cc.device == 'glamdring': 
			vin_dir = '%s/analysis_muse' % (cc.base_dir)
		else: 
			vin_dir = '%s/Data/muse/analysis' % (cc.base_dir)
	return vin_dir

def get_absorption(lines, pp=None, galaxy=None, bin=None, opt=None, 
	instrument='vimos', sigma=None, res=None):
	if (pp is None and galaxy is None) or (pp is not None and galaxy is not None):
		raise ValueError('Either pp or galaxy must be supplied: one or tother')
	if galaxy is not None and (bin is None or opt is None):
		raise ValueError('If galaxy is supplied, so to must bin and opt')

	if pp is None:
		vin_dir = get_vin_dir(instrument)
		vin_dir_gasMC = "%s/%s/%s/MC" % (vin_dir, galaxy, opt)
	
		if cc.device == 'glamdring':
			data_file = "%s/analysis/galaxies.txt" % (cc.base_dir)
		else:
			data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
		# different data types need to be read separetly
		z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		i_gal = np.where(galaxy_gals==galaxy)[0][0]
		z = z_gals[i_gal]

		vel = np.loadtxt("%s/%d.dat" % (vin_dir_gasMC, bin), unpack=True, 
			usecols=(1,))[0]
		lam = np.loadtxt("%s/lambda/%d.dat" % (vin_dir_gasMC, bin)) * (
			1 + z)/(1 + z + (vel/c))
		spectrum = np.loadtxt("%s/input/%d.dat" % (vin_dir_gasMC, bin))
		matrix = np.loadtxt("%s/bestfit/matrix/%d.dat" % (vin_dir_gasMC, bin), 
			dtype=str)
		temp_weights =  np.loadtxt("%s/temp_weights/%d.dat" % (vin_dir_gasMC, 
			bin), unpack=True, usecols=(1,))
		noise = np.loadtxt("%s/noise_input/%d.dat" % (vin_dir_gasMC, bin))
		bestfit = np.loadtxt("%s/bestfit/%d.dat" %(vin_dir_gasMC, bin))
		mpweight = np.loadtxt("%s/mpweights/%d.dat" %(vin_dir_gasMC, bin))
		e_lines = np.array([not s.isdigit() for s in matrix[:,0]])
		e_line_spec = matrix[e_lines,1:].astype(float)
		stellar_temps = matrix[~e_lines,0]
		e_temps = matrix[e_lines, 0]
		library = 'Miles'
	else:
		lam = pp.lam*(1+pp.z)/(1+pp.z+(pp.sol[0][0]/c))
		spectrum = pp.galaxy
		matrix = pp.matrix.T.astype(str)
		temp_weights = pp.weights
		noise = pp.noise 
		bestfit = pp.bestfit
		mpweight = pp.mpolyweights

		e_lines = pp.component != 0
		e_line_spec =  matrix[e_lines,:].astype(float)
		stellar_temps = pp.templatesToUse[~e_lines]
		e_temps = pp.templatesToUse[e_lines]

		library = pp.params.library

	e_line_spec = np.einsum('ij,i->ij',e_line_spec,temp_weights[e_lines])

	# Find and smooth residuals
	residuals = spectrum - bestfit
	_, residuals, _ = moving_weighted_average(lam, residuals, step_size=3., 
		interp=True)
	noise = np.sqrt(residuals**2 + noise**2)
	# Only use line if formally detected.
	if '[OIII]5007d' in e_temps:
		if np.max(e_line_spec[e_temps=='[OIII]5007d'])/np.median(
			noise[(np.log(lam) > np.log(5007) - 300/c) * 
			(np.log(lam) < np.log(5007) + 300/c)]) < 4:

			e_line_spec[:,:] = 0
		else:
			if '[NI]d' in e_temps:
				if np.max(e_line_spec[e_temps=='[NI]d'])/np.median(
					noise[(np.log(lam) > np.log(5200) - 300/c) * 
					(np.log(lam) < np.log(5200) + 300/c)]) < 4:

					e_line_spec[e_temps=='[NI]d',:] = 0
				else:
					print '[NI] detected'

			if 'Hbeta' in e_temps:
				if np.max(e_line_spec[e_temps=='Hbeta'])/np.median(
					noise[(np.log(lam) > np.log(4861) - 300/c) * 
					(np.log(lam) < np.log(4861) + 300/c)]) < 3:

					e_line_spec[e_temps=='Hbeta',:] = 0
	else:
		e_line_spec[:,:] = 0

	e_line_spec = np.sum(e_line_spec, axis=0)
	continuum = spectrum - e_line_spec

	if library == 'Miles':
		if cc.getDevice() == 'uni':
			files = glob('%s/Data/idl_libraries/ppxf/MILES_library/' % (cc.base_dir) +
				'm0[0-9][0-9][0-9]V')
		else:
			files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
		temp_res = 2.5 # A (FWHM)
		wav = np.loadtxt(files[0], usecols=(0,), unpack=True)

		a = [min(np.where(wav>=lam[0])[0]), max(np.where(wav<=lam[-1])[0])]
		unconvolved_spectrum = np.zeros(a[1]-a[0])

		for i, n in enumerate(stellar_temps):
			template =  np.loadtxt(files[int(n)], usecols=(1,), unpack=True)
			unconvolved_spectrum += template[a[0]:a[1]]*temp_weights[i]
		unconvolved_spectrum *= np.polynomial.legendre.legval(np.linspace(-1,1,
				len(unconvolved_spectrum)), np.append(1, mpweight))
		unconvolved_lam =  wav[a[0]:a[1]]
		CDELT = unconvolved_lam[1] - unconvolved_lam[0]
	elif library == 'vazdekis':
		templateFiles = glob(
			'%s/libraries/python/ppxf/spectra/Rbi1.30z*.fits' % (cc.home_dir))
		temp_res = 1.8 # A (FWHM)

		f = fits.open(templateFiles[0])
		unconvolved_lam = np.arange(f[0].header['NAXIS1'])*f[0].header['CDELT1'] + \
			f[0].header['CRVAL1']
		unconvolved_spectrum = np.zeros((f[0].header['NAXIS1'], len(templateFiles)))

		for i in range(len(templateFiles)):
			f = fits.open(templateFiles[i])
			unconvolved_spectrum[:,i] = f[0].data
		unconvolved_spectrum = unconvolved_spectrum.dot(
			temp_weights[~e_lines])
		unconvolved_spectrum *= np.polynomial.legendre.legval(
			np.linspace(-1,1,len(unconvolved_spectrum)), 
			np.append(1, mpweight))
		CDELT = f[0].header['CDELT1']
		f.close()

	# Convolve to required velocity dispersion
	if sigma is not None:
		sig_pix = np.median(unconvolved_lam)*(sigma/c)/CDELT
		convolved = gaussian_filter1d(unconvolved_spectrum, sig_pix)
	else:
		# bestfit is already convolved to velocty dispersion
		convolved = bestfit - e_line_spec
		

	# Continuum need be brought require resolution
	if res is not None:
		if instrument == 'vimos':
			instr_res = 2.5 # A (FWHM)
		elif instrument == 'muse':
			instr_res = 2.3 # A (FWHM)
		elif instrument == 'sauron':
			instr_res = 4.2

		if instr_res > res:
			import warnings
			warnings.warn('get_absorption cannot increase the ' +
				"resolution of %s from %.2f A to %.2f A: That's " % (
				instrument, instr_res, res) + 'impossible you chump!')
		else:
			sig_pix = np.sqrt(res**2 - instr_res**2) / 2.355 \
				/ np.median(np.diff(lam))
			continuum = gaussian_filter1d(continuum, sig_pix)


	# unc_spec and conv_spec need to be the same resolution
	# unc_spec is at temp_res
	conv_res = max((instr_res, temp_res))
	if conv_res < temp_res:
		if sigma is None:
			sig_pix = np.sqrt(temp_res**2 - conv_res**2) / 2.355 \
				/ np.median(np.diff(lam))
		else:
			sig_pix = np.sqrt(temp_res**2 - conv_res**2) / 2.355 \
				/ np.median(np.diff(unconvolved_lam))
		convolved = gaussian_filter1d(convolved, sig_pix)
	else:
		sig_pix = np.sqrt(conv_res**2 - temp_res**2) / 2.355 \
			/ np.median(np.diff(unconvolved_lam))
		unconvolved_spectrum = gaussian_filter1d(unconvolved_spectrum, 
			sig_pix)

	ab_lines = {}
	uncerts = {}
	for l in lines:
		if sigma is None:
			ab, uncert = absorption(l, lam, continuum, noise=noise,
				unc_lam=unconvolved_lam, unc_spec=unconvolved_spectrum, 
				conv_spec=convolved)
		else:
			ab, uncert = absorption(l, lam, continuum, noise=noise,
				unc_lam=unconvolved_lam, unc_spec=unconvolved_spectrum, 
				conv_spec=convolved, conv_lam=unconvolved_lam)
		ab_lines[l], uncerts[l] = ab[0], uncert[0]
	return ab_lines, uncerts

class population(object):

	def __init__(self, pp=None, galaxy=None, opt='pop', ab_index=None, 
		ab_uncert=None, instrument='vimos', method='mean'):
		self.pp = pp
		self.galaxy = galaxy
		self.instrument = instrument
		self.opt = opt

		if self.instrument == 'vimos':
			self.lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 
				'Fe5015', 'Mg_b']
		elif self.instrument == 'muse':
			self.lines = ['H_beta', 'Fe5015', 'Mg_b', 'Fe5270', 'Fe5335', 
				'Fe5406', 'Fe5709', 'Fe5782', 'NaD', 'TiO1', 'TiO2']
		elif self.instrument == 'sauron':
			self.lines = ['H_beta', 'Fe5015', 'Mg_b']

		grid_length = 40

		if ab_index is not None:
			if ab_uncert is None:
				raise ValueError('Uncertainty values must be supplied')
			self.ab_lines = ab_index
			self.uncerts = ab_uncert
			self.lines = self.ab_lines.keys()
		else:
			if self.pp is None:

				try:
					self.bin=int(sys.argv[4])
					self.instrument = str(sys.argv[1])
					self.i_gal=int(sys.argv[2])
					self.opt=str(sys.argv[3])
				except IndexError:
					self.i_gal=int(sys.argv[1])
					self.opt=str(sys.argv[2])
					self.bin=int(sys.argv[3])

				if self.instrument == 'vimos':
					galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 
						'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 
						'pks0718-34', 'eso443-g024']
					# self.lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
					# 	'H_beta', 'Fe5015', 'Mg_b']
				elif self.instrument == 'muse':
					galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
					# self.lines = ['H_beta', 'Fe5015', 'Mg_b', 'Fe5270', 
					# 	'Fe5335', 'Fe5406', 'Fe5709', 'Fe5782', 'NaD', 
					# 	'TiO1', 'TiO2']
				self.galaxy = galaxies[self.i_gal]

				self.ab_lines, self.uncerts = get_absorption(self.lines, 
					galaxy=self.galaxy, bin=self.bin, opt=self.opt, 
					instrument=self.instrument, res=2.5)
			else:
				self.ab_lines, self.uncerts = get_absorption(self.lines, 
					pp=self.pp, instrument=self.instrument, res=2.5)

		s=[grid_length,grid_length,grid_length]

		models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)
		titles = np.loadtxt(models_dir, dtype=str)[0]

		age_in, metallicity_in, alpha_in = np.loadtxt(models_dir, usecols=(0,1,2), 
			unpack=True, skiprows=35)

		self._age = np.logspace(np.log10(min(age_in)), np.log10(max(age_in)), num=s[0])
		self._age[np.argmax(self._age)] = max(age_in)
		self._metallicity = np.linspace(min(metallicity_in), max(metallicity_in), num=s[1])
		self._alpha = np.linspace(min(alpha_in), max(alpha_in), num=s[2])

		age_p=np.array([[m]*s[2]*s[1] for m in self._age]).flatten()
		metallicity_p=np.array([[m]*s[2] for m in self._metallicity]*s[0]).flatten()
		alpha_p=np.array(list(self._alpha)*s[1]*s[0])

		models = {}
		self.interp = {}
		for i, l in enumerate(titles):
			if l in self.lines:
				models[l] = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
					skiprows=35)
				self.interp[l] = LinearNDInterpolator(
					np.array([age_in, metallicity_in, alpha_in]).transpose(), 
					models[l])

		ndim, nwalkers, nsteps = 3, 200, 500
		self.samples = np.zeros((nwalkers*(nsteps-50), ndim))

		sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob)
		pos = [np.array([13,0.2,-0.1]) +1e-3*np.random.randn(ndim) for i in 
			range(nwalkers)]
		sampler.run_mcmc(pos, nsteps)

		# Discard initial steps away from initial point
		self.samples[:,:] = sampler.chain[:, 50:, :].reshape((-1, ndim))

		if method == 'mean':
			self.age = np.nanmean(self.samples[:,0])
			self.unc_age = np.nanstd(self.samples[:,0])
			self.metallicity = np.nanmean(self.samples[:,1])
			self.unc_met = np.nanstd(self.samples[:,1])
			self.alpha = np.nanmean(self.samples[:,2])
			self.unc_alp = np.nanstd(self.samples[:,2])
		elif method == 'mostlikely':
			print min(self.samples[:,0]), max(self.samples[:,0])
			hist = np.histogram(self.samples[:,0], bins=40)
			x = (hist[1][0:-1]+hist[1][1:])/2
			hist = hist[0]
			self.age = x[np.argmax(hist)]
			self.unc_age = np.nanstd(self.samples[:,0])
			
			# from tools import fwhm
			# self.unc_age = fwhm(x, hist)

			# gt_fwhm = hist >= (np.max(hist)/2)
			# # Not quite definition of fwhm, but close enough 
			# self.unc_age = (np.max(x[gt_fwhm]) - np.min(x[gt_fwhm]))/2

			hist = np.histogram(self.samples[:,1], bins=40)
			x = (hist[1][0:-1]+hist[1][1:])/2
			hist = hist[0]
			self.metallicity = x[np.argmax(hist)]
			self.unc_met = np.nanstd(self.samples[:,1])
			# self.unc_met = fwhm(x, hist)

			# gt_fwhm = hist >= (np.max(hist)/2)
			# self.unc_met = (np.max(x[gt_fwhm]) - np.min(x[gt_fwhm]))/2


			hist = np.histogram(self.samples[:,2], bins=40)
			x = (hist[1][0:-1]+hist[1][1:])/2
			hist = hist[0]
			self.alpha = x[np.argmax(hist)]
			self.unc_alp = np.nanstd(self.samples[:,2])
			# self.unc_alp = fwhm(x, hist)

			# gt_fwhm = hist >= (np.max(hist)/2)
			# self.unc_alp = (np.max(x[gt_fwhm]) - np.min(x[gt_fwhm]))/2
		else:
			raise ValueError("Method '%s' has not been programed in yet."%(
				method))

		self.prob = np.exp(self.lnprob((self.age, self.metallicity, self.alpha)))

		
		if self.pp is None and ab_index is None:
			vin_dir = get_vin_dir(self.instrument)
			self.vout_dir = '%s/%s/%s/pop' % (vin_dir, self.galaxy, self.opt)
			if not os.path.exists(self.vout_dir): os.makedirs(self.vout_dir)
			self.save()
			self.plot_probability_distribution(
				saveTo="%s/plots/%i.png" % (self.vout_dir, self.bin))
		else:
			print 'Age: ', self.age, '+/-', self.unc_age
			print 'Metallicity: ', self.metallicity, '+/-', self.unc_met
			print 'Alpha-enhancement: ', self.alpha, '+/-', self.unc_alp

#############################################################################

	def plot_probability_distribution(self, saveTo=None, f=None, ax_array=None,
		label='', legend =False):
		import matplotlib.pyplot as plt

		if f is None:
			f, ax_array = plt.subplots(2,2)
			alpha = 1
		else: alpha = 0.5

		if self.galaxy is not None:
			f.suptitle('%s Probability Distribution' % (self.galaxy.upper()))
		else:
			f.suptitle('Probability Distribution')
		ax_array[0,0].hist(self.samples[:,0],bins=40,histtype='step',normed=True, 
			alpha=alpha, label='Probability Distribution%s' % (label))
		ax_array[0,1].hist(self.samples[:,1],bins=40,histtype='step',normed=True, 
			alpha=alpha)
		ax_array[1,0].hist(self.samples[:,2],bins=40,histtype='step',normed=True, 
			alpha=alpha)

		ax_array[0,0].axvline(self.age - self.unc_age, color='r', alpha=alpha)
		ax_array[0,0].axvline(self.age + self.unc_age, color='r', alpha=alpha,
			label='Uncertainty of population fit%s' % (label))
		ax_array[0,0].axvline(self.age, alpha=alpha,
			label='Bestfitting population%s' % (label))
		ax_array[0,1].axvline(self.metallicity - self.unc_met, color='r', alpha=alpha)
		ax_array[0,1].axvline(self.metallicity + self.unc_met, color='r', alpha=alpha)
		ax_array[0,1].axvline(self.metallicity, alpha=alpha)
		ax_array[1,0].axvline(self.alpha - self.unc_alp, color='r', alpha=alpha)
		ax_array[1,0].axvline(self.alpha + self.unc_alp, color='r', alpha=alpha)
		ax_array[1,0].axvline(self.alpha, alpha=alpha)
		ax_array[0,0].set_title('Age')
		ax_array[0,1].set_title('Metallicity')
		ax_array[1,0].set_title('Alpha/Fe ratio')
		ax_array[1,1].axis('off')
		# plt.tight_layout()

		if legend:
			h, l = ax_array[0,0].get_legend_handles_labels()
			ax_array[1,1].legend(h,l)

		if saveTo is not None:# and self.pp is None:
			if not os.path.exists(os.path.dirname(saveTo)):
				os.makedirs(os.path.dirname(saveTo))
			f.savefig(saveTo)

		self.fig = f
		self.ax = ax_array

		# if not cc.remote:
		# 	plt.show()

#############################################################################

	def lnprob(self, theta):
		f_age,f_metal,f_alpha = theta
		if f_age<min(self._age) or f_age>max(self._age):
			return -np.inf
		if f_metal<min(self._metallicity) or f_metal>max(self._metallicity):
			return -np.inf
		if f_alpha<min(self._alpha) or f_alpha>max(self._alpha):
			return -np.inf

		chi2 = 0
		for l in self.lines:
			if ~np.isnan(self.ab_lines[l]) and ~np.isnan(self.uncerts[l]):
				chi2 += (self.ab_lines[l] - self.interp[l](theta))**2/\
					self.uncerts[l]**2
		return -chi2/2


	def save(self):
		if self.pp is not None:
			raise ValueError('population was run using the ppxf output object and ' +
				'therefore has no obvious place to save the output')
		file = "%s/%i.dat" % (self.vout_dir, self.bin)
		with open(file, 'w') as f:
			f.write('%f   %f   %f \n' % (self.age, self.metallicity, self.alpha))
			f.write('%f   %f   %f ' % (self.unc_age, self.unc_met, self.unc_alp))


		file = "%s/distribution/%i.dat" % (self.vout_dir, self.bin)
		if not os.path.exists(os.path.dirname(file)): 
			os.makedirs(os.path.dirname(file))
		s = self.samples.shape

		with open(file, 'w') as f:
			for j in range(s[0]):
				f.write('%f   %f   %f \n' % (self.samples[j,0], self.samples[j,1],
					self.samples[j,2]))

#############################################################################

if __name__=="__main__":
	p = population()