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
from prefig import Prefig 
Prefig()

c = 299792.458

def get_vin_dir(instrument):
	if instrument=='vimos':
		if cc.device == 'glamdring': vin_dir = '%s/analysis' % (cc.base_dir)
		else: vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	elif instrument=='muse':
		if cc.device == 'glamdring': vin_dir = '%s/analysis_muse' % (cc.base_dir)
		else: vin_dir = '%s/Data/muse/analysis' % (cc.base_dir)
	return vin_dir

def get_absorption(lines, pp=None, galaxy=None, bin=None, opt=None, 
	instrument='vimos', sigma=None):
	if (pp is None and galaxy is None) or (pp is not None and galaxy is not None):
		raise ValueError('Either pp or galaxy must be supplied: one or tother')
	if galaxy is not None and (bin is None or opt is None):
		raise ValueError('If galaxy is supplied, so to must bin and opt')

	vin_dir = get_vin_dir(instrument)
	if pp is None:
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

		library = pp.params.library

	e_line_spec = np.einsum('ij,i->j',e_line_spec,temp_weights[e_lines])

	continuum = spectrum - e_line_spec #np.nansum(e_line_spec,axis=0)

	if library == 'Miles':
		if cc.getDevice() == 'uni':
			files = glob('%s/Data/idl_libraries/ppxf/MILES_library/' % (cc.base_dir) +
				'm0[0-9][0-9][0-9]V')
		else:
			files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
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
		FWHM_tem = 1.8

		f = fits.open(templateFiles[0])
		unconvolved_lam = np.arange(f[0].header['NAXIS1'])*f[0].header['CDELT1'] + \
			f[0].header['CRVAL1']
		unconvolved_spectrum = np.zeros((f[0].header['NAXIS1'], len(templateFiles)))

		for i in range(len(templateFiles)):
			f = fits.open(templateFiles[i])
			unconvolved_spectrum[:,i] = f[0].data
		unconvolved_spectrum = unconvolved_spectrum.dot(temp_weights[~e_lines])
		unconvolved_spectrum *= np.polynomial.legendre.legval(np.linspace(-1,1,
			len(unconvolved_spectrum)), np.append(1, mpweight))
		CDELT = f[0].header['CDELT1']
		f.close()

	if sigma is not None:
		from scipy.ndimage.filters import gaussian_filter1d
		sig_pix = np.median(unconvolved_lam)*(sigma/c)/CDELT
		convolved = gaussian_filter1d(unconvolved_spectrum, sig_pix)
	else:
		convolved = bestfit - e_line_spec # np.nansum(e_line_spec, axis=0)
		

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
		ab_uncert=None, instrument='vimos'):
		self.pp = pp
		self.galaxy = galaxy
		self.instrument = instrument
		self.opt = opt


		self.lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
		grid_length = 40

		if ab_index is not None:
			if ab_uncert is None:
				raise ValueError('Uncertainty values must be supplied')
			self.ab_lines = ab_index
			self.uncerts = ab_uncert
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
					galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
						'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
				elif self.instrument == 'muse':
					galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
				self.galaxy = galaxies[self.i_gal]

				self.ab_lines, self.uncerts = get_absorption(self.lines, 
					galaxy=self.galaxy, bin=self.bin, opt=self.opt, 
					instrument=self.instrument)
			else:
				self.ab_lines, self.uncerts = get_absorption(self.lines, pp=self.pp,
					instrument=self.instrument)

		vin_dir = get_vin_dir(self.instrument)
		self.vout_dir = '%s/%s/%s/pop' % (vin_dir, self.galaxy, self.opt)
		if not os.path.exists(self.vout_dir): os.makedirs(self.vout_dir)

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

		self.age = np.nanmean(self.samples[:,0])
		self.unc_age = np.nanstd(self.samples[:,0])
		self.metallicity = np.nanmean(self.samples[:,1])
		self.unc_met = np.nanstd(self.samples[:,1])
		self.alpha = np.nanmean(self.samples[:,2])
		self.unc_alp = np.nanstd(self.samples[:,2])
		
		if self.pp is None and ab_index is None:
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
			h, l = ax_array[0,0].get_legend_handle_labels()
			ax_array[1,1].legend(h,l)

		if saveTo is not None and self.pp is None:
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