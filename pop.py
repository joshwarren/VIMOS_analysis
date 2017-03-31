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

c = 299792.458

if cc.device == 'glamdring': vin_dir = '%s/analysis' % (cc.base_dir)
else: vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)

class population(object):

	def __init__(self):#, i_gal=None, bin=None):
		#ab_lines, uncerts, interp=None, grid_length=40, previous=False, self.galaxy=None):
		self.i_gal=int(sys.argv[1])
		self.bin=int(sys.argv[2])

		self.lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']


		galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
		self.galaxy = galaxies[self.i_gal]

		self.vout_dir = '%s/%s/pop' % (vin_dir, self.galaxy)
		if not os.path.exists(self.vout_dir): os.makedirs(self.vout_dir)

		grid_length = 40

		vin_dir_gasMC = "%s/%s/pop_MC" % (vin_dir, self.galaxy)
		
		lam = np.loadtxt("%s/lambda/%d.dat" % (vin_dir_gasMC, self.bin))
		spectrum = np.loadtxt("%s/input/%d.dat" % (vin_dir_gasMC, self.bin))
		matrix = np.loadtxt("%s/bestfit/matrix/%d.dat" % (vin_dir_gasMC, self.bin), 
			dtype=str)
		e_lines = np.array([not s.isdigit() for s in matrix[:,0]])
		e_line_spec = matrix[e_lines,1:].astype(float)
		temp_weights =  np.loadtxt("%s/temp_weights/%d.dat" % (vin_dir_gasMC, self.bin), 
			unpack=True, usecols=(1,))
		e_line_spec = np.einsum('ij,i->ij',e_line_spec,temp_weights[e_lines])
		continuum = spectrum - np.nansum(e_line_spec,axis=0)
		noise = np.loadtxt("%s/noise_input/%d.dat" % (vin_dir_gasMC, self.bin))
		bestfit = np.loadtxt("%s/bestfit/%d.dat" %(vin_dir_gasMC, self.bin))
		convolved = bestfit - np.nansum(e_line_spec, axis=0) 

		if cc.getDevice() == 'uni':
			files = glob('%s/Data/idl_libraries/ppxf/MILES_library/' % (cc.base_dir) +
				'm0[0-9][0-9][0-9]V')
		else:
			files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
		wav = np.loadtxt(files[0], usecols=(0,), unpack=True)

		a = [min(np.where(wav>=lam[0])[0]), max(np.where(wav<=lam[-1])[0])]
		unconvolved_spectrum = np.zeros(a[1]-a[0])

		for i, n in enumerate(matrix[~e_lines,0]):
			template =  np.loadtxt(files[int(n)], usecols=(1,), unpack=True)
			unconvolved_spectrum += template[a[0]:a[1]]*temp_weights[i]
		mpweight = np.loadtxt("%s/mpweights/%d.dat" %(vin_dir_gasMC, self.bin))
		unconvolved_spectrum *= np.polynomial.legendre.legval(np.linspace(-1,1,
				len(unconvolved_spectrum)), np.append(1, mpweight))
		unconvolved_lam =  wav[a[0]:a[1]]

		self.ab_lines = {}
		self.uncerts = {}
		for l in self.lines:
			ab, uncert = absorption(l, lam, continuum, noise=noise,
				unc_lam=unconvolved_lam, unc_spec=unconvolved_spectrum, 
				conv_spec=convolved)
			self.ab_lines[l], self.uncerts[l] = ab[0], uncert[0]


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
		
		self.save()
		self.plot_probability_distribution()


#############################################################################

	def plot_probability_distribution(self):
		import matplotlib.pyplot as plt

		f, ax_array = plt.subplots(2,2)
		if self.galaxy is not None:
			f.suptitle('%s Probability Distribution' % (self.galaxy.upper()))
		else:
			f.suptitle('Probability Distribution')
		ax_array[0,0].hist(self.samples[:,0],bins=40,histtype='step',normed=True)
		ax_array[0,1].hist(self.samples[:,1],bins=40,histtype='step',normed=True)
		ax_array[1,0].hist(self.samples[:,2],bins=40,histtype='step',normed=True)

		ax_array[0,0].axvline(self.age - self.unc_age, color='r')
		ax_array[0,0].axvline(self.age + self.unc_age, color='r')
		ax_array[0,0].axvline(self.age)
		ax_array[0,1].axvline(self.metallicity - self.unc_met, color='r')
		ax_array[0,1].axvline(self.metallicity + self.unc_met, color='r')
		ax_array[0,1].axvline(self.metallicity)
		ax_array[1,0].axvline(self.alpha - self.unc_alp, color='r')
		ax_array[1,0].axvline(self.alpha + self.unc_alp, color='r')
		ax_array[1,0].axvline(self.alpha)
		ax_array[0,0].set_title('Age')
		ax_array[0,1].set_title('Metallicity')
		ax_array[1,0].set_title('Alpha/Fe ratio')
		ax_array[1,1].axis('off')
		plt.tight_layout()

		if not os.path.exists("%s/plots/" % (self.vout_dir)):
			os.makedirs("%s/plots/" % (self.vout_dir))
		file = "%s/plots/%i.png" % (self.vout_dir, self.bin)
		f.savefig(file)

		if not cc.remote:
			plt.show()

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



	# import cPickle as pickle

	# gal='ngc3100'
	# self.bin=25
	# wav_range='4200-'

	# out_pickle = '%s/Data/vimos/analysis/%s/results/%s/pickled' % (cc.base_dir,gal,
	# wav_range)
	# pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	# D = pickle.load(pickleFile)
	# pickleFile.close()

	# lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
	# 	'Mg_b']

	# line_dir = {}
	# uncert_dir = {}
	# for l in lines:
	# 	ab, uncert = D.absorption_line(l, uncert=True)
	# 	line_dir[l] = np.array([ab[self.bin]])
	# 	uncert_dir[l] = np.array([uncert[self.bin]])

	# pop = population(line_dir, uncert_dir, grid_length=40)