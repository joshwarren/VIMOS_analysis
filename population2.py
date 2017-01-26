## ==================================================================
## 		Stellar population
## ==================================================================
## A method to impliment MCMC routines to find the bestfit model to the 
## indices strengths given.
## warrenj 20170126 Added previous keyword to load from pop directory that 
## is created in the save routine. This was an attempt to create something a 
## little quicker than pickling (not sure it has worked...)
##
#######################################################################
# Keywords:
# ab_lines: 			Dict of absorption line strength. NB: entry can be 
#						array to fit more than one bin simulataniously.
# uncerts:				Dict of uncertainty in the absorption line strenghts.
# interp:		None 	Dict of interpolated models. If None, models are read 
#						and interpolated everytime routine is called.
# grid_length	40		Int of the length of a side of the cube of interpolated 
#						models. NB: This is not used if interp is given.
#######################################################################
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import emcee
from tools import length as len
from checkcomp import checkcomp
cc = checkcomp()

class population(object):

	def __init__(self, ab_lines, uncerts, interp=None, grid_length=40, previous=False, 
		galaxy=None):
		# Absorption lines provided:
		self.lines = ab_lines.keys()
		self.ab_lines = ab_lines
		self.uncerts = uncerts
		self.nbins = len(ab_lines[self.lines[0]])
		if self.nbins==1 and (type(ab_lines[self.lines[0]])!=list or \
			type(ab_lines[self.lines[0]])!=np.ndarray):
			for key in ab_lines:
				ab_lines[key] = np.array([ab_lines[key]])
		if interp is None:
			s=[grid_length,grid_length,grid_length]
		else:
			s = interp[self.lines[0]].shape

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

		# if interp is None:
		models = {}
		self.interp = {}
		# 	print '    Interpolating models'
		for i, l in enumerate(titles):
			if l in self.lines:
				models[l] = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
					skiprows=35)
				self.interp[l] = LinearNDInterpolator(
					np.array([age_in, metallicity_in, alpha_in]).transpose(), 
					models[l])

		ndim, nwalkers, nsteps = 3, 200, 500
		self.samples = np.zeros((self.nbins, nwalkers*(nsteps-50), ndim))
		if previous:
			print '    Loading'
			for bin in range(self.nbins):
				a,m,al = np.loadtxt('%s/Data/vimos/analysis/%s/' % (cc.base_dir, galaxy) +
					'results/4200-/pickled/pop/%i.dat' % (bin), unpack=True)
				self.samples[bin,:,0] = a 
				self.samples[bin,:,1] = m 
				self.samples[bin,:,2] = al
		else:
			print '    Fitting'
			for bin in range(self.nbins):

				sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, args=[bin])
				pos = [np.array([13,0.2,-0.1]) +1e-3*np.random.randn(ndim) for i in 
					range(nwalkers)]
				sampler.run_mcmc(pos, 500)

				# Discard initial steps away from initial point
				self.samples[bin,:,:] = sampler.chain[:, 50:, :].reshape((-1, ndim))

				if galaxy is not None:
					self.save(galaxy)

		self.age = np.nanmean(self.samples[:,:,0], axis=1)
		self.unc_age = np.nanstd(self.samples[:,:,0], axis=1)
		self.metallicity = np.nanmean(self.samples[:,:,1], axis=1)
		self.unc_met = np.nanstd(self.samples[:,:,1], axis=1)
		self.alpha = np.nanmean(self.samples[:,:,2], axis=1)
		self.unc_alp = np.nanstd(self.samples[:,:,2], axis=1)

#############################################################################

	def plot_probability_distribution(self, galaxy=None):
		import matplotlib.pyplot as plt

		f, ax_array = plt.subplots(2,2)
		if galaxy is not None:
			f.suptitle('%s Probability Distribution' % (galaxy.upper()))
		else:
			f.suptitle('Probability Distribution')
		for b in range(self.nbins):
			ax_array[0,0].hist(self.samples[b,:,0],bins=40,histtype='step',normed=True)
			ax_array[0,1].hist(self.samples[b,:,1],bins=40,histtype='step',normed=True)
			ax_array[1,0].hist(self.samples[b,:,2],bins=40,histtype='step',normed=True)

			ax_array[0,0].axvline(self.age[b])
		ax_array[0,0].set_title('Age')
		ax_array[0,1].set_title('Metallicity')
		ax_array[1,0].set_title('Alpha/Fe ratio')
		ax_array[1,1].axis('off')
		plt.show()


#############################################################################

	def lnprob(self, theta, bin):
		f_age,f_metal,f_alpha = theta
		if f_age<min(self._age) or f_age>max(self._age):
			return -np.inf
		if f_metal<min(self._metallicity) or f_metal>max(self._metallicity):
			return -np.inf
		if f_alpha<min(self._alpha) or f_alpha>max(self._alpha):
			return -np.inf

		chi2 = 0
		for l in self.lines:
			if ~np.isnan(self.ab_lines[l][bin]) and ~np.isnan(self.uncerts[l][bin]):
				chi2 += (self.ab_lines[l][bin] - self.interp[l](theta))**2/\
					self.uncerts[l][bin]**2
		return -chi2/2


	def save(self, galaxy):
		import os
		if not os.path.exists("%s/Data/vimos/analysis/%s/results/4200-/pickled/pop/" % (
			cc.base_dir, galaxy)):
			os.makedirs("%s/Data/vimos/analysis/%s/results/4200-/pickled/pop/" % (
				cc.base_dir, galaxy))
		for i in range(self.nbins):
			file = "%s/Data/vimos/analysis/%s/results/4200-/pickled/pop/%i.dat" % (
				cc.base_dir, galaxy, i)

			f = open(file, 'w')
			s = self.samples.shape
			for j in range(s[1]):
				f.write('%f   %f   %f \n' % (self.samples[i,j,0], self.samples[i,j,1],
					self.samples[i,j,2]))

#############################################################################

if __name__=="__main__":
	import cPickle as pickle

	gal='ngc3100'
	bin=25
	wav_range='4200-'

	out_pickle = '%s/Data/vimos/analysis/%s/results/%s/pickled' % (cc.base_dir,gal,wav_range)
	pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		'Mg_b']

	line_dir = {}
	uncert_dir = {}
	for l in lines:
		ab, uncert = D.absorption_line(l, uncert=True)
		line_dir[l] = np.array([ab[bin]])
		uncert_dir[l] = np.array([uncert[bin]])

	pop = population(line_dir, uncert_dir, grid_length=40)