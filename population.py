## ==================================================================
## 		Stellar population
## ==================================================================
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
from scipy.interpolate import griddata
from tools import length as len
from checkcomp import checkcomp
cc = checkcomp()

class population(object):

	def __init__(self, ab_lines, uncerts, interp=None, grid_length=40):
		# Absorption lines provided:
		lines = ab_lines.keys()
		n = len(ab_lines[lines[0]])
		if n==1 and (type(ab_lines[lines[0]])!=list or type(ab_lines[lines[0]])!=np.ndarray):
			for key in ab_lines:
				ab_lines[key] = np.array([ab_lines[key]])
		if interp is None:
			s=[grid_length,grid_length,grid_length]
		else:
			s = interp[lines[0]].shape

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

		if interp is None:
			models = {}
			interp = {}
			print '    Interpolating models'
			for i, l in enumerate(titles):
				if l in lines:
					models[l] = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
						skiprows=35)
					interp[l] = griddata(
						np.array([age_in, metallicity_in, alpha_in]).transpose(), 
						models[l], np.array([age_p, metallicity_p, alpha_p]).transpose()
						).reshape((s[0],s[1],s[2]))
			s = interp[lines[0]].shape

		print '    Fitting'
		chi2 = np.zeros((s[0], s[1], s[2], n))
		n_lines =  np.zeros(n)
		for line in lines:
			ab = ab_lines[line]
			uncert = uncerts[line]

			d = np.array([~np.isnan(ab), ~np.isnan(uncert)]).all(axis=0)
			n_lines += d
			for i in range(s[0]): 			# age
				for j in range(s[1]): 		# metallicity
					for k in range(s[2]):	# alpha
						chi2[i,j,k,d] += (ab[d] - interp[line][i,j,k])**2/uncert[d]**2
		chi2[chi2==0] = np.nan

		chi2_m = np.array(chi2)
		for b in range(n):
			chi2_m[:,:,:,b] -= np.nanmin(chi2_m[:,:,:,b])

		self.age_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(1,2))
		self.metal_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(0,2))
		self.alpha_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(0,1))
		i = np.argmax(self.age_prob_dist, axis=0)
		j = np.argmax(self.metal_prob_dist,axis=0)
		k = np.argmax(self.alpha_prob_dist,axis=0)

		self.age = self._age[i]
		self.metallicity = self._metallicity[j]
		self.alpha = self._alpha[k]
		self.red_chi2 = chi2[i,j,k,range(n)]/(n_lines-3)

		self.unc_age = np.sqrt(np.sum(self.age_prob_dist*(self._age[:,np.newaxis] - 
			self.age)**2,axis=0)/np.sum(self.age_prob_dist,axis=0))
		self.unc_met = np.sqrt(np.sum(self.metal_prob_dist*(self._metallicity[:,np.newaxis] - 
			self.metallicity)**2,axis=0)/np.sum(self.metal_prob_dist,axis=0))
		self.unc_alp = np.sqrt(np.sum(self.alpha_prob_dist*(self._alpha[:,np.newaxis] - 
			self.alpha)**2,axis=0)/np.sum(self.alpha_prob_dist,axis=0))

#############################################################################

	def plot_probability_distribution(self, galaxy=None):
		import matplotlib.pyplot as plt

		f, ax_array = plt.subplots(2,2)
		if galaxy is not None:
			f.suptitle('%s Probability Distribution' % (galaxy.upper()))
		else:
			f.suptitle('Probability Distribution')
		for b in range(len(self.age)):
			ax_array[0,0].plot(self._age, self.age_prob_dist[:,b]/np.nansum(
				self.age_prob_dist[:,b]))
			ax_array[0,1].plot(self._metallicity, self.metal_prob_dist[:,b]/np.nansum(
				self.metal_prob_dist[:,b]))
			ax_array[1,0].plot(self._alpha, self.alpha_prob_dist[:,b]/np.nansum(
				self.alpha_prob_dist[:,b]))
		ax_array[0,0].set_title('Age')
		ax_array[0,1].set_title('Metallicity')
		ax_array[1,0].set_title('Alpha/Fe ratio')
		ax_array[1,1].axis('off')
