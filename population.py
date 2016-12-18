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

def population(ab_lines, uncerts, interp=None, grid_length=40):
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

	age = np.logspace(np.log10(min(age_in)), np.log10(max(age_in)), num=s[0])
	age[np.argmax(age)] = max(age_in)
	metallicity = np.linspace(min(metallicity_in), max(metallicity_in), num=s[1])
	alpha = np.linspace(min(alpha_in), max(alpha_in), num=s[2])

	age_p=np.array([[m]*s[2]*s[1] for m in age]).flatten()
	metallicity_p=np.array([[m]*s[2] for m in metallicity]*s[0]).flatten()
	alpha_p=np.array(list(alpha)*s[1]*s[0])

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

	age_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(1,2))
	metal_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(0,2))
	alpha_prob_dist = np.nansum(np.exp(-np.square(chi2_m)/2),axis=(0,1))
	i = np.argmax(age_prob_dist, axis=0)
	j = np.argmax(metal_prob_dist,axis=0)
	k = np.argmax(alpha_prob_dist,axis=0)

	return age[i], metallicity[j], alpha[k], chi2[i,j,k,range(n)]/(n_lines-3)
#############################################################################