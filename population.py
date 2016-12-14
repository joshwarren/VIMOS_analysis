## ==================================================================
## 		Stellar population
## ==================================================================
import numpy as np
from scipy.interpolate import griddata
from checkcomp import checkcomp
cc = checkcomp()

def population(ab_lines, uncerts, interp=None):
	grid_length = 40
	# Find lines:
	lines = ab_lines.keys()

	models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)
	titles = np.loadtxt(models_dir, dtype=str)[0]

	age_in, metallicity_in, alpha_in = np.loadtxt(models_dir, usecols=(0,1,2), 
		unpack=True, skiprows=35)

	age = np.logspace(np.log10(min(age_in)), np.log10(max(age_in)), 
		num=grid_length)
	age[np.argmax(age)] = max(age_in)
	metallicity = np.linspace(min(metallicity_in), max(metallicity_in), 
		num=grid_length)
	alpha = np.linspace(min(alpha_in), max(alpha_in), num=grid_length)

	age_p=np.array([[m]*grid_length**2 for m in age]).flatten()
	metallicity_p=np.array([[m]*grid_length for m in metallicity]*grid_length
		).flatten()
	alpha_p=np.array(list(alpha)*grid_length**2)

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
					).reshape((grid_length,grid_length,grid_length))

	s = interp[lines[0]].shape
	chi2 = np.zeros((s[0], s[1], s[2]))
	for line in lines:
		d = np.array([~np.isnan(ab_lines[line]), ~np.isnan(uncerts[line])]).all(axis=0)
		for i, ag in enumerate(age):
			for j, me in enumerate(metallicity):
				for k, al in enumerate(alpha):
					chi2[i,j,k] += (ab_lines[line] - interp[line][i,j,k])**2/\
						uncerts[line][d]**2

	chi2_m = chi2/np.min(chi2[:,:,:])

	i = np.argmax(np.nansum(np.exp(-np.square(chi2_m[:,:,:])/2),axis=(1,2)))
	j = np.argmax(np.nansum(np.exp(-np.square(chi2_m[:,:,:])/2),axis=(0,2)))
	k = np.argmax(np.nansum(np.exp(-np.square(chi2_m[:,:,:])/2),axis=(0,1)))

	return age[i], metallicity[j], alpha[k], chi2[i,j,k]
#############################################################################