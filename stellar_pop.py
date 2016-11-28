## ==================================================================
## 		Stellar population
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
#from sauron_colormap2 import sauron2 as sauron
import numpy as np
import os
from scipy.interpolate import griddata
from checkcomp import checkcomp
cc = checkcomp()

def stellar_pop(galaxy, wav_range="", vLimit=0, D=None):
	grid_length = 40
	# Find lines:
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		#'Mg_1', 'Mg_2', 
		'Mg_b']

	print 'Stellar populations'

	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""

	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_plots = "%splots" % (output)

	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
		#pickleFile = open("%s/dataObj_%s.pkl" % (cc.home_dir, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)
	titles = np.loadtxt(models_dir, dtype=str)[0]

	age3, metallicity3, alpha3 = np.loadtxt(models_dir, usecols=(0,1,2), 
		unpack=True, skiprows=36)

	age2 = np.unique(age3)
	metallicity2 = np.unique(metallicity3)
	alpha2 = np.unique(alpha3)

	age = np.logspace(np.log10(min(age2)), np.log10(max(age2)), 
		num=grid_length)
	age[np.argmax(age)] = max(age2)
	metallicity = np.linspace(min(metallicity2), max(metallicity2), 
		num=grid_length)
	alpha = np.linspace(min(alpha2), max(alpha2), num=grid_length)

	age_p=np.array([[m]*grid_length**2 for m in age]).flatten()
	metallicity_p=np.array([[m]*grid_length for m in metallicity]*grid_length
		).flatten()
	alpha_p=np.array(list(alpha)*grid_length**2)

	models = {}
	interp = {}
	n_lines = np.zeros(D.number_of_bins)
	print '    Interpolating models'
	for i, l in enumerate(titles):
		if l in lines:
			models[l] = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
				skiprows=36)
			interp[l] = griddata(np.array([age3, metallicity3, alpha3]).transpose(), 
				models[l], np.array([age_p, metallicity_p, alpha_p]).transpose()
				).reshape((grid_length,grid_length,grid_length))
			n_lines += (~np.isnan(D.absorption_line(l))).astype(int)

	
	
	chi2 = np.zeros((grid_length,grid_length,grid_length, D.number_of_bins))
	for line in lines:
		print '    Fitting ' + line
		ab_line, uncert = D.absorption_line(line, uncert=True)
		d = np.array([~np.isnan(ab_line), ~np.isnan(uncert)]).all(axis=0)
		for i, ag in enumerate(age):
			for j, me in enumerate(metallicity):
				for k, al in enumerate(alpha):
					chi2[i,j,k,d] += np.square(ab_line[d] -	
						interp[line][ag,me,al])/((uncert[d]**2)*n_lines[d])

	# f = 'test.pkl'
	# p = open(f, 'wb')
	# pickle.dump(chi2,p)
	# p.close()

	chi2[chi2==0] = np.nan
	age_map=np.zeros(D.number_of_bins)
	metal_map=np.zeros(D.number_of_bins)
	alpha_map=np.zeros(D.number_of_bins)
	chi2_map=np.zeros(D.number_of_bins)
	for bin in range(D.number_of_bins):
		i = np.argmax(np.nansum(np.exp(-chi2[:,:,:,bin]**2/2),axis=(1,2)))
		j = np.argmax(np.nansum(np.exp(-chi2[:,:,:,bin]**2/2),axis=(0,2)))
		k = np.argmax(np.nansum(np.exp(-chi2[:,:,:,bin]**2/2),axis=(0,1)))
		age_map[bin] = age[i]
		metal_map[bin] = metallicity[j]
		alpha_map[bin] = alpha[k]
		chi2_map[bin] = chi2[i,j,k,bin]





	# Produce and Save Plots
	d = {'chi2':chi2_map, 'age':age_map,'metallicity':metal_map,'alpha':alpha_map}
	f, ax_array = plt.subplots(2, 2, sharex='col', sharey='row')
	i=0
	print '    Plotting and saving'
	for plot, values in d.iteritems():
		vmin = sorted(values[~np.isnan(values)])[vLimit]
		vmax = sorted(values[~np.isnan(values)])[-1-vLimit]

		ax_array[i%2,np.floor(i/2)] = plot_velfield_nointerp(D.x, D.y, 
			D.bin_num, D.xBar, D.yBar, values, vmin=vmin, vmax=vmax,
			nodots=True, colorbar=True, title=plot, 
			ax=ax_array[i%2,np.floor(i/2)], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)
		i+=1
	saveTo = "%s/stellar_pop_%s.pdf" % (out_plots, wav_range)
	#saveTo = '%s/test.pdf' % (cc.home_dir)
	f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	f.suptitle(galaxy.upper())
	if not os.path.exists(os.path.dirname(saveTo)):
		os.makedirs(os.path.dirname(saveTo))  
	f.savefig(saveTo, bbox_inches="tight")

	return D



#############################################################################

# Use of stellar_pop.py

if __name__ == '__main__':
	stellar_pop('ic1459', wav_range='4200-', vLimit=2)

