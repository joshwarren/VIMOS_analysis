## ==================================================================
## 		Stellar population
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
#from sauron_colormap2 import sauron2 as sauron
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from checkcomp import checkcomp
cc = checkcomp()

def stellar_pop(galaxy, wav_range="", vLimit=0):
	grid_length = 40
	# Find lines:
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		'Mg_1', 'Mg_2', 'Mg_b']

	print 'Stellar populations'

	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""

	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_plots = "%splots" % (output)
	out_pickle = '%s/pickled' % (output)
	pickleFile = open("%s/dataObj_%s.pkl" % (out_pickle, wav_range), 'rb')
	#pickleFile = open("%s/dataObj_%s.pkl" % (cc.home_dir, wav_range), 'rb')

	D = pickle.load(pickleFile)
	pickleFile.close()

	models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)
	titles = np.loadtxt(models_dir, dtype=str)[0]

	age, metalicity, alpha = np.loadtxt(models_dir, usecols=(0,1,2), 
		unpack=True, skiprows=36)

	age2 = np.unique(age)
	metalicity2 = np.unique(metalicity)
	alpha2 = np.unique(alpha)

	age = np.logspace(np.log10(min(age2)), np.log10(max(age2)), 
		num=grid_length)
	age[np.argmax(age)] = max(age2)
	metalicity = np.linspace(min(metalicity2), max(metalicity2), 
		num=grid_length)
	alpha = np.linspace(min(alpha2), max(alpha2), num=grid_length)

	models = {}
	interp = {}
	n_lines = np.zeros(D.number_of_bins)
	print '    Interpolating models'
	for i, l in enumerate(titles):
		if l in lines:
			models[l] = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
				skiprows=36).reshape((len(age2),len(metalicity2),len(alpha2)))
			interp[l] = RegularGridInterpolator((age2, metalicity2, alpha2), 
				models[l])
			n_lines += (~np.isnan(D.absorption_line(l))).astype(int)

	
	
	chi2 = np.zeros((grid_length,grid_length,grid_length, D.number_of_bins))
	for line in lines:
		print '    Fitting ' + line
		ab_line, uncert = D.absorption_line(line, uncert=True)
		d = np.array([~np.isnan(ab_line), ~np.isnan(uncert)]).all(axis=0)
		for i, ag in enumerate(age):
			for j, me in enumerate(metalicity):
				for k, al in enumerate(alpha):
					chi2[i,j,k,d] += np.square(ab_line[d] -	
						interp[line]([ag,me,al]))/(uncert[d]*n_lines[d])

	# f = 'test.pkl'
	# p = open(f, 'wb')
	# pickle.dump(chi2,p)
	# p.close()

	chi2[chi2==0] = np.nan
	#chi2[chi2 > 100] = np.nan

	# Finding locations of minimum chi2 for each bin
	#a = [np.unravel_index(np.nanargmin(chi2[:,:,:,i]),chi2[:,:,:,i].shape)
	#	for i in range(D.number_of_bins)]
	a = []
	#nans = []
	for i in range(D.number_of_bins):
		try:
			a.append(np.unravel_index(np.nanargmin(chi2[:,:,:,i]),chi2[:,:,:,i].shape))
	#		nans.append(False)
		except ValueError:
			print i
			a.append(np.unravel_index(np.argmin(chi2[:,:,:,i]),chi2[:,:,:,i].shape))
	#		nans.append(True)
	#nans = np.ravel(nans)


	chi2 = np.array([chi2[i[0],i[1],i[2],j] for j, i in enumerate(a)])
	age = np.array([age[i[0]] for i in a])
	metalicity = np.array([metalicity[i[1]] for i in a])
	alpha = np.array([alpha[i[2]] for i in a])

	nans = chi2>100
	chi2[nans] = np.nan
	age[nans] = np.nan
	metalicity[nans] = np.nan
	alpha[nans] = np.nan

	# Produce and Save Plots
	d = {'chi2':chi2, 'age':age,'metalicity':metalicity,'alpha':alpha}
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
	f.savefig(saveTo, bbox_inches="tight")





#############################################################################

# Use of stellar_pop.py

if __name__ == '__main__':
	stellar_pop('ic1459', wav_range='4200-', vLimit=2)
