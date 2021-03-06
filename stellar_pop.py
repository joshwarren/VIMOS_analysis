## ==================================================================
## 		Stellar population
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
#from sauron_colormap2 import sauron2 as sauron
import numpy as np
import os
from population2 import population
from checkcomp import checkcomp
cc = checkcomp()

def stellar_pop(galaxy, wav_range="", vLimit=0, D=None, previous=False):
	print 'Stellar populations'

	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range)
	out_plots = "%s/plots" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)

	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	if previous:
		out_pickle = '%s/pickled' % (output)
		print "    Load Stellar population object"
		if not os.path.exists(out_pickle):
			os.makedirs(out_pickle) 
		pickleFile = open("%s/stellarPopObj_%s.pkl" % (out_pickle, wav_range), 'rb')
		pop = pickle.load(pickleFile)
		pickleFile.close()

	else:
		grid_length = 40
		# Find lines:
		lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
			'Mg_b']
		#lines = ['H_beta']
		
		line_dir = {}
		uncert_dir = {}
		for l in lines:
			ab, uncert = D.absorption_line(l, uncert=True)
			line_dir[l] = ab
			uncert_dir[l] = uncert

		pop = population(line_dir, uncert_dir, grid_length=grid_length)


	# Produce and Save Plots
	d = {'age':pop.age,'metallicity':pop.metallicity,'alpha':pop.alpha}
	d_uncert = {'age':pop.unc_age,'metallicity':pop.unc_met,
		'alpha':pop.unc_alp}
	c_label = {'chi2':'', 'age':'Gyrs','metallicity':'[Z/H]','alpha':'[alpha/Fe]'}
	f, ax_array = plt.subplots(2, 2, sharex='col', sharey='row')
	f2, ax_array2 = plt.subplots(2, 2, sharex='col', sharey='row')
	i=0
	print '    Plotting and saving'
	for plot, values in d.iteritems():
		
		if plot=='chi2':
			vmin = sorted(values[~np.isnan(values)])[vLimit]
			vmax = sorted(values[~np.isnan(values)])[-1-vLimit]
			
			mean = np.nanmedian(values)
			std = np.nanstd(values)

			vmin = max([vmin, mean-2*std])
			vmax = min([vmax, mean+2*std])
		elif plot=='age':
			vmin,vmax=0,15
		elif plot=='metallicity':
			vmin,vmax=-2.25,0.67
		elif plot=='alpha':
			vmin,vmax=-0.3,0.5

		ax_array[i%2,int(np.floor(i/2))] = plot_velfield_nointerp(D.x, D.y, 
			D.bin_num, D.xBar, D.yBar, values, vmin=vmin, vmax=vmax,
			nodots=True, colorbar=True, title=plot, label=c_label[plot],
			ax=ax_array[i%2,int(np.floor(i/2))], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)

		ax_array2[i%2,int(np.floor(i/2))] = plot_velfield_nointerp(D.x, D.y, 
			D.bin_num, D.xBar, D.yBar, d_uncert[plot], #vmin=vmin, vmax=vmax,
			nodots=True, colorbar=True, title='Uncertainty in '+plot, label=c_label[plot],
			ax=ax_array2[i%2,int(np.floor(i/2))], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)
		i+=1
	saveTo = "%s/stellar_pop_%s.pdf" % (out_plots, wav_range)
	f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	f.suptitle(galaxy.upper())
	if not os.path.exists(os.path.dirname(saveTo)):
		os.makedirs(os.path.dirname(saveTo))  
	f.savefig(saveTo, bbox_inches="tight")

	saveTo = "%s/stellar_pop_%s_uncert.pdf" % (out_plots, wav_range)
	f2.tight_layout()
	ax_array2[0,1].set_xlabel('')
	ax_array2[0,0].set_xlabel('')
	ax_array2[0,1].set_ylabel('')
	ax_array2[1,1].set_ylabel('')
	f2.suptitle(galaxy.upper())
	f2.savefig(saveTo, bbox_inches="tight")


# ------------================ Pickling =================----------
	if not previous:
		print "    Pickling Stellar population object"
		if not os.path.exists(out_pickle):
			os.makedirs(out_pickle) 
		pickleFile = open("%s/stellarPopObj_%s.pkl" % (out_pickle, wav_range), 'wb')
		pickle.dump(pop,pickleFile)
		pickleFile.close()

	return D



#############################################################################

# Use of stellar_pop.py

if __name__ == '__main__':
	stellar_pop('ic1531', wav_range='4200-', vLimit=2, previous=True)

