## ==================================================================
## Plot Absorption Line Strength for most-likely Stellar population
## ==================================================================
## Routine to return the absorption line strenghts for the most likely stellar 
## Populations as found by populations.py

import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
from scipy.interpolate import griddata
#from sauron_colormap2 import sauron2 as sauron
import numpy as np
import os
from population2 import population
from plot_results import set_lims
from tools import nearest
from checkcomp import checkcomp
cc = checkcomp()

def stellar_pop_self_consis(galaxy, wav_range="", vLimit=0, interp=None, D=None):
	print 'Stellar Populations Self Consistancy Check'

	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range)
	out_plots = "%s/plots" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)
	titles = np.loadtxt(models_dir, dtype=str)[0]
	age_in, metallicity_in, alpha_in = np.loadtxt(models_dir, usecols=(0,1,2), 
		unpack=True, skiprows=35)

	

	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	pickleFile = open("%s/stellarPopObj_%s.pkl" % (out_pickle, wav_range), 'rb')
	pop = pickle.load(pickleFile)
	pickleFile.close()
	
	grid_length = 40
		# Find lines:
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
	
	## Alternative method for loading pop
	# line_dir = {}
	# uncert_dir = {}
	# for l in lines:
	# 	ab, uncert = D.absorption_line(l, uncert=True)
	# 	line_dir[l] = ab
	# 	uncert_dir[l] = uncert
	# pop = population(line_dir, uncert_dir, previous=True, galaxy=galaxy)

	s = [grid_length,grid_length,grid_length]

	# Produce and Save Plots
	f, ax_array = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, sharex='col', 
		sharey='row')
	for i, line in enumerate(lines):
		print "    " + line

		# import corner

		ab_line = []
		f2,ax2=plt.subplots(1,3)
		for j in range(D.number_of_bins):
			# fig = corner.corner(pop.samples[j,:,:], labels=["age", "metalicity", "alpha"])#, truths=[m_true, b_true, np.log(f_true)])
			# plt.show()
			ax2[0].hist(pop.samples[j,:,0])
			ax2[1].hist(pop.samples[j,:,1])
			ax2[2].hist(pop.samples[j,:,2])

			a = pop._age[nearest(pop._age, pop.age[j])]
			b = pop._metallicity[nearest(pop._metallicity, pop.metallicity[j])]
			c = pop._alpha[nearest(pop._alpha, pop.alpha[j])]
			ab_line.append(pop.interp[line](a,b,c))
		plt.show()

		ab_line = np.array(ab_line)

		abmin, abmax = set_lims(D.absorption_line(line))

		ax_array[int(np.floor(i/2)),i%2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, ab_line, vmin=abmin, vmax=abmax,
			nodots=True, colorbar=True, label='Index strength ('+r'$\AA$'+')', 
			title=line, ax=ax_array[int(np.floor(i/2)),i%2], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)#, signal_noise=D.SNRatio)

		


	f.set_size_inches(8.5,int(np.ceil(len(lines)/2.0))*1.8)

	saveTo = "%s/stellar_pop_self_consis_%s.pdf" % (out_plots, wav_range)
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
	stellar_pop_self_consis('ic1531', wav_range='4200-', vLimit=2)

