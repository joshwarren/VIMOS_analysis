## ==================================================================
## 		Plot the absorption indices
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
#from sauron_colormap2 import sauron2 as sauron
import numpy as np 
from checkcomp import checkcomp
cc = checkcomp()

def plot_absorption(galaxy, wav_range="", vLimit=0):
	# Find lines:
	lines = ['Fe5015', 'H_beta', 'Ca4455', 'Mg_b']

	print 'Absorption lines'

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
	D = pickle.load(pickleFile)
	pickleFile.close()

	# Set up figure and subplots
	f, ax_array = plt.subplots(2, 2, sharex='col', sharey='row')#np.ceil(len(lines)/2.0), sharex='col', sharey='row')
	for i, line in enumerate(lines):
		print "    " + line

		# Remove big outliers 
		ab_line = D.absorption_line(galaxy,line)
		std = np.nanstd(ab_line)
		mean = np.nanmedian(ab_line)
		ab_line = ab_line[ab_line<mean+3*std]
		ab_line = ab_line[ab_line>mean-3*std]
		std = np.nanstd(ab_line)
		mean = np.nanmedian(ab_line)

		a_sorted = np.array(sorted(np.unique(D.absorption_line(galaxy, line))))
		if len(a_sorted) < 2*vLimit:
			a_sorted= np.array(sorted(D.absorption_line(galaxy, line)))
		a_sorted = a_sorted[~np.isnan(a_sorted)]
		abmin, abmax = a_sorted[vLimit], a_sorted[-vLimit-1]
		abmax = min([abmax, mean + std])
		abmin = max([abmin, mean - std])


		ax_array[i%2,np.floor(i/2)] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
			D.yBar, D.absorption_line(galaxy, line), vmin=abmin, vmax=abmax,
			nodots=True, colorbar=True, label='Index strength ('+r'$\AA$'+')', 
			title=line, ax=ax_array[i%2,np.floor(i/2)], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)



	saveTo = "%s/absorption_%s.pdf" % (out_plots, wav_range)
	f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	f.suptitle(galaxy.upper())
	f.savefig(saveTo, bbox_inches="tight")








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1459', wav_range='4200-', vLimit=2)

