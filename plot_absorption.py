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

def plot_absorption(galaxy, wav_range="", vLimit=0, D=None, uncert=True):
	# Find lines:
	#lines = ['Fe5015', 'H_beta', 'Ca4455', 'Mg_b']
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		#'Mg_1', 'Mg_2', 
		'Mg_b']
	limits = {#'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
		'H_beta':[1.0,2.9], 'Fe5015':[3.5,5.9], 'Mg_b':[3.1,4.7]}

	print 'Absorption lines'

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
		D = pickle.load(pickleFile)
		pickleFile.close()

	# Set up figure and subplots
	f, ax_array = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, sharex='col', 
		sharey='row')
	if uncert:
		f_uncert, ax_array_uncert = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, 
			sharex='col', sharey='row')

	for i, line in enumerate(lines):
		print "    " + line

		if uncert:
			ab_line, ab_uncert = D.absorption_line(line, uncert=True)
		else:
			ab_line = D.absorption_line(line)
		
		# Remove big outliers 
		ab_line = D.absorption_line(line)
		std = np.nanstd(ab_line)
		mean = np.nanmedian(ab_line)
		ab_line = ab_line[ab_line<mean+3*std]
		ab_line = ab_line[ab_line>mean-3*std]
		std = np.nanstd(ab_line)
		mean = np.nanmedian(ab_line)

		a_sorted = np.array(sorted(np.unique(D.absorption_line(line))))
		a_sorted = a_sorted[~np.isnan(a_sorted)]
		if len(a_sorted) < 2*vLimit:
			a_sorted= np.array(sorted(D.absorption_line(line)))
		abmin, abmax = a_sorted[vLimit], a_sorted[-vLimit-1]
		abmax = min([abmax, mean + std])
		abmin = max([abmin, mean - std])

		if line in limits.keys():
			abmin = limits[line][0]
			abmax = limits[line][1]

		ax_array[np.floor(i/2),i%2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
			D.yBar, D.absorption_line(line), vmin=abmin, vmax=abmax,
			nodots=True, colorbar=True, label='log Index strength ('+r'$\AA$'+')', 
			title=line, ax=ax_array[np.floor(i/2),i%2], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux)

		if uncert:
			a_sorted = np.array(sorted(np.unique(ab_uncert)))
			a_sorted = a_sorted[~np.isnan(a_sorted)]
			if len(a_sorted) < 2*vLimit:
				a_sorted= np.array(sorted(ab_uncert))
			abmin, abmax = a_sorted[vLimit], a_sorted[-vLimit-1]

			ax_array_uncert[np.floor(i/2),i%2] = plot_velfield_nointerp(D.x, D.y, 
				D.bin_num, D.xBar, D.yBar, ab_uncert, vmin=abmin, vmax=abmax,
				nodots=True, colorbar=True, label='log Index strength ('+r'$\AA$'+')', 
				title=line, ax=ax_array_uncert[np.floor(i/2),i%2], cmap='gnuplot2', 
				flux_unbinned=D.unbinned_flux)


	f.set_size_inches(8.5,int(np.ceil(len(lines)/2.0))*1.8)

	saveTo = "%s/absorption_%s.pdf" % (out_plots, wav_range)
	f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	f.suptitle(galaxy.upper())
	f.savefig(saveTo, bbox_inches="tight")


	if uncert:
		f_uncert.set_size_inches(8.5,int(np.ceil(len(lines)/2.0))*1.8)

		saveTo = "%s/absorption_uncert_%s.pdf" % (out_plots, wav_range)
		f_uncert.tight_layout()
		ax_array_uncert[0,1].set_xlabel('')
		ax_array_uncert[0,0].set_xlabel('')
		ax_array_uncert[0,1].set_ylabel('')
		ax_array_uncert[1,1].set_ylabel('')
		f_uncert.suptitle(galaxy.upper() + ' Uncertainties')
		f_uncert.savefig(saveTo, bbox_inches="tight")

	return D








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1459', wav_range='4200-', vLimit=2)

