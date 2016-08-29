import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
from checkcomp import checkcomp
cc = checkcomp()

def plot_absorption(galaxy, wav_range=""):
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

	lines = ['Fe5015', 'H_beta', 'Ca4455', 'Mg_b']
	f, ax_array = plt.subplots(2, 2)#np.ceil(len(lines)/2), sharex='col', sharey='row')

	for i, line in enumerate(lines):
		print "    " + line
		ax_array[i%2,np.floor(i/2)] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
			D.yBar, D.absorption_line(galaxy, line), #vmin=vmin, vmax=vmax,
			nodots=True, colorbar=True, label='Index strength ('+r'$\AA$'+')', 
			title=line, ax=ax_array[i%2,np.floor(i/2)])
	saveTo = "%s/absorb_%s.pdf" % (out_plots, wav_range)		
	f.savefig(saveTo, bbox_inches="tight")








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1459', wav_range='4200-')

