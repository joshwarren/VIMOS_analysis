## ==================================================================
## 		Measure Lick indices for absorption lines
## ==================================================================
import cPickle as pickle
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from checkcomp import checkcomp
cc = checkcomp()

def absorption(galaxy, line_name, D):

	c = 299792.458 # speed of light in km/s
# ------------====== Read absorption line file  =========----------

	ab_file = '/home/HOME/Documents/useful_files/ab_linelist.dat'
	i1, i2, b1, b2, r1, r2, units = np.genfromtxt(ab_file, unpack=True, 
		usecols=(1,2,3,4,5,6,7), skip_header=2, skip_footer=2)
	lines = np.genfromtxt(ab_file, unpack=True, dtype=str, usecols=(8), 
		skip_header=2, skip_footer=2)
	line = np.where(lines==line_name)[0][0]
	index = [i1[line],i2[line]]
	blue = [b1[line],b2[line]]
	red =[r1[line],r2[line]]
# ------------========== Find line strength  ============----------

	line_strength = []
	for bin in D.bin:
		lam = bin.lam/(1+bin.components['stellar'].vel/c) # move to rest frame
		# remove dispersion
		# pixel locations
		indexx = [np.abs(lam-index[0]).argmin(), np.abs(lam-index[1]).argmin()]
		redx = [np.abs(lam-red[0]).argmin(), np.abs(lam-red[1]).argmin()]
		bluex = [np.abs(lam-blue[0]).argmin(), np.abs(lam-blue[1]).argmin()]
		# average flux in side bands
		F_red = np.mean(bin.continuum[redx[0]:redx[1]])
		F_blue = np.mean(bin.continuum[bluex[0]:bluex[1]])
		# Gradient of staight line representing continuum
		m = 2 * (F_red - F_blue)/(redx[0]+redx[1] - bluex[0]-bluex[1])
		# Staight line representing continuum
		F_c = m *np.arange(len(bin.continuum)) + F_red - m * (redx[0]+redx[1])/2
		# Indice value
		I = np.trapz(1-bin.continuum[indexx[0]:indexx[1]]/F_c[indexx[0]:indexx[1]], 
			x=lam[indexx[0]:indexx[1]])
		line_strength.append(I)

	return line_strength






##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
	galaxy = 'ic1459'
	line = 'Fe5015'

	wav_range="4200-"
	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""
	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_pickle = '%s/pickled' % (output)
	pickleFile = open("%s/dataObj_%s.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	s = absorption(galaxy, line, D)

	plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, s, nodots=True,
		colorbar=True, galaxy = galaxy.upper(), save='/home/HOME/test.png')