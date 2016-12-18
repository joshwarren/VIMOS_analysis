## ==================================================================
## 		Measure Lick indices for absorption lines
## ==================================================================
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from tools import length as len
from checkcomp import checkcomp
cc = checkcomp()

c = 299792.458 # speed of light in km/s

##############################################################################
def calc(lam, spectrum, index, blue, red, uncert=None):
# ------------========== Find line strength  ============----------
	# pixel locations
	if blue[0] < lam[0] or red[1] > lam[-1]:
		line_strength = np.nan
		uncert_out = np.nan
	else:
		indexx = [np.abs(lam-index[0]).argmin(), np.abs(lam-index[1]).argmin()]
		redx = [np.abs(lam-red[0]).argmin(), np.abs(lam-red[1]).argmin()]
		bluex = [np.abs(lam-blue[0]).argmin(), np.abs(lam-blue[1]).argmin()]
		# average flux in side bands
		F_red = np.nanmean(spectrum[redx[0]:redx[1]])
		F_blue = np.nanmean(spectrum[bluex[0]:bluex[1]])
		# Gradient of staight line representing continuum
		m = 2 * (F_red - F_blue)/(red[0]+red[1] - blue[0]-blue[1])
		# Staight line representing continuum
		F_c = m*lam + F_red - m * (red[0]+red[1])/2
		# Indice value
		line_strength = np.trapz(1-spectrum[indexx[0]:indexx[1]]/
			F_c[indexx[0]:indexx[1]], x=lam[indexx[0]:indexx[1]])
		if uncert is not None:
			# uncertainty in continuum
			b = np.sqrt(np.nanmean(uncert[bluex[0]:bluex[1]]**2))/np.sqrt(bluex[1]-bluex[0])
			r = np.sqrt(np.nanmean(uncert[redx[0]:redx[1]]**2))/np.sqrt(redx[1]-redx[0])
			F_c_uncert = ((2*lam-blue[1]-blue[0])/(red[0]+red[1]-blue[0]-blue[1])
				)*np.sqrt(b**2+r**2)

			# Uncertainty in absorption line
			integ_var = np.sqrt(F_c_uncert**2+uncert**2)
			uncert_out = np.sqrt(((lam[indexx[0]+1]-lam[indexx[0]])*integ_var[indexx[0]])**2
				+ ((lam[indexx[1]]-lam[indexx[1]-1])*integ_var[indexx[1]])**2
				+ np.sum(((lam[indexx[0]+2:indexx[1]]-lam[indexx[0]:indexx[1]-2])
					*integ_var[indexx[0]+1:indexx[1]-1])**2))/2

	if uncert is not None:
		return line_strength, uncert_out
	else:
		return line_strength
##############################################################################


##############################################################################
def absorption(line_name, lam, spec, unc_lam=None, unc_spec=None, conv_spec=None, 
	noise=None):
# Keywords:
# line_name: 			Name of line to be found
# lam:					Wavelngth array of observed spectrum and conv_spec if supplied 
# spec:					Observed spectrum
# unc_lam:		None	Wavelength of unconvolved spectrum - from the templates used by pPXF
# unc_spec:		None	Unconcolved spectum
# conv_spec:	None	Convolved spectrum i.e. pPXF bestfit
# noise:		None	Observed noise from reduction pipeline
##############################################################################
	if len(unc_lam) != len(unc_spec):
		raise ValueError('Length of unconvoled spectrum and corresponding '+\
			'wavelength should be the same.')
	if conv_spec is not None:
		if len(lam) != len(conv_spec) != len(spec):
			raise ValueError('Length of spectrum (bestfit and observed) and '+\
				'corresponding wavelength should be the same.')
	if noise is not None:
		if len(noise) != len(spec):
			raise ValueError('Length of noise and observed spectrum '+\
				'should be the same.')
# ------------====== Read absorption line file ==========----------

	ab_file = '%s/Documents/useful_files/ab_linelist.dat' % (cc.home_dir)
	i1, i2, b1, b2, r1, r2, units = np.genfromtxt(ab_file, unpack=True, 
		usecols=(1,2,3,4,5,6,7), skip_header=2, skip_footer=2)
	lines = np.genfromtxt(ab_file, unpack=True, dtype=str, usecols=(8), 
		skip_header=2, skip_footer=2)
	line = np.where(lines==line_name)[0][0]
	index = [i1[line],i2[line]]
	blue = [b1[line],b2[line]]
	red =[r1[line],r2[line]]

# ------------========= Find line strenghts ==========----------
	if unc_lam is not None and unc_spec is not None and conv_spec is not None:
		# Line strength of unconvolved (/convolved to LICK resolution) spectrum.
		sig_pix = 200*np.mean(index)/c/(unc_lam[1]-unc_lam[0])
		lick_spec = gaussian_filter1d(unc_spec, sig_pix)
		line_strength_uncon = calc(unc_lam, lick_spec, index, blue, red)
		# Line strength of convolved spectrum (bestfit - emission lines)
		line_strength_con = calc(lam, conv_spec, index, blue, red)
		# LOSVD correction (From SAURON VI: Section 4.2.2)
		corr = line_strength_uncon/line_strength_con
	else:
		corr = 1
	if noise is not None:
		line_strength, uncert = calc(lam, spec, index, blue, red, uncert=noise)
		line_strength *= corr
		uncert *= corr
		return line_strength, uncert
	else:
		line_strength = corr * calc(lam, spec, index, blue, red)
		return line_strength
##############################################################################



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
	import cPickle as pickle
	from plot_velfield_nointerp import plot_velfield_nointerp


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
	pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	s, uncert = D.absorption_line(line, uncert=True)

	plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, s, nodots=True,
		colorbar=True, galaxy = galaxy.upper(), save='/home/HOME/test.png')