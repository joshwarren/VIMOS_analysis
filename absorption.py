## ==================================================================
## 		Measure Lick indices for absorption lines
## ==================================================================
import cPickle as pickle
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from checkcomp import checkcomp
from glob import glob
cc = checkcomp()

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
#		m = 2 * (F_red - F_blue)/(redx[0]+redx[1] - bluex[0]-bluex[1])
		m = 2 * (F_red - F_blue)/(red[0]+red[1] - blue[0]-blue[1])
		# Staight line representing continuum
#		F_c = m *np.arange(len(spectrum)) + F_red - m * (redx[0]+redx[1])/2
		F_c = m*lam + F_red - m * (red[0]+red[1])/2
		# Indice value
		line_strength = np.trapz(1-spectrum[indexx[0]:indexx[1]]/
			F_c[indexx[0]:indexx[1]], x=lam[indexx[0]:indexx[1]])
		if uncert is not None:
			# uncertainty in continuum
			m_uncert = np.sqrt(np.nanmean(uncert[bluex[0]:bluex[1]]**2) +
				np.nanmean(uncert[redx[0]:redx[1]]**2))*\
				2/(red[0]+red[1] - blue[0]-blue[1])

			c_uncert = np.sqrt(np.nanmean(uncert[redx[0]:redx[1]]**2) +
				(m_uncert*(red[0]+red[1])/2)**2)
			F_c_uncert = np.sqrt((lam*m_uncert)**2 + c_uncert**2)
			b = np.sqrt(np.nanmean(uncert[bluex[0]:bluex[1]]**2))/np.sqrt(bluex[1]-bluex[0])
			r = np.sqrt(np.nanmean(uncert[redx[0]:redx[1]]**2))/np.sqrt(redx[1]-redx[0])
			F_c_uncert = ((2*lam-blue[1]-blue[0])/(red[0]+red[1]-blue[0]-blue[1])
				)*np.sqrt(b**2+r**2)
			# print b/F_blue
			# print r/F_red
			# print m_uncert/m
			# print c_uncert/np.abs(F_red - m * (red[0]+red[1])/2)
			# print F_c_uncert[indexx[0]:indexx[1]]/F_c[indexx[0]:indexx[1]]
			# print uncert[indexx[0]:indexx[1]]/spectrum[indexx[0]:indexx[1]]
			# print ''
			# print np.sqrt(b**2+r**2)/(F_red-F_blue)
			# import matplotlib.pyplot as plt
			# plt.close('all') 
			# f,ax=plt.subplots()
			# ax.errorbar(lam, F_c, yerr=F_c_uncert)
			# ax.errorbar(lam,spectrum,yerr=uncert, color='g')
			# ax.errorbar([np.mean(red),np.mean(blue)],[F_red,F_blue],yerr=[r,b],fmt='.')
			# ax.axvline(blue[0],color='b')
			# ax.axvline(red[0],color='r')
			# ax.axvline(index[0],color='g')
			# ax.axvline(blue[1],color='b')
			# ax.axvline(red[1],color='r')
			# ax.axvline(index[1],color='g')
			# ax.set_xlim([blue[0]-5,red[1]+5])
			# ax.set_ylim([0.1,0.25])
			# plt.show()
			# kjshkd

			# F_c_uncert = np.sqrt(
			# 	(np.trapz((spectrum[bluex[0]:bluex[1]] - F_c[bluex[0]:bluex[1]])**2,x=lam[bluex[0]:bluex[1]]) +
			# 	np.trapz((spectrum[redx[0]:redx[1]] - F_c[redx[0]:redx[1]])**2,x=lam[redx[0]:redx[1]]))
			# 	/(redx[1]-redx[0] + bluex[1]-bluex[0]))

			# Uncertainty in absorption line
			integ_var = np.sqrt(F_c_uncert**2+uncert**2)
			uncert_out = np.sqrt(((lam[indexx[0]+1]-lam[indexx[0]])*integ_var[indexx[0]])**2
				+ ((lam[indexx[1]]-lam[indexx[1]-1])*integ_var[indexx[1]])**2
				+ np.sum(((lam[indexx[0]+2:indexx[1]]-lam[indexx[0]:indexx[1]-2])
					*integ_var[indexx[0]+1:indexx[1]-1])**2))/2

			# uncert_out = np.sqrt(np.trapz(np.square(uncert[indexx[0]:indexx[1]]/
			# 	F_c[indexx[0]:indexx[1]]),x=lam[indexx[0]:indexx[1]]))
	if uncert is not None:
		return line_strength, uncert_out
	else:
		return line_strength




def absorption(line_name, D, uncert=False):
	_uncert=uncert
	c = 299792.458 # speed of light in km/s

# ------------====== Templates for unconvolved =======----------
	#files = glob('%s/Data/idl_libraries/ppxf/MILES_library/m0[0-9][0-9][0-9]V' %
	#	(cc.base_dir))
	files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
	wav = np.loadtxt(files[0], usecols=(0,), unpack=True)
	templates = {}
	for template in D.bin[0].temp_weight.keys():
		if template.isdigit():
			templates[template] = np.loadtxt(files[int(template)], usecols=(1,), 
				unpack=True) 

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
	line_map = []
	uncert_map = []
	for i, bin in enumerate(D.bin):
		lam, spec = bin.unconvolved_spectrum(wav, templates)
		# Line strength of unconvolved spectrum.
		line_strength_uncon = calc(lam, spec, index, blue, red)
		# Line strength of convolved spectrum (bestfit - emission lines)
		convolved = bin.bestfit - np.nansum([line.spectrum_nomask for key, line in 
			bin.e_line.iteritems()], axis=0)
		# move observed spectrum to rest frame (z is already accounted for in 
		# errors2.py)
		lam = bin.lam/(1+bin.components['stellar'].vel/c)
		line_strength_con = calc(lam, convolved, index, blue, red)
		# LOSVD correction (From SAURON VI: Section 4.2.2)
		corr = line_strength_uncon/line_strength_con
		if _uncert:
			line_strength, uncert = calc(lam, bin.continuum, index, blue, red, 
				uncert=bin.noise)
			line_strength * corr
			uncert_map.append(uncert*corr)
		else:
			line_strength = corr * calc(lam, bin.continuum, index, blue, red)
		line_map.append(line_strength)

	if _uncert:
		return np.array(line_map), np.array(uncert_map)	
	else:
		return np.array(line_map)









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
	pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	s, uncert = absorption(line, D, uncert=True)
	for i in range(len(s)):
		print s[i], uncert[i]

	plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, s, nodots=True,
		colorbar=True, galaxy = galaxy.upper(), save='/home/HOME/test.png')