## ==================================================================
## 		Plot the absorption indices
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from astropy.io import fits
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from plot_results import set_lims, add_
from checkcomp import checkcomp
cc = checkcomp()
from errors2 import get_dataCubeDirectory
from prefig import Prefig
Prefig()

def plot_absorption(galaxy, D=None, uncert=True, opt='pop', overplot={}):
	# Find lines:
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		#'Mg_1', 'Mg_2', 
		'Mg_b']
	limits = {#'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
		'H_beta':[1.0,2.9], 'Fe5015':[3.5,5.9], 'Mg_b':[3.1,4.7]}

	print 'Absorption lines'


	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots/absorption" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header
	f.close()

	data_file =  "%s/galaxies.txt" % (out_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]
	z_gals, x_cent_gals, y_cent_gals, SN_target_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,4,5,col), dtype='float,int,int,float')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	SN_target=SN_target_gals[i_gal]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])


	for i, line in enumerate(lines):
		print "    " + line

		if uncert:
			ab_line, ab_uncert = D.absorption_line(line, uncert=True)
		else:
			ab_line = D.absorption_line(line)

		abmin, abmax = set_lims(ab_line, positive=True)

		# if line in limits.keys():
		# 	abmin = limits[line][0]
		# 	abmax = limits[line][1]
		

		saveTo = '%s/%s.png' % (out_plots, line)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, D.absorption_line(line), header, vmin=abmin, vmax=abmax,
			nodots=True, colorbar=True, label='Index strength ('+r'$\AA$'+')', 
			title=line, cmap='gnuplot2', redshift=z,
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center, save=saveTo)
		ax.saveTo = saveTo
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		
		if uncert:
			abmin, abmax = set_lims(ab_uncert)

			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, ab_uncert, 
				header, vmin=abmin, vmax=abmax, nodots=True, colorbar=True, 
				label='Index strength ('+r'$\AA$'+')', title=line, cmap='gnuplot2', 
				flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
				signal_noise_target=SN_target, center=center, 
				save='%s/%s_uncert.png' % (out_plots, line))

	return D








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1531', uncert=False)

