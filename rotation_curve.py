# Routine to plot rotation curve
import numpy as np
import matplotlib.pyplot as plt
from checkcomp import checkcomp
cc =checkcomp()
import prefig
prefig.Prefig(transparent=False)
from classify import get_R_e
from plot_results_muse import set_lims


def rotation_curve(galaxy, opt='kin', D=None):
	res = 0.67
	if D is None:
		import cPickle as pickle
		pickleFile = open("%s/Data/vimos/analysis/%s/%s/pickled/dataObj.pkl" % (
			cc.base_dir, galaxy, opt), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	galaxy_file = '%s/Data/vimos/analysis/galaxies.txt' % (cc.base_dir)
	galaxy_file2 = '%s/Data/vimos/analysis/galaxies2.txt' % (cc.base_dir)
	galaxy_gals = np.loadtxt(galaxy_file, skiprows=1, unpack=True, usecols=(0,), 
		dtype=str)
	x_gals, y_gals = np.loadtxt(galaxy_file, skiprows=1, unpack=True, usecols=(4,5), 
		dtype=int)

	i_gal = np.where(galaxy_gals == galaxy)[0][0]
	x_cent, y_cent = x_gals[i_gal], y_gals[i_gal]

	galaxy_gals = np.loadtxt(galaxy_file2, skiprows=1, unpack=True, usecols=(0,), 
		dtype=str)
	pa_gals = np.loadtxt(galaxy_file2, skiprows=1, unpack=True, usecols=(4,))
	i_gal2 = np.where(galaxy_gals == galaxy)[0][0]
	pa = pa_gals[i_gal]





	R_e = get_R_e(galaxy)
	# R = np.sqrt((D.xBar - x_cent)**2 + (D.yBar - y_cent)**2)*res/R_e

	# Distance from axis of rotation to bin
	R = np.abs(np.sin(pa)*D.xBar - np.cos(pa)*D.yBar + x_cent*np.cos(pa) - 
		y_cent*np.sin(pa))*res/R_e




	fig, ax = plt.subplots()
	ax.scatter(R, np.abs(D.components['stellar'].plot['vel']))
	lims = set_lims(np.abs(D.components['stellar'].plot['vel']))
	ax.set_ylim(lims)

	ax.set_xlabel(r'Distance ($R_e$)')
	ax.set_ylabel('Velocity (km/s)')

	fig.savefig('%s/Data/vimos/analysis/%s/%s/plots/rotation_curve.png' % (cc.base_dir,
		galaxy, opt))





if __name__=='__main__':
	# for g in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
	# 	rotation_curve(g, opt='kin')
	# rotation_curve('ngc3100', opt='kin')
	rotation_curve('ngc3557', opt='kin')