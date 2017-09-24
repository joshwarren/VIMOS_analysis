import numpy as np 
from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp
import cPickle as pickle
from plot_results_muse import set_lims
from astropy.io import fits
from errors2 import get_dataCubeDirectory
import disk_fit_functions_binned as dfn
from sauron_colormap import sauron
from prefig import Prefig 



def fit_disk(galaxy, D=None, opt='kin'):
	leeway  = 1.
	sigclip = None#3.

	pa = {'ngc0612':136.238, 'ic1459':None}
	galaxy_gals, pa_gals = np.loadtxt('%s/Data/vimos/analysis/galaxies2.txt' % 
		(cc.base_dir), usecols=(0,3), skiprows=1, unpack=True, dtype=str)
	i_gal = np.where(galaxy_gals == galaxy)[0][0]
	pa = float(pa_gals[i_gal])

	Prefig(subplots=(4,2))
	fig, ax = plt.subplots(2,4)

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header
	f.close()

	if D is None:
		pickleFile = open('%s/Data/vimos/analysis/%s/%s/pickled/dataObj.pkl' % (
			cc.base_dir, galaxy,opt))
		D = pickle.load(pickleFile)
		pickleFile.close()
	D.sauron = True

	# Use gas disk
	vel = D.components['[OIII]5007d'].plot['vel']
	vel_err = D.components['[OIII]5007d'].plot['vel'].uncert
	# vel_err = vel*0+1

	disk,pars=dfn.disk_fit_exp(D.xBar,D.yBar,vel.copy(),vel_err.copy(),leeway=leeway, 
		sigclip=sigclip, grid_length=40, pa=pa)

	ax[1,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['[OIII]5007d'].plot['vel'], header, #vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[1,0],
		title=r'Gas observed velocity (v$_\mathrm{gas}$)')

	vmin,vmax = set_lims(D.components['[OIII]5007d'].plot['vel'].uncert, positive=True)
	ax[1,1] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['[OIII]5007d'].plot['vel'].uncert, header, vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[1,1],
		title=r'Gas observed velocity Uncertainty')

	ax[1,2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, disk, header, 
		#vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[1,2],
		title=r'Gas model velocity (v$_\mathrm{g,mod})$')

	plot = D.components['[OIII]5007d'].plot['vel'] - disk
	vmin,vmax = set_lims(plot, symmetric=True)

	ax[1,3] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, plot, header, 
		vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, 
		ax=ax[1,3], title=r'Residuals: v$_\mathrm{gas}$ - v$_\mathrm{g,mod}$')


	# Use stellar disk
	vel = D.components['stellar'].plot['vel']
	vel_err = D.components['stellar'].plot['vel'].uncert
	# vel_err = vel*0+1

	disk,pars=dfn.disk_fit_exp(D.xBar, D.yBar, vel.copy(),vel_err.copy(),leeway=leeway, 
		sigclip=sigclip, grid_length=40, pa=pa)

	
	ax[0,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['stellar'].plot['vel'], header, #vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[0,0],
		title=r'Stellar observed velocity (v$_\mathrm{\ast}$)')

	vmin,vmax = set_lims(D.components['stellar'].plot['vel'].uncert, positive=True)
	ax[0,1] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['stellar'].plot['vel'].uncert, header, vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[0,1],
		title=r'Stellar observed velocity uncertainty')

	ax[0,2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, disk, header, 
		#vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[0,2],
		title=r'Stellar model velocity (v$_\mathrm{\ast,mod})$')


	plot = D.components['stellar'].plot['vel']-disk
	vmin,vmax = set_lims(plot, symmetric=True)

	ax[0,3] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, plot, header, 
		vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, 
		ax=ax[0,3], title=r'Residuals: v$_\mathrm{\ast}$ - v$_\mathrm{\ast,mod}$',
		signal_noise=D.SNRatio, signal_noise_target=30)

	# plot = D.components['[OIII]5007d'].plot['vel']-disk
	# ax[0,3] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, plot, header, 
	# 	# vmin=vmin, vmax=vmax, 
	# 	flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, 
	# 	ax=ax[0,3], title=r'Residuals: v$_\mathrm{gas}$ - v$_\mathrm{\ast,mod}$')


	fig.suptitle('Outflows in %s' % (galaxy))

	fig.savefig('%s/Data/vimos/analysis/%s/%s/plots/outflows.png' % (cc.base_dir, galaxy,
		opt))

	return D


if __name__=='__main__':
	# fit_disk('ngc0612')
	fit_disk('ic1459')
	# for g in ['ngc0612', 'eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc1399', 
	# 	'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']:
	# 	fit_disk(g)