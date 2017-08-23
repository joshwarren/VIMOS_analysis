import numpy as np 
from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp
import cPickle as pickle
from astropy.io import fits
from errors2 import get_dataCubeDirectory
import disk_fit_functions as dfn
from sauron_colormap import sauron
from prefig import Prefig 



def fit_disk(galaxy, D=None, opt='kin'):
	Prefig(subplots=(3,2))
	fig, ax = plt.subplots(2,3)

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header
	f.close()

	if D is None:
		pickleFile = open('%s/Data/vimos/analysis/%s/%s/pickled/dataObj.pkl' % (cc.base_dir,
			galaxy,opt))
		D = pickle.load(pickleFile)
		pickleFile.close()

	# Use gas disk
	vel = D.components['Hbeta'].plot['vel'].unbinned
	vel_err = D.components['Hbeta'].plot['vel'].uncert.unbinned

	disk,pars=dfn.disk_fit_exp(vel.copy(),vel_err.copy(),sigclip=3.0,leeway=2.)

	

	ax[1,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['Hbeta'].plot['vel'], header, #vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[1,0],
		title=r'Gas observed velocity (v$_\mathrm{gas}$)')

	ax[1,1].imshow(disk, cmap=sauron)
	ax[1,1].set_title(r'Gas model velocity (v$_\mathrm{g,mod})$')

	plot = D.components['Hbeta'].plot['vel']-D.rebin(disk, flux_weighted=True)
	# vmin,vmax = set_lims(plot)

	ax[1,2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, plot, header, 
		# vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, 
		ax=ax[1,2], title=r'Residuals: v$_\mathrm{gas}$ - v$_\mathrm{g,mod}$')



	# Use stellar disk
	vel = D.components['stellar'].plot['vel'].unbinned
	vel_err = D.components['stellar'].plot['vel'].uncert.unbinned

	disk,pars=dfn.disk_fit_exp(vel.copy(),vel_err.copy(),sigclip=3.0,leeway=2.)

	
	ax[0,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		D.components['stellar'].plot['vel'], header, #vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, ax=ax[0,0],
		title=r'Stellar observed velocity (v$_\mathrm{\ast}$)')

	ax[0,1].imshow(disk, cmap=sauron)
	ax[0,1].set_title(r'Gas model velocity (v$_\mathrm{\ast,mod}$)')


	plot = D.components['Hbeta'].plot['vel']-D.rebin(disk, flux_weighted=True)
	# vmin,vmax = set_lims(plot)

	ax[0,2] = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, plot, header, 
		# vmin=vmin, vmax=vmax, 
		flux_unbinned=D.unbinned_flux, nodots=True, colorbar=True, 
		ax=ax[0,2], title=r'Residuals: v$_\mathrm{gas}$ - v$_\mathrm{\ast,mod}$')

	fig.suptitle('Outflows in %s' % (galaxy))

	fig.savefig('%s/Data/vimos/analysis/%s/%s/plots/outflows.png' % (cc.base_dir, galaxy,
		opt))

	return D


if __name__=='__main__':
	for g in ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100',
		'ngc3557', 'ngc7075', 'pks0718-34']:
		fit_disk(g)