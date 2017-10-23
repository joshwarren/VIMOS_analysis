# Routine to show global Mg-sigma relation

import numpy as np # for array handling
# import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
# from scipy import ndimage # for gaussian blur
# from scipy.spatial.distance import cdist
# import os
# import sys
from checkcomp import checkcomp
cc = checkcomp()
if cc.device == -1:
	cc = checkcomp(override='glamdring')
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
# import of matplotlib.pyplot is within errors routine
from errors2 import run_ppxf, set_params, apply_range, get_dataCubeDirectory
from classify import get_R_e
# from sigma_e import in_aperture
from pop import get_absorption
from disk_fit_functions_binned import rebin

# pixel units
def in_aperture(x_cent, y_cent, r, instrument='vimos', os=20):
	if instrument=='vimos':
		frac = np.zeros((40,40))
	elif instrument=='muse':
		frac = np.zeros((150,150))

	s = frac.shape

	sampler = np.zeros((s[0]*os,s[1]*os))
	index = np.array(np.meshgrid(np.arange(s[0]*os)-x_cent*os, 
		np.arange(s[1]*os)-y_cent*os))
	sampler[index[0]**2 + index[1]**2 < (r*os)**2] = 1

	i = np.arange(s[0])
	j = np.arange(s[1])
	i, j = np.meshgrid(i,j)
	i, j = i.flatten(), j.flatten()
	frac[i,j] = rebin(sampler, i*os+os/2, j*os+os/2, statistic='mean')

	return frac


# aperture in arcsec
def mg_sigma(galaxy, aperture=1.0):
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
	params = set_params(reps=10, produce_plot=False, opt='pop')
	
	galaxies = np.array(['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024'])

	i_gal = np.where(galaxies == galaxy)[0][0]

	if cc.device == 'glamdring':
		dir = cc.base_dir
	else:
		dir = '%s/Data/vimos' % (cc.base_dir)


	data_file = dir + "/analysis/galaxies.txt"
	# different data types need to be read separetly
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent_pix = x_cent_gals[i_gal]
	y_cent_pix = y_cent_gals[i_gal]

	
	if aperture == 'R_e':
		ap = get_R_e(galaxy)/header['CDELT1']
	else: ap = aperture

## ----------===============================================---------
## ----------=============== Run analysis  =================---------
## ----------===============================================---------
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = get_dataCubeDirectory(galaxy)
		
	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	# Normalise each spaxel for population pipeline
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	s = galaxy_data.shape
## ----------========== Spatially Integrating =============---------
	frac_in_ap = in_aperture(x_cent_pix, y_cent_pix, ap)
	galaxy_data = np.einsum('ijk,jk->ijk', galaxy_data, frac_in_ap)
	galaxy_noise = np.einsum('ijk,jk->ijk', galaxy_noise**2, 
		frac_in_ap)
	bin_lin = np.nansum(galaxy_data, axis=(1,2))
	bin_lin_noise = np.sqrt(np.nansum(galaxy_noise, axis=(1,2)))
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = apply_range(bin_lin, lam=lam, 
		set_range=params.set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])
	bin_lin_noise = bin_lin_noise[cut]
	print 'here'
	pp = run_ppxf(galaxy, bin_lin, bin_lin_noise, lamRange, CDELT_spec, 
		params)
	print 'here2'

## ----------=============== Find sigma_e  =================---------
	sigma_e = pp.sol[0][1]
	unc_sigma_e = np.std(pp.MCstellar_kin[:,1])

	if aperture == 'R_e':
		area = np.sum(frac_in_ap)*header['CDELT1']*header['CDELT2']
		if area < 0.97 * np.pi * R_e**2:
			R = np.sqrt(area/np.pi)

			sigma_e = sigma_e * (R_e/R)**-0.066
			unc_sigma_e = np.sqrt(unc_sigma_e**2 + 
				((R_e/R)**-0.066 * np.log(R_e/R) * 0.035)**2)

# ## ----------============ Find dynamical mass ===============---------
# 		G = 4.302*10**-6 # kpc (km/s)^2 M_odot^-1 
# 		M = 5.0 * R_e * sigma_e**2/G

## ----------============ Find dynamical mass ===============---------
	mg, mg_uncert = get_absorption(['Mg_b'], pp=pp, instrument='vimos')

	return mg['Mg_b'], mg_uncert['Mg_b'], sigma_e, unc_sigma_e

if __name__=='__main__':
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']

	m = []
	u_m = []
	s = []
	u_s = []
	for g in galaxies:
		mg, mg_uncert, sigma_e, unc_sigma_e = mg_sigma(g)
		m = np.append(m, mg)
		u_m = np.append(u_m, mg_uncert)
		s = np.append(s, sigma_e)
		u_s = np.append(u_s, unc_sigma_e)

	fig, ax = plt.subplots()
	ax.errorbar(np.log10(s), m, xerr=u_s/s/np.log(10), yerr=u_m, fmt='.', 
		color='k')

	params, cov = np.polyfit(np.log10(s), m, 1, 
		w=np.sqrt(1/(u_s/s/np.log(10))**2 + 1/u_m**2), cov=True)

	lims = ax.get_xlim()
	ax.plot(lims, np.poly1d(params)(lims), 'k')
	ax.plot(lims, 2.7*lims - 1.65, '--', label='Ziegler1997')
	ax.set_xlim(lims)
	fig.text(0.15, 0.8, r'grad: %.3f $\pm$ %.3f'%(params[0], 
		np.sqrt(np.diag(cov))[0]), color='k')

	fig.savefig('%s/Data/vimos/analysis/Mg_sigma.png' % (cc.base_dir))

	f = f('%s/Data/vimos/analysis/galaxies_Mg_sigma.txt')
	temp = "{0:12}{1:6}{2:6}{3:6}{4:6}\n"
	with open(f, 'w') as f:
		f.write(temp.format('Galaxy', 'Mg_b', 'e_Mg_b', 'sigma', 'e_sigma'))
		for i in range(len(galaxies)):
			f.write(temp.format(galaxies[i], str(round(m[i],3)), 
				str(round(u_m[i],3)), str(round(s[i],3)), 
				str(round(u_s[i],3))))