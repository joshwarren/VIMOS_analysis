## ==================================================================
## Propergate uncertainty
## ==================================================================
## warrenj 20150216 Process to progerate the uncertainty using Monty
## Carlo methods to get uncertainty in velocity space.
## warrenj 20170405 Addapted to get sigma_e


import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
from scipy import ndimage # for gaussian blur
from scipy.spatial.distance import cdist
import os
import sys
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
from errors2 import run_ppxf, set_params, apply_range
from ppxf import ppxf
import ppxf_util as util
from classify import get_R_e


#-----------------------------------------------------------------------------
class in_aperture(object):
	def __init__(self, x_cent, y_cent, r, x, y):
		self.x_cent = x_cent
		self.y_cent = y_cent
		self.r = r

		self.x = x
		self.y = y

	@property
	def contains(self):
		import matplotlib.path as mplPath
		aperture = mplPath.Path.circle([self.x_cent, self.y_cent], self.r)
		return aperture.contains_points(zip(self.x,self.y))

	@property
	def fraction(self):
		x_sample= np.linspace(min(self.x)-0.5, max(self.x)+0.5, 
			np.ceil(np.sqrt(len(self.x))*10)).repeat(np.ceil(np.sqrt(len(self.y))*10))
		y_sample= np.tile(np.linspace(min(self.y)-0.5, max(self.y)+0.5, 
			np.ceil(np.sqrt(len(self.y))*10)), int(np.ceil(np.sqrt(len(self.x))*10)))

		xdelt = np.subtract.outer(x_sample,self.x)
		ydelt = np.subtract.outer(y_sample,self.y)
		sample_ownership = np.zeros(len(x_sample))
		# De-vectorised due to MemoryError
		for i in xrange(len(sample_ownership)):
			sample_ownership[i] = np.argmin(xdelt[i,:]**2+ydelt[i,:]**2)

		contained = in_aperture(self.x_cent, self.y_cent, self.r, x_sample, 
			y_sample).contains

		inside, counts_in = np.unique(sample_ownership[contained], return_counts=True)
		total, counts_tot = np.unique(sample_ownership, return_counts=True)

		frac = np.zeros(len(self.x))
		frac[inside.astype(int)] = counts_in
		frac /= counts_tot
		return frac
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def sigma_e(i_gal=None):
	if i_gal is None: i_gal=int(sys.argv[1])
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
	params = set_params()
	params.reps = 3
	
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[i_gal]

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

	R_e = get_R_e(galaxy)

## ----------===============================================---------
## ----------=============== Run analysis  =================---------
## ----------===============================================---------
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = "%s/cubes/%s.cube.combined.corr.fits" % (dir,galaxy)
		
	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	# Normalise each spaxel for population pipeline
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	R_e_pix = R_e/header['CDELT1']
	s = galaxy_data.shape

	rows_to_remove = range(params.discard)
	rows_to_remove.extend([s[1]-1-i for i in range(params.discard)])
	cols_to_remove = range(params.discard)
	cols_to_remove.extend([s[2]-1-i for i in range(params.discard)])

	galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
	galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
	galaxy_noise = np.delete(galaxy_noise, rows_to_remove, axis=1)
	galaxy_noise = np.delete(galaxy_noise, cols_to_remove, axis=2)
	galaxy_badpix = np.delete(galaxy_badpix, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_badpix, cols_to_remove, axis=2)

	s = galaxy_data.shape
## ----------========== Spatially Integrating =============---------
	frac_in_ap = in_aperture(x_cent_pix, y_cent_pix, R_e, np.arange(s[1]).repeat(s[2]), 
		np.tile(np.arange(s[2]),s[1])).fraction.reshape(s[1],s[2])
	galaxy_data = np.einsum('ijk,jk->ijk', galaxy_data, frac_in_ap)
	galaxy_noise = np.einsum('ijk,jk->ijk', galaxy_noise, frac_in_ap)
	bin_lin = np.nansum(galaxy_data,axis=(1,2))
	bin_lin_noise = np.nansum(galaxy_noise**2, axis=(1,2))
	bin_lin_noise = np.sqrt(bin_lin_noise)
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = apply_range(bin_lin, lam=lam, 
		set_range=params.set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])
	bin_lin_noise = bin_lin_noise[cut]

	pp = run_ppxf(galaxy, bin_lin, bin_lin_noise, lamRange, CDELT_spec, params)

## ----------=============== Find sigma_e  =================---------
	sigma_e = pp.sol[0][1]
	unc_sigma_r = np.std(pp.stellar_output[:,1])

	area = np.sum(frac_in_ap)*0.67**2 
	if area < 0.97 * np.pi * R_e**2:
		R = np.sqrt(area/np.pi)

		sigma_e = sigma_e * (R_e/R)**-0.066
		unc_sigma_e = np.sqrt(unc_sigma_r**2 + 
			((R_e/R)**-0.066 * np.log(R_e/R) * 0.035)**2)

## ----------============ Find dynamical mass  ==============---------
	G = 4.302*10**-6 # kpc (km/s)^2 M_odot^-1 
	M = 5.0 * R_e * sigma_e**2/G

## ----------============ Write ouputs to file =============---------
	data_file = "%s/analysis/galaxies_sigma_e.txt" % (dir)
	try:
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		sigma_e_gals, unc_sigma_e_gals, mass_gals = np.loadtxt(data_file, skiprows=1, 
			usecols=(1,2,3), unpack=True, dtype=float)

		if sigma_e_gals.shape == ():
			sigma_e_gals = np.array([sigma_e_gals])
			unc_sigma_e_gals = np.array([unc_sigma_e_gals])
			mass_gals = np.array([mass_gals])
			galaxy_gals = np.array([galaxy_gals])
	except IOError:
		galaxy_gals = np.array([galaxy])
		sigma_e_gals = np.array([sigma_e])
		mass_gals = np.array([M])
		unc_sigma_e_gals = np.array([unc_sigma_e])

	i_gal = np.where(galaxy_gals==galaxy)[0]

	if len(i_gal) == 0:
		galaxy_gals = np.append(galaxy_gals, galaxy)
		sigma_e_gals = np.append(sigma_e_gals, sigma_e)
		unc_sigma_e_gals = np.append(unc_sigma_e_gals, unc_sigma_e)
		mass_gals = np.append(mass_gals, M)
	else:
		sigma_e_gals[i_gal] = sigma_e
		unc_sigma_e_gals[i_gal] = unc_sigma_e
		mass_gals[i_gal] = M


	temp = "{0:12}{1:9}{2:15}{3:9}\n"
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "sigma_e", "uncert sigma_e", "Dyn mass"))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(sigma_e_gals[i],4)),
				str(round(unc_sigma_e_gals[i],4)), '{:0.2e}'.format(mass_gals[i])))



##############################################################################


# Use of plot_results.py


if __name__ == '__main__':
	if len(sys.argv)<2:
		for i in range(10):
			sigma_e(i) 
	else: sigma_e()




