## ==================================================================
## 		Find absorption line values for given apature
## ==================================================================
## warrenj 20170206 A routine to return absorption line strengths for a 
##	given appature

import numpy as np
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
from errors2 import determine_goodpixels, remove_anomalies, get_stellar_templates, \
	get_emission_templates
import ppxf_util as util
from ppxf import ppxf
from scipy import ndimage # for gaussian blur
from absorption import absorption

c = 299792.458 # km/s

stellar_moments = 4 
gas_moments = 2
quiet = True
gas = 0


def compare_absorption(galaxy):

	from tools import funccontains
## ----------============== Load galaxy info ================---------
	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]

	f = fits.open('%s/Data/vimos/cubes/%s.cube.combined.corr.fits' % (cc.base_dir, 
		galaxy))
	lam = np.arange(f[0].header['NAXIS3'])*f[0].header['CDELT3'] + f[0].header['CRVAL3']
## ----------============ Find data within slit =============---------

	x = (np.arange(f[0].header['NAXIS1']) * f[0].header['CDELT1']).repeat(
		f[0].header['NAXIS2'])
	y = np.tile(np.arange(f[0].header['NAXIS2']) * f[0].header['CDELT2'],
		f[0].header['NAXIS1'])

	slit_x, slit_y = get_slit(galaxy, 10, 0.9, 30)

	frac_image = np.zeros((f[0].header['NAXIS1'], f[0].header['NAXIS2']))
	# frac = funccontains(slit, (slit_x,slit_y), x=x, y=y).astype(int)#, fraction=True)
	frac = funccontains(slit, (slit_x,slit_y), x=x, y=y, fraction=True)
	frac_image[np.arange(f[0].header['NAXIS1']).repeat(f[0].header['NAXIS2']),np.tile(
		np.arange(f[0].header['NAXIS2']),f[0].header['NAXIS1'])] = frac

	cube = f[0].data
	cube[f[3].data==1] = 0
	noise_cube = f[1].data
	noise_cube[f[3].data==1] = 0.000000001
	cube[~np.isfinite(noise_cube)] = 0
	noise_cube[~np.isfinite(noise_cube)] = 0.000000001


	gal_spec = np.einsum('ijk,jk->i', cube, frac_image)
	gal_noise = np.sqrt(np.einsum('ijk,jk->i', noise_cube**2, frac_image**2))

	gal_spec, lam, cut = remove_anomalies(gal_spec, window=201, repeats=3, 
		lam=lam, set_range=np.array([4200,10000]), return_cuts=True)
	gal_noise = gal_noise[cut]
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
## ----------================= Templates ====================---------
	FWHM_gal = 2.5 # VIMOS documentation (and fits header)
	FWHM_gal = FWHM_gal/(1+z) # Adjust resolution in Angstrom

	stellar_templates = get_stellar_templates(galaxy, FWHM_gal)
	velscale = stellar_templates.velscale

	e_templates = get_emission_templates(gas, lamRange, 
		stellar_templates.logLam_template, FWHM_gal)

	if gas:
		templates = np.column_stack((stellar_templates.templates, e_templates.templates))
	else:
		templates = stellar_templates.templates
	component = [0]*len(stellar_templates.templatesToUse) + e_templates.component
	templatesToUse = np.append(stellar_templates.templatesToUse, 
		e_templates.templatesToUse)
	element = ['stellar'] + e_templates.element

	start = [[vel, sig]] * (max(component) + 1)
	moments = [stellar_moments] + [gas_moments] * max(component)
## ----------============== Final calibrations ==============---------
	## smooth spectrum to fit with templates resolution
	if FWHM_gal < stellar_templates.FWHM_tem:
		sigma = stellar_templates.FWHM_dif/2.355/f[0].header['CDELT3']
		gal_spec = ndimage.gaussian_filter1d(gal_spec, sigma)
		gal_noise = np.sqrt(ndimage.gaussian_filter1d(gal_noise**2, sigma))
	
	## rebin spectrum logarthmically
	gal_spec_log, logLam_bin, _ = util.log_rebin(lamRange, gal_spec, velscale=velscale)
	gal_noise_log, logLam_bin, _ = util.log_rebin(lamRange, gal_noise**2, 
		velscale=velscale)
	gal_noise_log = np.sqrt(gal_noise_log)

	gal_noise_log = gal_noise_log + 0.0000000000001



	dv = (stellar_templates.logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodPixels = determine_goodpixels(logLam_bin,stellar_templates.lamRange_template,
		vel, z, gas=gas!=0) 
	lambdaq = np.exp(logLam_bin)
## ----------=================== pPXF =======================---------
	pp = ppxf(templates, gal_spec_log, gal_noise_log, velscale, start, 
		goodpixels=goodPixels, moments=moments, degree=-1, vsyst=dv, 
		component=component, lam=lambdaq, plot=True, 
		quiet=quiet, mdegree=10)
## ----------============== Absorption Line =================---------
	
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
	results = {}
	uncert = {}
	for line in lines:
		### NEED TO REMOVE EMISSION LINES
		ab, uncert = absorbsion(line, lam, gal_spec_log, unc_lam=stellar_templates.wav,
			unc_spec=stellar_templates.templates.dot(
			pp.weights[:stellar_templates.ntemp]), conv_spec=pp.bestfit, 
			noise=gal_noise_log)

		result[line] = ab
		uncert[line] = uncert

	for line in lines:
		print '%s: %f +/- %f' % (line, result[line], round(uncert[line],2))









## is actually an arbitary quadrilateral
## NB: not sure how this would perform with a concave shape
def slit(x, args):
	from tools import length as len
	x_points,y_points = args
	result = np.zeros((len(x),2))

	x_points = np.asarray(x_points).astype(float)
	y_points = np.asarray(y_points).astype(float)
	x = np.asarray(x).astype(float)

	# Bin is inside slit
	above = x_points < np.outer(x, np.ones(4))
	inside = (0 < np.sum(above, axis=1)) * (4 > np.sum(above, axis=1))

	# A Reference Index array
	ref = np.outer(np.ones(len(x)),np.arange(4))
	ref[above] = np.nan
	ref[~inside,:] = np.nan # If outside of the quadrilateral.


	# Lines which are interected by vertical lines at x
	# NB: use_lines is only transfered to int type later in order to retain the nan's for
	#	as long as possible
	use_lines = np.dstack([np.vstack([np.nanmax(ref, axis=1), np.nanmin(ref, axis=1)]).T,
		np.vstack([np.nanmax(ref, axis=1)+1, np.nanmin(ref,axis=1)-1]).T])
	use_lines[use_lines==4] = 0 # Cyclic boundary
	use_lines[use_lines==-1] = 3

	# Calc
	x0 = np.full((len(x),2),np.nan)
	x0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(
		x_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape(
		(np.sum(inside),2))
	x1 = np.full((len(x),2),np.nan)
	x1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(
		x_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape(
		(np.sum(inside),2))

	y0 = np.full((len(x),2),np.nan)
	y0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(
		y_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape(
		(np.sum(inside),2))
	y1 = np.full((len(x),2),np.nan)
	y1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(
		y_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape(
		(np.sum(inside),2))

	result = y0 + (y1 - y0)/(x1 - x0) * (np.outer(x, np.ones(2)) - x0)

	return result[:,0], result[:,1]



# warrenj 20170206 Routine to get the coordinates of a slit from hight (h), width(w) 
#	and postion angle (pa). NB: pa must be in degrees
def get_slit(galaxy,h,w,pa):
	pa = np.radians(pa)
	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	res=0.67
	x_cent, y_cent = x_gals[i_gal]*res, y_gals[i_gal]*res

	# NB: RA (x) is increasing right to left
	x = [x_cent + w/2, x_cent - w/2, x_cent - w/2, x_cent + w/2]
	y = [y_cent + h/2, y_cent + h/2, y_cent - h/2, y_cent - h/2]
	# Apply pa (angle N to E or top to left)
	rotx = x_cent + np.cos(pa) * (x - x_cent) - np.sin(pa) * (y - y_cent)
	roty = y_cent + np.sin(pa) * (x - x_cent) + np.cos(pa) * (y - y_cent)

	return rotx,roty




if __name__ == '__main__':
	# a=slit([3.,4.,10.,1.5],([1.,3.,4.,2.],[3.,5.,3.,2.]))
	# print a
	galaxy = 'ngc3557'
	# galaxy = 'ic1459'
	compare_absorption(galaxy)