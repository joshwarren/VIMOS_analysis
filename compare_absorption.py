## ==================================================================
## 		Find absorption line values for given apature
## ==================================================================
## warrenj 20170206 A routine to return absorption line strengths for a 
##	given appature

import numpy as np
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
from errors2 import set_lines, use_templates, determine_goodpixels, \
	remove_anomalies, get_stellar_templates

def compare_absorption(galaxy):
	from tools import funccontains

	f = fits.open('%s/Data/vimos/cubes/%s.cube.combined.corr.fits' % (cc.base_dir, galaxy))
	
	x = (np.arange(f[0].header['NAXIS1']) - f[0].header['CRPIX1']) * f[0].header['CDELT1'] \
		+ f[0].header['CRVAL1']

	y = (np.arange(f[0].header['NAXIS2']) - f[0].header['CRPIX2']) * f[0].header['CDELT2'] \
		+ f[0].header['CRVAL2']

	slit_x, slit_y = get_slit(galaxy, 5, 0.2, 30)

	frac_image = np.zeros((f[0].header['NAXIS1'], f[0].header['NAXIS2']))
	# frac = funccontains(slit, (slit_x,slit_y), x=x, y=y, fraction=True)
	frac_image[np.arange(f[0].header['NAXIS1']),np.arange(f[0].header['NAXIS2'])] #= frac

	gal_spec = np.einsum('ijk,jk->i', f[0].data, frac_image)
	gal_noise = np.sqrt(np.einsum('ijk,jk->i', f[1].data**2, frac_image**2))

	lam = np.arange(f[0].header['NAXIS3'])*f[0].header['CDELT3'] + f[0].header['CRVAL3']
	gal_spec, lam, cut = remove_anomalies(gal_spec, window=201, repeats=3, 
		lam=lam, set_range=np.array([4200,10000]), return_cuts=True)
	gal_noise = gal_noise[cut]



	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]

	lamRange = np.array([lam[0],lam[-1]])/(1+z)


	FWHM_gal = 2.5 # VIMOS documentation (and fits header)
	templates, logLam_template, templatesToUse = get_stellar_templates(galaxy, FWHM_gal)
	gas = 3
	e_templates, e_component, e_templatesToUse = emission_templates(gas, lamRange, 
		logLam_template, templates, FWHM_gal)

	templates.append(e_templates)
	componet = [0].append(e_component)
	templatesToUse.append(e_templatesToUse)

	# goodPixels = determine_goodpixels(_,lamRange_template,
	# 	vel, z, gas=True)












def emission_templates(gas, lamRange, logLam_template, templates, FWHM_gal, quiet=True):
	import ppxf_util as util
## ----------===============================================---------
## ----------=============== Emission lines ================---------
## ----------===============================================---------
	# moments = stellar_moments
		# moments = [stellar_moments]
		# start_sav = start
	element = []
	component = []
	templatesToUse = []


	## ----------============ All lines together ===============---------
	if gas == 1:
		emission_lines, line_name, line_wav = util.emission_lines(
			logLam_template, lamRange, FWHM_gal, quiet=quiet)

		templatesToUse = np.append(templatesToUse, line_name)

		component = component + [1]*len(line_name)
		templates = np.column_stack((templates, emission_lines))
	   
		# start = [start_sav,start_sav]
		# moments.append(gas_moments)
		element.append('gas')
	## ----------=============== SF and shocks lines ==============---------
	if gas == 2:
		emission_lines, line_name, line_wav = util.emission_lines(
			logLam_template, lamRange, FWHM_gal, quiet=quiet)

		for i in range(len(line_name)):
			

			if 'H' in line_name[i]:
				templatesToUse = np.append(templatesToUse, line_name[i])
				templates = np.column_stack((templates, emission_lines[:,i]))
				component = component + [1]
				element.append['SF']
			else:
				templatesToUse = np.append(templatesToUse, line_name[i])
				templates = np.column_stack((templates, emission_lines[:,i]))
				component = component + [2] 
				element.append['shocks']      

		# start = [start, start_sav, start_sav]
		# moments = [stellar_moments, gas_moments, gas_moments]
	## ----------=========== All lines inderpendantly ==============---------
	if gas == 3:
		emission_lines, line_name, line_wav = util.emission_lines(
			logLam_template, lamRange, FWHM_gal, quiet=quiet)

		aph_lin = np.sort(line_name)

		for i in range(len(line_name)):
			templatesToUse = np.append(templatesToUse, line_name[i])

			# line listed alphabetically
			component = component + [np.where(line_name[i] == aph_lin)[0][0]+1]
			templates = np.column_stack((templates, emission_lines[:,i]))
			# moments.append(gas_moments)
			element.append(aph_lin[i])
		# Bit of a fudge for start (limits ability to set different start for gas)
		# start = [start_sav]*(len(line_name)+1)
		

	return templates, component, templatesToUse







## is actually an arbitary quadrilateral
## NB: not sure how this would perform with a concave shape
def slit(x, args):
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

	# Lines which are interected by vertical lines at x
	use_lines = np.dstack([np.vstack([np.nanmax(ref, axis=1), np.nanmin(ref, axis=1)]).T,
		np.vstack([np.nanmax(ref, axis=1)+1, np.nanmin(ref,axis=1)-1]).T])
	use_lines[use_lines==4] = 0 # Cyclic boundary

	# Calc
	x0 = np.full((len(x),2),np.nan)
	x0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(x_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape((np.sum(inside),2))
	x1 = np.full((len(x),2),np.nan)
	x1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(x_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape((np.sum(inside),2))

	y0 = np.full((len(x),2),np.nan)
	y0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(y_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape((np.sum(inside),2))
	y1 = np.full((len(x),2),np.nan)
	y1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(y_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape((np.sum(inside),2))

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
	x = [x_cent + w/2, x_cent - w/2, x_cent + w/2, x_cent - w/2]
	y = [y_cent + h/2, y_cent + h/2, y_cent - h/2, y_cent - h/2]
	# Apply pa (angle N to E or top to left)
	rotx = x_cent + np.cos(pa) * (x - x_cent) - np.sin(pa) * (y - y_cent)
	roty = y_cent + np.sin(pa) * (x - x_cent) + np.cos(pa) * (y - y_cent)

	return rotx,roty




if __name__ == '__main__':
	# a=slit([3.,4.,10.,1.5],([1.,3.,4.,2.],[3.,5.,3.,2.]))
	# print a
	galaxy = 'ngc3557'
	compare_absorption(galaxy)