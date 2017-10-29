## ==================================================================
## 		Plot stellar population maps
## ==================================================================
## warrenj 20170331 Routine to plot the stellar populations found by pop.py
## on Glamdring. 

# import cPickle as pickle
# import matplotlib.pyplot as plt 
import numpy as np 
import os
from astropy.io import fits
from errors2 import get_dataCubeDirectory, apply_range, set_params, run_ppxf
from pop import population
# from plot_results import add_, set_lims
from checkcomp import checkcomp
cc = checkcomp()
# from plot_velfield_nointerp import plot_velfield_nointerp
from classify import get_R_e
# from prefig import Prefig
from global_mg_sigma import in_aperture
# Prefig()

def stellar_pop_grad(galaxy, method='median', opt='pop'):
	print 'Finding population gradients'

	vin_dir = '%s/Data/vimos/analysis/%s/%s/pop' % (cc.base_dir, 
		galaxy, opt)

	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots/pop_distributions" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header

	data_file =  "%s/galaxies.txt" % (out_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(4,5), dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	R_e = get_R_e(galaxy)
	apertures = np.array([R_e/8, R_e/2, R_e])
	str_ap = ['R_e_8', 'R_e_2', 'R_e']

	for i, a in enumerate(apertures):
		params = set_params(reps=100, opt='pop', gas=1, produce_plot=False)

		mask = in_aperture(center[0], center[1], a/header['CDELT1'],
			instrument='vimos')

		spec = np.einsum('ijk,jk->ijk', f[0].data, mask)
		noise = np.einsum('ijk,jk->ijk', f[1].data**2, mask)

		spec = np.nansum(spec, axis=(1,2))
		noise = np.sqrt(np.nansum(noise, axis=(1,2)))

		lam = np.arange(len(spec))*header['CDELT3'] + header['CRVAL3']
		spec, lam, cut = apply_range(spec, lam=lam, 
			set_range=params.set_range, return_cuts=True)
		lamRange = np.array([lam[0],lam[-1]])
		noise = noise[cut]

		pp = run_ppxf(galaxy, spec, noise, lamRange, header['CDELT3'], 
			params)

		pop = population(pp=pp, instrument='vimos', method=method)

		pop.plot_probability_distribution(saveTo='%s/%s.png' % (
			out_plots, str_ap[i]))


		file = '%s/pop_sigma.txt' % (out_dir)
		value = np.loadtxt(file, unpack=True, skiprows=2, dtype=str)
		i_gal = np.where(value[0]==galaxy)[0][0]
		
		value[i*8+1, i_gal] = str(round(pp.sol[0][1], 2))
		value[i*8+2, i_gal] = str(round(np.std(pp.MCstellar_kin[:,1]), 2))
		value[i*8+3, i_gal] = str(round(pop.age, 2))
		value[i*8+4, i_gal] = str(round(pop.unc_age, 2))
		value[i*8+5, i_gal] = str(round(pop.metallicity, 2))
		value[i*8+6, i_gal] = str(round(pop.unc_met, 2))
		value[i*8+7, i_gal] = str(round(pop.alpha, 2))
		value[i*8+8, i_gal] = str(round(pop.unc_alp, 2))


		temp = '{0:12}{1:7}{2:7}{3:7}{4:7}{5:7}{6:7}{7:7}{8:7}{9:7}{10:7}'+\
			'{11:7}{12:7}{13:7}{14:7}{15:7}{16:7}{17:7}{18:7}{19:7}{20:7}'+\
			'{21:7}{22:7}{23:7}{24:7}\n'

		with open(file, 'w') as out:
			out.write(temp.format('Galaxy', 'sig', 'e_sig', 'age', 'e_age', 
				'met', 'e_met', 'alpha', 'e_alpha', 'sig', 'e_sig', 
				'age', 'e_age', 'met', 'e_met', 'alpha', 'e_alpha', 
				'sig', 'e_sig', 'age', 'e_age', 'met', 'e_met', 'alpha', 
				'e_alpha'))
			out.write(temp.format('', 'R_e/8', 'R_e/8', 'R_e/8', 'R_e/8',
				'R_e/8', 'R_e/8', 'R_e/8', 'R_e/8', 'R_e/2', 'R_e/2', 'R_e/2', 
				'R_e/2', 'R_e/2', 'R_e/2', 'R_e/2', 'R_e/2', 'R_e', 'R_e', 
				'R_e', 'R_e', 'R_e', 'R_e', 'R_e', 'R_e'))
			for j in range(len(value[0])):
				out.write(temp.format(*value[:,j]))
		








##############################################################################

# Use of plot_stellar_pop.py

if __name__ == '__main__':
	galaxies = [
			'eso443-g024',
			'ic1459',
			'ic1531', 
			'ic4296',
			'ngc0612',
			'ngc1399',
			'ngc3100',
			'ngc3557',
			'ngc7075',
			'pks0718-34'
			]
	for galaxy in galaxies:
		stellar_pop_grad(galaxy, method='mostlikely')

