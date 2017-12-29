# Routine to create a 'bin' which contains the entire image.
import numpy as np
import os
from astropy.io import fits
from errors2 import run_ppxf, apply_range, set_params, get_dataCubeDirectory
from checkcomp import checkcomp
cc = checkcomp()
from Bin import trapz_uncert
from global_mg_sigma import in_aperture
from tools import moving_weighted_average

n_e = 100 # cm^-3
c = 299792.458 # speed of light in km/s
H0 = 70 # km/s/Mpc

# Solar masses
def get_Mass(Ha_flux, D, instrument='vimos'):
	if instrument == 'vimos':
		return 280 * (D/10)**2 * (Ha_flux*10**-15/10**-14) * (1000/n_e) 
	elif instrument == 'muse':
		return 280 * (D/10)**2 * (Ha_flux*10**-20/10**-14) * (1000/n_e) 

def whole_image(galaxy, verbose=True):
	print galaxy
	if cc.device == 'glamdring':
		data_file = "%s/analysis/galaxies.txt" % (cc.base_dir)
	else:
		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1, 4, 5), dtype='float,int,int')
	galaxy_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(0,), dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	D = z*c/H0 # Mpc
	centre = (x_gals[i_gal], y_gals[i_gal])

	limits_file = '%s/Data/vimos/analysis/galaxies_gasMass.txt' %(cc.base_dir)
	galaxy_gals, mass, e_mass = np.loadtxt(limits_file, unpack=True, 
		dtype=str, skiprows=1)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	max_radius = 29 # in spaxels
	radius = float(max_radius)
	Mass_sav = 0
	f = fits.open(get_dataCubeDirectory(galaxy))
	while radius > 2:
		mask = in_aperture(centre[0], centre[1], radius, instrument='vimos')

		spec = f[0].data
		noise = f[1].data

		spec[np.isnan(spec)] = 0
		noise[np.isnan(noise)] = 0

		spec = np.einsum('ijk,jk->i', spec, mask)#/np.sum(mask)
		noise = np.sqrt(np.einsum('ijk,jk->i', noise**2, mask))#/np.sum(mask)

		params = set_params(opt='pop', reps=4, temp_mismatch=True, 
			produce_plot=False, gas=1)

		lam = (np.arange(len(spec)) - (f[0].header['CRPIX3'] - 1)) * \
			f[0].header['CDELT3'] + f[0].header['CRVAL3']
		spec, lam, cut = apply_range(spec, lam=lam, return_cuts=True, 
			set_range=params.set_range)
		lamRange = np.array([lam[0],lam[-1]])
		noise = noise[cut]

		pp = run_ppxf(galaxy, spec, noise, lamRange, f[0].header['CDELT3'], 
			params)

		# pp.noise = np.min([pp.noise, np.abs(pp.galaxy-pp.bestfit)],axis=0)
		residuals = pp.galaxy - pp.bestfit
		_, residuals, _ = moving_weighted_average(pp.lam, residuals, step_size=3., 
			interp=True)
		noise = np.sqrt(residuals**2 + pp.noise**2)

		OIII_spec = pp.matrix[:, pp.templatesToUse==\
			'[OIII]5007d'].flatten()*\
			pp.weights[pp.templatesToUse=='[OIII]5007d']

		Hb_spec = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() * \
			pp.weights[pp.templatesToUse=='Hbeta']
		Hb_flux = np.trapz(Hb_spec, x=pp.lam)
		Ha_flux = 2.86 * Hb_flux

		Hb_spec_uncert = pp.MCgas_uncert_spec[
			pp.templatesToUse[pp.component!=0]=='Hbeta', :].flatten()
		Hb_flux_uncert = trapz_uncert(Hb_spec_uncert, x=pp.lam)
		Ha_flux_uncert = 2.86 * Hb_flux_uncert

		Mass = get_Mass(Ha_flux, D)
		e_Mass = get_Mass(Ha_flux_uncert, D)

		Hb_spec2 = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() \
			/ np.max(pp.matrix[:, pp.templatesToUse=='Hbeta']) \
			* np.median(noise[(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c))
			* (pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))])
		Hb_flux2 = np.trapz(Hb_spec2, x=pp.lam)
		Ha_flux2 = 2.86 * Hb_flux2
		Mass2 = get_Mass(Ha_flux2, D)
		

		if max(OIII_spec)/np.median(noise[
			(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
			(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))]) > 4:

			if max(Hb_spec)/np.median(noise[
				(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c)) *
				(pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:
				
				mass[i_gal] = str(round(np.log10(Mass),4))
				e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
					np.log(10)), 4))
				if verbose:
					print '%s +/- %s log10(Solar Masses)' % (
						mass[i_gal], e_mass[i_gal])

					from ppxf import create_plot
					fig, ax = create_plot(pp).produce 
					ax.set_xlim([4800, 4900])
					ax.legend()
					fig1, ax1 = create_plot(pp).produce 
					ax1.legend()
					import matplotlib.pyplot as plt
					plt.show()

				radius = -1
			else:
				Mass_sav = max(Mass, Mass2, Mass_sav)
				if Mass_sav == Mass2: e_Mass = np.nan

				# if radius == max_radius:
				mass[i_gal] = '<'+str(round(np.log10(Mass_sav),4))
				e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
					np.log(10)), 4))
				if verbose:
					print '<%s +/- %s log10(Solar Masses)' % (
						mass[i_gal], e_mass[i_gal])

					from ppxf import create_plot
					fig, ax = create_plot(pp).produce 
					ax.set_xlim([4800, 4900])
					ax.legend()
					import matplotlib.pyplot as plt
					plt.show()
		else:
			Mass_sav = max(Mass, Mass2, Mass_sav)
			if Mass_sav == Mass2: e_Mass = np.nan

			# if radius == max_radius:
			mass[i_gal] = '<'+str(round(np.log10(Mass_sav),4))
			e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
					np.log(10)), 4))
			if verbose:
				print '<%s +/- %s log10(Solar Masses)' % (
					mass[i_gal], e_mass[i_gal])
				from ppxf import create_plot
				fig, ax = create_plot(pp).produce 
				ax.set_xlim([4800, 4900])
				ax.legend()
				import matplotlib.pyplot as plt
				plt.show()
		radius -= 1

	temp = "{0:12}{1:10}{2:10}\n"
	with open(limits_file, 'w') as l:
		l.write(temp.format('Galaxy', 'Mass', 'e_Mass'))
		for i in range(len(galaxy_gals)):
			l.write(temp.format(galaxy_gals[i], mass[i], e_mass[i]))



if __name__=='__main__':
	galaxies = ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']
	# galaxies = ['pks0718-34']
	# for g in galaxies:
	# 	whole_image(g, verbose=False)
	whole_image('pks0718-34', verbose=False)