# Routine to create a 'bin' which contains the entire image.
import numpy as np
import os
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
from Bin import trapz_uncert
from global_mg_sigma import in_aperture
from tools import moving_weighted_average
from classify import get_R_e

n_e = 100 # cm^-3
c = 299792.458 # speed of light in km/s
H0 = 70 # km/s/Mpc

# Solar masses
def get_Mass(Ha_flux, D, instrument='vimos'):
	if instrument == 'vimos':
		return 280 * (D/10)**2 * (Ha_flux*10**-15/10**-14) * (1000/n_e) 
	elif instrument == 'muse':
		return 280 * (D/10)**2 * (Ha_flux*10**-20/10**-14) * (1000/n_e) 

def whole_image(galaxy, verbose=True, instrument='vimos'):
	print galaxy

	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	galaxy_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(0,), dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	if instrument=='vimos':
		from errors2 import run_ppxf, apply_range, set_params, get_dataCubeDirectory
		res = 0.67 # arcsec/pix
		z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
			usecols=(1, 4, 5), dtype='float,int,int')
		i_gal2 = i_gal
		fits_ext = 0
	elif instrument == 'muse':
		from errors2_muse import run_ppxf, apply_range, set_params, \
			get_dataCubeDirectory
		res = 0.2 # arcsec/pix
		z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,))
		data_file2 = "%s/Data/muse/analysis/galaxies.txt" % (cc.base_dir)
		x_gals, y_gals = np.loadtxt(data_file2, unpack=True, skiprows=1, 
			usecols=(1,2), dtype=int)
		galaxy_gals = np.loadtxt(data_file2, unpack=True, skiprows=1, 
			usecols=(0,), dtype=str)
		i_gal2 = np.where(galaxy_gals==galaxy)[0][0]
		fits_ext = 1
	
	z = z_gals[i_gal]
	D = z*c/H0 # Mpc
	centre = (x_gals[i_gal2], y_gals[i_gal2])

	limits_file = '%s/Data/%s/analysis/galaxies_gasMass.txt' %(cc.base_dir, 
		instrument)
	galaxy_gals, mass, e_mass = np.loadtxt(limits_file, unpack=True, dtype=str, 
		skiprows=1, usecols=(0,1,2))
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	if instrument == 'muse':
		bd, e_bd = np.loadtxt(limits_file, unpack=True, dtype=str, skiprows=1, 
			usecols=(3,4))

		if 'ic' in galaxy:
			f = fits.open(get_dataCubeDirectory(galaxy))			
		elif 'ngc' in galaxy:
			f = fits.open(get_dataCubeDirectory(galaxy)[:-5]+'2.fits')
	
		f[fits_ext].header['CDELT3'] = f[fits_ext].header['CD3_3']

	elif instrument == 'vimos':
		f = fits.open(get_dataCubeDirectory(galaxy))
	R_e = get_R_e(galaxy)/res

	mask = in_aperture(centre[0], centre[1], R_e/2, instrument=instrument)
	spec = f[fits_ext].data
	noise = f[fits_ext+1].data

	spec[np.isnan(spec)] = 0
	noise[np.isnan(noise)] = 0

	spec = np.einsum('ijk,jk->i', spec, mask)
	noise = np.sqrt(np.einsum('ijk,jk->i', noise**2, mask))

	if galaxy == 'ngc1316':
		params = set_params(opt='pop', reps=4, temp_mismatch=True, 
			produce_plot=False, gas=1, set_range_star=np.array([2000, 5800]))
	else:
		params = set_params(opt='pop', reps=4, temp_mismatch=True, 
			produce_plot=False, gas=1)

	lam = (np.arange(len(spec)) - (f[fits_ext].header['CRPIX3'] - 1)) * \
		f[fits_ext].header['CDELT3'] + f[fits_ext].header['CRVAL3']
	spec, lam, cut = apply_range(spec, lam=lam, return_cuts=True, 
		set_range=params.set_range)
	lamRange = np.array([lam[0],lam[-1]])
	noise = noise[cut]

	if instrument == 'vimos':
		pp = run_ppxf(galaxy, spec, noise, lamRange, f[fits_ext].header['CDELT3'], 
			params)
	elif instrument=='muse':
		pp = run_ppxf(galaxy, spec, noise, f[fits_ext].header['CDELT3'],
			 f[fits_ext].header['CRVAL3'], params)

	residuals = pp.galaxy - pp.bestfit
	_, residuals, _ = moving_weighted_average(pp.lam, residuals, step_size=3., 
		interp=True)
	residuals[np.isnan(residuals)] = 0
	noise = np.sqrt(residuals**2 + pp.noise**2)
	# noise = np.array(residuals)
	

	OIII_spec = pp.matrix[:, pp.templatesToUse=='[OIII]5007d'].flatten()*\
		pp.weights[pp.templatesToUse=='[OIII]5007d']
	
	# These results for Hb are used later, even for MUSE
	Hb_spec = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() * \
		pp.weights[pp.templatesToUse=='Hbeta']
	Hb_flux = np.trapz(Hb_spec, x=pp.lam)
	ANR = max(Hb_spec)/np.median(
		noise[(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c)) * (pp.lam > 4861.
		/(1 + (pp.sol[1][0] + 300)/c))]) 
	
	Hb_spec_uncert = pp.MCgas_uncert_spec[
		pp.templatesToUse[pp.component!=0]=='Hbeta', :].flatten()

	Hb_spec_norm = Hb_spec/np.max(Hb_spec)
	Hb_spec_uncert = np.sqrt(Hb_spec_uncert**2 + (noise*Hb_spec_norm)**2)

	Hb_spec_uncert_plus = Hb_spec + Hb_spec_uncert
	Hb_spec_uncert_minus = Hb_spec - Hb_spec_uncert

	Hb_flux_uncert_plus = np.trapz(Hb_spec_uncert_plus, x=pp.lam)
	Hb_flux_uncert_minus = np.trapz(Hb_spec_uncert_minus, x=pp.lam)

	Hb_flux_uncert = np.mean([abs(Hb_flux_uncert_plus - Hb_flux), 
		abs(Hb_flux - Hb_flux_uncert_minus)])

	if instrument == 'vimos':
		Ha_flux = 2.86 * Hb_flux
		Ha_flux_uncert = 2.86 * Hb_flux_uncert
		Ha_flux_uncert_plus = 2.86 * Hb_flux_uncert_plus
		Ha_flux_uncert_minus = 2.86 * Hb_flux_uncert_minus

	elif instrument == 'muse':
		Ha_spec = pp.matrix[:, pp.templatesToUse=='Halpha'].flatten() * \
			pp.weights[pp.templatesToUse=='Halpha']
		Ha_flux = np.trapz(Ha_spec, x=pp.lam)

		Ha_spec_uncert = pp.MCgas_uncert_spec[
			pp.templatesToUse[pp.component!=0]=='Halpha', :].flatten()
		Ha_spec_norm = Ha_spec/np.max(Ha_spec)
		Ha_spec_uncert = np.sqrt(Ha_spec_uncert**2 + (noise*Ha_spec_norm)**2)
		# Ha_flux_uncert = trapz_uncert(Ha_spec_uncert, x=pp.lam)

		Ha_spec_uncert_plus = Ha_spec + Ha_spec_uncert
		Ha_spec_uncert_minus = Ha_spec - Ha_spec_uncert

		Ha_flux_uncert_plus = np.trapz(Ha_spec_uncert_plus, x=pp.lam)
		Ha_flux_uncert_minus = np.trapz(Ha_spec_uncert_minus, x=pp.lam)

		Ha_flux_uncert = np.mean([abs(Ha_flux_uncert_plus - Ha_flux), 
			abs(Ha_flux - Ha_flux_uncert_minus)])

		Hb_ANR = np.array(ANR)
		ANR = max(Ha_spec)/np.median(
			noise[(pp.lam < 6563./(1 + (pp.sol[1][0] - 300)/c)) * (pp.lam > 6563.
			/(1 + (pp.sol[1][0] + 300)/c))]) 

	Mass = get_Mass(Ha_flux, D, instrument=instrument)
	e_Mass = get_Mass(Ha_flux_uncert, D, instrument=instrument)
	# e_Mass = np.mean([
	# 	abs(get_Mass(Ha_flux_uncert_plus, D, instrument=instrument) - Mass),
	# 	abs(Mass - get_Mass(Ha_flux_uncert_minus, D, instrument=instrument))])


	OIII_ANR = max(OIII_spec)/np.median(noise[
		(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
		(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))])

	if OIII_ANR > 4 and ANR > 2.5:
		mass[i_gal] = str(round(np.log10(Mass),4))
		e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
			np.log(10)), 4))

		if verbose:
			print '%s +/- %s log10(Solar Masses)' % (
				mass[i_gal], e_mass[i_gal])

			import matplotlib.pyplot as plt
			fig, ax = plt.subplots(2)
			from ppxf import create_plot
			plot = create_plot(pp)
			plot.ax = ax[0]
			plot.produce 
			ax[0].set_xlim([4840, 4880])
			ax[0].ax2.set_ylim([0,500000])
			ax[0].ax2.plot(pp.lam, residuals, 'k')
			ax[0].legend()
			plot.ax = ax[1]
			plt.show() if 'home' in cc.device else fig.savefig('whole_image2.png')

	else:
		if instrument == 'vimos':
			Hb_spec2 = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() \
				/ np.max(pp.matrix[:, pp.templatesToUse=='Hbeta']) \
				* np.median(noise[(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c))
				* (pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) * 2.5
			Hb_flux2 = np.trapz(Hb_spec2, x=pp.lam)
			Ha_flux2 = 2.86 * Hb_flux2
		elif instrument == 'muse':
			Ha_spec2 = pp.matrix[:, pp.templatesToUse=='Halpha'].flatten() \
				/ np.max(pp.matrix[:, pp.templatesToUse=='Halpha']) \
				* np.median(noise[(pp.lam < 6563./(1 + (pp.sol[1][0] - 300)/c))
				* (pp.lam > 6563./(1 + (pp.sol[1][0] + 300)/c))]) * 2.5
			Ha_flux2 = np.trapz(Ha_spec2, x=pp.lam)
		Mass2 = get_Mass(Ha_flux2, D, instrument=instrument)

		mass[i_gal] = '<'+str(round(np.log10(Mass2),4))
		e_mass[i_gal] = '-'

		if verbose:
			print '<%s +/- %s log10(Solar Masses)' % (
				mass[i_gal], e_mass[i_gal])

			import matplotlib.pyplot as plt
			fig, ax = plt.subplots(2)
			from ppxf import create_plot
			plot = create_plot(pp)
			plot.ax = ax[0]
			plot.produce 
			ax[0].set_xlim([4800, 4900])
			ax[0].ax2.plot(pp.lam, residuals, 'k')
			ax[0].legend()
			plot.ax = ax[1]
			plot.produce
			plt.show() if 'home' in cc.device else fig.savefig('whole_image2.png')

	if instrument == 'muse':
		if OIII_ANR > 4 and ANR > 2.5 and Hb_ANR > 2.5:
			bd[i_gal] = str(round(Ha_flux/Hb_flux, 4))
			e_bd[i_gal] = str(round(Ha_flux/Hb_flux * np.sqrt(
				(Ha_flux_uncert/Ha_flux)**2 + (Hb_flux_uncert/Hb_flux)**2), 4))
			print 'Balmer dec uncert', Ha_flux_uncert/Ha_flux, \
				Hb_flux_uncert/Hb_flux, \
				np.sqrt((Ha_flux_uncert/Ha_flux)**2 + (Hb_flux_uncert/Hb_flux)**2)
		elif OIII_ANR > 4 and ANR > 2.5:
			Hb_spec2 = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() \
				/ np.max(pp.matrix[:, pp.templatesToUse=='Hbeta']) \
				* np.median(noise[(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c))
				* (pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) * 2.5
			Hb_flux2 = np.trapz(Hb_spec2, x=pp.lam)

			bd[i_gal] = '<' + str(round(Ha_flux/Hb_flux2, 4))
			e_bd[i_gal] = '-'
		else:
			Hb_spec2 = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() \
				/ np.max(pp.matrix[:, pp.templatesToUse=='Hbeta']) \
				* np.median(noise[(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c))
				* (pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) * 2.5
			Hb_flux2 = np.trapz(Hb_spec2, x=pp.lam)
			Ha_spec2 = pp.matrix[:, pp.templatesToUse=='Halpha'].flatten() \
				/ np.max(pp.matrix[:, pp.templatesToUse=='Halpha']) \
				* np.median(noise[(pp.lam < 6563./(1 + (pp.sol[1][0] - 300)/c))
				* (pp.lam > 6563./(1 + (pp.sol[1][0] + 300)/c))]) * 2.5
			Ha_flux2 = np.trapz(Ha_spec2, x=pp.lam)

			bd[i_gal] = '>' + str(round(Ha_flux2/Hb_flux2, 4))
			e_bd[i_gal] = '-'

	adskj
	if instrument == 'vimos':
		temp = "{0:12}{1:10}{2:10}\n"
		with open(limits_file, 'w') as l:
			l.write(temp.format('Galaxy', 'Mass', 'e_Mass'))
			for i in range(len(galaxy_gals)):
				l.write(temp.format(galaxy_gals[i], mass[i], e_mass[i]))
	elif instrument =='muse':
		temp = "{0:12}{1:10}{2:10}{3:10}{4:10}\n"
		with open(limits_file, 'w') as l:
			l.write(temp.format('Galaxy', 'Mass', 'e_Mass', 'Bul_dec', 'e_Bul_dec'))
			for i in range(len(galaxy_gals)):
				l.write(temp.format(galaxy_gals[i], mass[i], e_mass[i], bd[i], 
					e_bd[i]))





if __name__=='__main__':
	if cc.device == 'uni':
		# whole_image('ngc1399', verbose=False, instrument='muse')
		galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
		for g in galaxies:
			whole_image(g, verbose=True, instrument='muse')

	elif 'home' in cc.device:
		# whole_image('ngc1399', verbose=False, instrument='vimos')
		galaxies = ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
			'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']
		# galaxies = ['ic1459']
		for g in galaxies:
			whole_image(g, verbose=False, instrument='vimos')
		# whole_image('pks0718-34', verbose=False)