# New routine to compare VIMOS observations with Ogando

import numpy as np 
from checkcomp import checkcomp
cc=checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
from prefig import Prefig
Prefig()
from astropy.io import fits
from errors2 import apply_range, get_dataCubeDirectory,run_ppxf, set_params
from pop import get_absorption
# from ppxf import create_plot
from compare_absorption2 import Lick_to_LIS
from scipy.ndimage.interpolation import rotate
from disk_fit_functions_binned import rebin
from spectools import spectrum




c = 299792.458 # speed of light in km/s
H = 70 # Mpc/(km/s)

# Works in pixel units
# c: center, l: length, w:width, pa:position angle (anticlockwise from 
# 	vertical?)
def slitFoV(c, l, w, pa, instrument='vimos', os=20):
	if instrument=='vimos':
		frac = np.zeros((40,40))
	elif instrument=='muse':
		frac = np.zeros((150,150))
	# pa = np.radians(pa)

	s = frac.shape

	# sampler = np.array(np.meshgrid(np.arange(s[0]*100), np.arange(s[1]*100))
	# 	).T/100.

	sampler = np.zeros((s[0]*os,s[1]*os))
	sampler[int(np.round(c[0] - w/2))*os: int(np.round(c[0] + w/2))*os, 
		int(np.round(c[1] - l/2))*os: int(np.round(c[1] + l/2))*os] = 1

	sampler = rotate(sampler, pa, reshape=False)

	i = np.arange(s[0])
	j = np.arange(s[1])
	i, j = np.meshgrid(i,j)
	i, j = i.flatten(), j.flatten()
	frac[j, i] = rebin(sampler, i*os+os/2, j*os+os/2, statistic='mean')
	
	return frac




def compare_absortion(galaxy, O_sig=False, corr_lines='all'):
	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header

	# Load VIMOS values
	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals, x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,4,5), dtype='float,int,int')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	center = np.array([x_cent_gals[i_gal], y_cent_gals[i_gal]])

	data_file =  "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)
	pa_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(3,))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	pa = pa_gals[i_gal]


	lines = ['H_beta', 'Fe5015', 'Mg_b']
	e_lines = ['e_H_beta', 'e_Fe5015', 'e_Mg_b']

	# load Ogando values
	cols = Ogando_data = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), dtype=str)[0]
	cols = [i for i, co in enumerate(cols) if co in lines or co in e_lines]
	Ogando_data = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), unpack=True, skiprows=2, 
		usecols=np.append([1,2],cols))
	galaxies = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), unpack=True, skiprows=2, usecols=(0,), dtype=str)
	i_gal = np.where(galaxies==galaxy)[0][0]
	O_sigma, O_sigma_err = 10**Ogando_data[0, i_gal], np.abs(
		10**Ogando_data[0, i_gal] * Ogando_data[0, i_gal] * 
		Ogando_data[1, i_gal]/10)

	O_val = {}
	O_err = {}
	for i in range(0, 2*len(lines), 2):
		O_val[lines[i/2]] = Ogando_data[i, i_gal]
		O_err[lines[i/2]] = Ogando_data[i+1, i_gal]
	
	params = set_params(reps=0, opt='pop', gas=1, lines=corr_lines, 
		produce_plot=False)

	mask = slitFoV(center, 4.1/header['CDELT1'], 2.5/header['CDELT2'], pa, 
		instrument='vimos')

	ifu = np.array(f[0].data)
	ifu[np.isnan(ifu)] = 0
	spec = np.einsum('ijk,jk->i',ifu, mask)

	ifu = np.array(f[1].data)
	ifu[np.isnan(ifu)] = 0
	noise = np.sqrt(np.einsum('ijk,jk->i',ifu**2, mask))

	lam = np.arange(len(spec))*header['CDELT3'] + header['CRVAL3']
	spec, lam, cut = apply_range(spec, lam=lam, set_range=params.set_range, 
		return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])
	noise = noise[cut]

	pp = run_ppxf(galaxy, spec, noise, lamRange, header['CDELT3'], params)

	

	if O_sig:
		if isinstance(O_sig, bool):
			absorp, uncert = get_absorption(lines, pp=pp, sigma=O_sigma(a),
				instrument='vimos')
			sigma = O_sigma(a)
		else:
			absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos',
				sigma=O_sigma(a)+O_sig*(pp.sol[0][1]-O_sigma(a)))
			sigma = O_sigma(a)+O_sig*(pp.sol[0][1]-O_sigma(a))
	else:
		absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos')
		sigma = pp.sol[0][1]



	my = []
	e_my = []
	og = []
	e_og = []
	sig = []
	lin = []

	for i, l in enumerate(lines):
		lin = np.append(lin, l)
		sig = np.append(sig, sigma)
		# Aperture correction:
		r_ab = 1.025*np.sqrt(4.1 * 2.5 / np.pi) # arcsec
		r_ab = np.radians(r_ab/60**2) * z * c / H * 1000 # kpc
		if l=='H_beta' or l=='Hbeta':
			l2='hb'
			beta = 0.002 # from table 5 in Ogando '08
			e_beta = 0.027
		elif l == 'Fe5015':
			l2 = l 
			beta = -0.012
			e_beta = 0.027
		elif l=='Mg_b':
			l2='mgb'
			beta = -0.031
			e_beta = 0.034
		elif l=='NaD':
			l2='nad'
			beta = -0.034
			e_beta = 0.022
		# elif l=='TiO1':
		# 	l2='tio1'
		# elif l=='TiO2':
		# 	l2='tio2'
		elif l=='Fe5270':
			l2='fe52'
			beta = -0.016
			e_beta = 0.025
		elif l=='Fe5335':
			l2='fe53'
			beta = -0.012
			e_beta = 0.027
		elif l=='Fe5406':
			l2=l 
			beta = -0.015
			e_beta = 0.029
		elif l=='Fe5702':
			l2=l 
			beta = 0
			e_beta = 0.036
	
		# Change to mag units	
		s = spectrum(lam=pp.lam, lamspec=pp.galaxy)
		I = -2.5 * np.log10(1 - absorp[l]/np.diff(getattr(s,l2))[0])
		e_I = np.abs(2.5/np.log(10) * 
			uncert[l]/(np.diff(getattr(s,l2))[0] - absorp[l]))

		I = I - beta * np.log10(r_ab/1.19) # Choosen from Davies 1987
		e_I = np.sqrt(e_I**2 + (e_beta**2 * np.log10(r_ab/1.19)))

		absorp[l] = (1 - 10**(-I/2.5))*np.diff(getattr(s,l2))[0] # Back to A
		uncert[l] = np.abs(2.5 * absorp[l] * np.log(10) * e_I)



		lit_value = Ogando_data[i*2+2, i_gal]
		e_lit_value = Ogando_data[i*2+3, i_gal]
		
		e_lit_value = np.mean([np.abs(Lick_to_LIS(l, lit_value + e_lit_value) 
			- Lick_to_LIS(l, lit_value)), 
			np.abs(Lick_to_LIS(l, lit_value - e_lit_value) 
			- Lick_to_LIS(l, lit_value))]) 
		lit_value = Lick_to_LIS(l, lit_value)


		my.append(absorp[l])
		e_my.append(uncert[l])

		og.append(lit_value)
		e_og.append(e_lit_value)

	return my, e_my, og, e_og, lin

if __name__=='__main__':
	lines = ['H_beta', 'Fe5015', 'Mg_b']
	color = ['r',        'b',       'g']
	# compare_absortion('ic1459', corr_lines = ['Hbeta','[OIII]5007d'])
	# fig, ax = plt.subplots()
	my = []
	e_my = []
	og = []
	e_og = []
	line = []
	for gal in ['ic1459', 'ic4296', 'ngc1399', 'ngc3557', 'ngc7075']:
		m, e_m, o, e_o, lin = compare_absortion(gal, corr_lines = None)
		my = np.append(my, m)
		e_my = np.append(e_my, e_m)
		og = np.append(og, o)
		e_og = np.append(e_og, e_o)
		line = np.append(line, lin)

	fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})

	ax.errorbar(my, og, xerr=e_my, yerr=e_og, c='k', fmt='.')
	
	xlims = ax.get_xlim()
	ylims = ax.get_ylim()

	ax.plot([0,10],[0,10], 'k')

	ax.set_xlim([np.min((xlims, ylims)), np.max((xlims, ylims))])
	ax.set_ylim([np.min((xlims, ylims)), np.max((xlims, ylims))])

	fig.savefig('%s/Data/lit_absorption/Ogando_comparison.png' % (
		cc.base_dir))

	for l in np.unique(line):
		m = line == l

		print l, 'Offset:', np.mean(my[m] - og[m]), 'dispersion:', \
			np.std(my[m] - og[m])

