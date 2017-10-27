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
import cPickle as pickle
from astropy.io import fits
from classify import get_R_e
from errors2 import apply_range, get_dataCubeDirectory,run_ppxf, set_params
from glob import glob
from pop import get_absorption
from ppxf import create_plot
from spectools import spectrum
from scipy.interpolate import interp1d
from compare_absorption2 import Lick_to_LIS
from scipy.ndimage.interpolation import rotate
from scipy.stats import binned_statistic
from scipy.spatial import cKDTree as KDTree
from disk_fit_functions_binned import rebin



c = 299792.458 # speed of light in km/s

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

	lines = ['H_beta', 'Fe5015', 'Mg_b']
	e_lines = ['e_H_beta', 'e_Fe5015', 'e_Mg_b']

	# load Ogando values
	cols = Ogando_data = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), dtype=str)[0]
	cols = [i for i, c in enumerate(cols) if c in lines or c in e_lines]
	Ogando_data = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), unpack=True, skiprows=2, 
		usecols=np.append([1,2,3],cols))
	galaxies = np.loadtxt('%s/Data/lit_absorption/Ogando.txt' % (
		cc.base_dir), unpack=True, skiprows=2, usecols=(0,), dtype=str)
	i_gal = np.where(galaxies==galaxy)[0][0]
	O_sigma, O_sigma_err, rad = 10**Ogando_data[0, i_gal], np.abs(
		10**Ogando_data[0, i_gal] * Ogando_data[0, i_gal] * 
		Ogando_data[1, i_gal]/10), 10**Ogando_data[2, i_gal]

	O_val = {}
	O_err = {}
	for i in range(0, 2*len(lines), 2):
		O_val[lines[i/2]] = Ogando_data[i, i_gal]
		O_err[lines[i/2]] = Ogando_data[i+1, i_gal]

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

	index = np.zeros((40,40,2))
	for i in range(40):
		for j in range(40):
			index[i,j,:] = np.array([i,j]) - center
	
	params = set_params(reps=0, opt='pop', gas=1, lines=corr_lines, 
		produce_plot=False)

	# **********************
	mask = slitFoV(center, 4.1/header['CDELT1'], 2.5/header['CDELT2'], pa, 
		instrument='vimos')

	ifu = np.array(f[0].data)
	ifu[np.isnan(ifu)] = 0
	spec = np.einsum('ijk,jk->i',ifu, mask)

	ifu = np.array(f[1].data)
	ifu[np.isnan(ifu)] = 0
	spec = np.sqrt(np.einsum('ijk,jk->i',ifu**2, mask))
	# **********************
	lkadsjklas
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
			sigma = np.append(sigma, O_sigma(a))
		else:
			absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos',
				sigma=O_sigma(a)+O_sig*(pp.sol[0][1]-O_sigma(a)))
			sigma = np.append(sigma, O_sigma(a)+O_sig*(pp.sol[0][1]-O_sigma(a)))
	else:
		absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos')
		sigma = np.append(sigma, pp.sol[0][1])

	for l in lines:
		if a == min(apertures):
			my_values[l] = np.array([])
			my_errors[l] = np.array([])
		my_values[l] = np.append(my_values[l], absorp[l])
		my_errors[l] = np.append(my_errors[l], uncert[l])
	
	for i, l in enumerate(lines):
		ax.errorbar(a, absorp[l], yerr=uncert[l], color=color[i], fmt='x')
	for i, l in enumerate(lines):
	  ax.errorbar(np.nan, np.nan, color=color[i], fmt='x', label=l)
	ax.legend(facecolor='w')


	Rampazzo_file = '%s/Data/lit_absorption/Rampazzo_aperture.txt' % (cc.base_dir)
	file_headings = np.loadtxt(Rampazzo_file, dtype=str)[0]

	for i, l in enumerate(lines):
		col = np.where(file_headings==l)[0][0]
		try:
			col2 = np.where(file_headings==l)[0][1]
		except:
			try:
				col2 = np.where(file_headings=='_'+l)[0][0]
			except:
				col2 = np.where(file_headings=='e_'+l)[0][0]
		R_obs, R_err = np.loadtxt(Rampazzo_file, unpack=True, skiprows=1, 
			usecols=(col,col2))
		R_galaxies = np.loadtxt(Rampazzo_file, unpack=True, skiprows=1, usecols=(0,), 
			dtype=str)

		mask = R_galaxies==galaxy

		order = np.argsort(apertures)

		lit_value = Lick_to_LIS(l, R_obs[mask][order])
		err = np.mean([np.abs(Lick_to_LIS(l, R_obs[mask][order] + 
			R_err[mask][order]) - Lick_to_LIS(l, R_obs[mask][order])), 
			np.abs(Lick_to_LIS(l, R_obs[mask][order] - R_err[mask][order]) -
			Lick_to_LIS(l, R_obs[mask][order]))], axis=0) 

		ax.errorbar(apertures[order], lit_value, yerr=err, color=color[i])

		if l=='H_beta' or l=='Hbeta':
			l2='hb'
		elif l=='Mg_b':
			l2='mgb'
		elif l=='NaD':
			l2='nad'
		elif l=='TiO1':
			l2='tio1'
		elif l=='TiO2':
			l2='tio2'
		elif l=='Fe5270':
			l2='fe52'
		elif l=='Fe5335':
			l2='fe53'
		else:
			l2=l

		ax3.scatter(
		  np.abs(my_values[l][order] - lit_value)/my_values[l][order],
		  np.abs(my_values[l][order] - lit_value)/np.sqrt(err**2 + 
		  my_errors[l][order]**2), color=color[i], s=4*apertures[order]**2,
		  label=l)
		ax4.scatter(sigma[order], 
		  np.abs(my_values[l][order] - lit_value)/my_values[l][order],
		  color=color[i], s=4*apertures[order]**2, label=l)
		ax5.scatter(sigma[order], 
		  np.abs(my_values[l][order] - lit_value)/np.sqrt(err**2 + 
		  my_errors[l][order]**2), marker='x',
		  color=color[i], s=4*apertures[order]**2, label=l)
		t.extend(sigma[order])
		g.extend(my_values[l][order]) 
		h.extend(lit_value)
		j.extend(err)
		e.extend(my_errors[l][order])
		r.extend([lines[i]]*len(sigma))
		w.extend(4*apertures[order]**2)


	ax.set_ylabel(r'Index strength, $\AA$')
	ax.set_xlabel('Radius, arcsec')
	fig.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_%i.png' % (
	  cc.base_dir, galaxy, params.gas))
	plt.close(fig)


	ax3.set_ylabel(r'Sigma difference ((Mine - Ramp)/Combined Uncert)')
	ax3.set_xlabel('Fractional difference ((Mine - Ramp)/Mine)')
	ax3.legend()
	fig3.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_fractional.png' % (
	  cc.base_dir, galaxy))
	plt.close(fig3)


	ax4.set_xlabel('Vel dispersion')
	ax3.set_ylabel('Fractional difference ((Mine - Ramp)/Mine)')
	ax5.set_ylabel(r'Sigma difference ((Mine - Ramp)/Combined Uncert)')
	ax4.legend()
	fig4.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_sigma.png' % (
	  cc.base_dir, galaxy))
	plt.close(fig4)


	return t, g, h, j, e, r, w

if __name__=='__main__':
	lines = ['H_beta', 'Fe5015', 'Mg_b']
	color = ['r',        'b',       'g']
	slitFoV((20,20), 10, 3, 45, instrument='vimos')
	# compare_absortion('ic1459', corr_lines = ['Hbeta','[OIII]5007d'])
	# fig, ax = plt.subplots()
	# for gal in ['ic1459','ic4296','ngc3557']:
	# 	s, f, g, h, j, r, z = compare_absortion(gal, 
	# 		corr_lines = ['Hbeta','[OIII]5007d'])
		

