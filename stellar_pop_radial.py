## ==================================================================
## 		Stellar population
## ==================================================================
## warrenj 20161214 Routine to bin in ellipses
import cPickle as pickle
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from absorption import absorption
from population2 import population
from tools import funccontains
from checkcomp import checkcomp
cc = checkcomp()

c = 299792.458 # speed of light in km/s

def stellar_pop(galaxy, wav_range="", vLimit=0, D=None):
	discard=2
	step_size=0.5

	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		'Mg_b']

	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range)
	out_plots = "%splots" % (output)

	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
		#pickleFile = open("%s/dataObj_%s.pkl" % (cc.home_dir, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()


	dataCubeDirectory = "%s/Data/vimos/cubes/%s.cube.combined.fits" % (cc.base_dir, 
		galaxy)
	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	s = galaxy_data.shape

	rows_to_remove = range(discard)
	rows_to_remove.extend([s[1]-1-i for i in range(discard)])
	cols_to_remove = range(discard)
	cols_to_remove.extend([s[2]-1-i for i in range(discard)])

	galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
	galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
	galaxy_noise = np.delete(galaxy_noise, rows_to_remove, axis=1)
	galaxy_noise = np.delete(galaxy_noise, cols_to_remove, axis=2)
	galaxy_badpix = np.delete(galaxy_badpix, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_badpix, cols_to_remove, axis=2)

	# Check for nan is data set.
	galaxy_badpix[np.isnan(galaxy_data)] = 1
	galaxy_data[galaxy_badpix==1] = 0
	galaxy_noise[galaxy_badpix==1] = 0.000000001
	s = galaxy_data.shape

	data_file = '%s/galaxies.txt' % (out_dir)
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]
	x_cent = x_gals[i_gal]
	y_cent = y_gals[i_gal]

	x = np.zeros(s[1]*s[2],dtype=int)
	y = np.zeros(s[1]*s[2],dtype=int)
	for i in range(s[1]):
		for j in range(s[2]):
			# Assign x and y
			x[i*s[1]+j] = i
			y[i*s[1]+j] = j


	lam = np.arange(s[0])*CDELT_spec+CRVAL_spec
	lam /= (1+z)
	lam /= (1+vel/c)
	ab_strength = {}
	uncerts = {}
	for l in lines:
		ab_strength[l] = []
		uncerts[l] = []

	for i in np.arange(0,s[1]/2,step_size):
		a = funccontains(annulus, x_cent ,y_cent, i, step_size, x=x,y=y)
		spec = np.nansum(galaxy_data[:,x[a],y[a]],axis=(1))
		noise = np.sqrt(np.nansum(galaxy_noise[:,x[a],y[a]]**2,axis=(1)))
		
		for l in lines:
			ab, uncert = absorption(l, lam, spec, noise=noise)
			ab_strength[l].append(ab)
			uncerts[l].append(uncert)
	for l in lines:
		ab_strength[l] = np.array(ab_strength[l]).flatten()
		uncerts[l] = np.array(uncerts[l]).flatten()







	x = np.arange(len(ab_strength[lines[0]]))*header['CDELT1']
	f, ax_array = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, sharex='col', 
		sharey='row')
	for i, l in enumerate(lines):
		ax_array[int(np.floor(i/2)),i%2].plot(x,ab_strength[l])
		ax_array[int(np.floor(i/2)),i%2].errorbar(x,ab_strength[l],yerr=uncerts[l])
		ax_array[int(np.floor(i/2)),i%2].set_title(l)
		ax_array[int(np.floor(i/2)),i%2].set_ylim([min(ab_strength[l])-2, max(ab_strength[l])+2])
	plt.show()
	lkjasdlk









	pop = population(ab_strength, uncerts)

	f,ax=plt.subplots(2,2)
	x = np.arange(len(pop.age))*header['CDELT1']
	ax[0,1].plot(x,pop.age)
	ax[0,1].errorbar(x,pop.age,yerr=pop.unc_age)
	ax[0,1].set_title('Age')
	ax[0,1].set_ylim([0,15.5])
	ax[1,0].plot(x,pop.metallicity)
	ax[1,0].errorbar(x,pop.metallicity, yerr=pop.unc_met)
	ax[1,0].set_title('Metallicity')
	ax[1,0].set_ylim([-2.3,0.7])
	ax[1,1].plot(x,pop.alpha)
	ax[1,1].set_title('Alpha')
	ax[1,1].errorbar(x,pop.alpha,yerr=pop.unc_alp)
	ax[1,1].set_ylim([-0.4,0.6])
	ax[0,0].plot(x,pop.red_chi2)
	ax[0,0].set_title('Reduced Chi^2')
	f.suptitle(galaxy.upper())

	
	pop.plot_probability_distribution(galaxy)
	

# need to consider angle
#def annulus(x, e, a1, d, t, x0, y0):
	# b1 = np.sqrt(a1**2 *(1-e**2))
	# a2 = a1+d
	# b2 = np.sqrt(a2**2 *(1-e**2))

def annulus(x, args):
	x0, y0, r1, d = args
	r2 = r1 + d

	a1 = np.sqrt(r1**2-(x-x0)**2)
	a2 = np.sqrt(r2**2-(x-x0)**2)

	return y0-a2, y0-a1, y0+a1, y0+a2






if __name__=="__main__":
	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
		'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxies=['ngc3100']
	for gal in galaxies:
		print gal
		stellar_pop(gal, wav_range='4200-', vLimit=2)
	plt.show()