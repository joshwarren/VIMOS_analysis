# ==================================================================
# 		MCMC to find systematic v and sigma
# ==================================================================
# warrenj 20160210 Changed to fit the whole galaxy spectrum and have
# combined mcmc and mcmc_fit_bin into one routine. 
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy import ndimage # for gaussian blur
from astropy.io import fits

from ppxf import ppxf
import ppxf_util as util
from errors2 import apply_range, use_templates, determine_goodpixels
from find_template import setup
from checkcomp import checkcomp
cc=checkcomp()

c = 299792.458
repeats = 100
threshold =0.85
local_step = 30

## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
def mcmc(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000]):
	print '     MCMC'
	
	results = np.zeros((2,repeats))

	for i in range(2):
		templates, bin_log, noise, velscale, start, goodpixels,	moments, degree, dv, \
			lambdaq, plot, quiet = setup(galaxy, z=z, vel=vel, sig=sig, discard=discard, 
			set_range=set_range)

		pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodpixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=plot, quiet=quiet)
		z = (z + 1)*np.sqrt((1 + pp.sol[0]/c)/(1 - pp.sol[0]/c)) - 1 

	templates, bin_log, noise, velscale, start, goodpixels,	moments, \
		degree, dv, lambdaq, plot, quiet = 	setup(galaxy, z=z, vel=vel, 
		sig=sig, discard=discard, set_range=set_range)

	chi = 1000000
	for i in range(repeats):
		v_sav = vel 
		sigma_sav = sig
		chi_sav = chi
		if not quiet:
			print 'i=',i
			print "input v: ",v_sav
			print "input sig: ",sigma_sav
			print "input chi2: ",chi_sav

		start = [vel,sig]

		if i != repeats-1:
			pp = ppxf(templates, bin_log, noise, velscale, start, 
				goodpixels=goodpixels, moments=moments, degree=degree, 
				vsyst=dv, lam=lambdaq, plot=not quiet, quiet=quiet)
		else:
			pp = ppxf(templates, bin_log, noise, velscale, start, 
				goodpixels=goodpixels, moments=moments, degree=degree, 
				vsyst=dv, lam=lambdaq, plot=False, quiet=quiet)

		vel = pp.sol[0]
		sig = pp.sol[1]
		chi = pp.chi2

		results[:,i] = [vel,sig]
		if abs(chi) > abs(chi_sav) or np.random.uniform() > threshold:
			vel = v_sav + np.random.uniform(low=-1, high=1)*local_step 
			sig = sigma_sav + np.random.uniform(low=-1, high=1)*local_step 
			chi = chi_sav

	vel = np.mean(results[0,:])
	sig = np.mean(results[1,:])

	if not quiet:
		print "Mean vel: :",vel
		print "MEAN vel dispersion: ",sig

# ----------===============================================---------
# ----------================= Save Result =================---------
# ----------===============================================---------

	dir = '%s/Data/vimos'  %(cc.base_dir)

	fig, ax = plt.subplots()
	ax.scatter(results[0,:],results[1,:])
	ax.scatter(vel,sig, marker='*')
	ax.set_xlabel("velocity")
	ax.set_ylabel("velocity dispersion")
	ax.set_title("MCMC for initial conditions")

	fig.savefig('%s/analysis/%s/MCMC_initial_fit.png' %(dir,galaxy))
	if plot: 
		plt.show()



	data_file = "%s/analysis/galaxies.txt" % (dir)
	d = np.loadtxt(data_file, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	z_gals, vel_gals, sig_gals = d[1][1:].astype(float), d[2][1:].astype(float), \
		d[3][1:].astype(float),
	x_gals, y_gals = d[4][1:].astype(int), d[5][1:].astype(int)
	SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(6,len(d))}

	# If galaxy is already in galaxies.txt file
	try:
		i_gal = np.where(galaxy_gals == galaxy)
	except:
		i_gal = -1
		galaxy_gals = [galaxy_gals, galaxy]
		z_gals = [z_gals, z]
		vel_gals = [vel_gals, vel]
		sig_gals = [sig_gals, sig]
	else:
		z_gals[i_gal] = z
		vel_gals[i_gal] = vel
		sig_gals[i_gal] = sig

	temp = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}" + ''.join(['{%i:%i}'%(i+6,len(t)+1) 
		for i, t in enumerate(SN_gals.keys())]) + "\n"

	SN_titles = list(SN_gals.keys())

	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "z", "velocity", "sigma", "x", "y", 
			*(s for s in SN_titles)))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(z_gals[i],7)), 
				str(round(vel_gals[i],4)), str(round(sig_gals[i],4)), 
				str(int(x_gals[i])), str(int(y_gals[i])), 
				*(str(round(SN_gals[s][i],2)) for s in SN_titles)))