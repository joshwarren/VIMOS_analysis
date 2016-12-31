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
from errors2 import remove_anomalies, use_templates, determine_goodpixels
from checkcomp import checkcomp
cc=checkcomp()

quiet = True
plot = False
c = 299792.458

repeats = 100
threshold =0.85
local_step = 30

def setup(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000]):
# ----------===============================================---------
# ----------============= Input parameters ================---------
# ----------===============================================---------
	FWHM_gal = 2.5/(1+z) # VIMOS documentation (and fits header)
	moments = 4 # number of componants to calc with ppxf (see 
				# keyword moments in ppxf.pro for more details)
	degree = 4 # order of addative Legendre polynomial used to 
				# correct the template continuum shape during the fit 	
	dir = '%s/Data/vimos/'  %(cc.base_dir)
	templatesDirectory = '%s/models/miles_library' % (cc.home_dir)

# ----------===============================================---------
# ----------=============== Run analysis  =================---------
# ----------===============================================---------


## ----------=============== Miles library =================---------
	# Finding the template files
	templateFiles = glob.glob(templatesDirectory + \
		'/m0[0-9][0-9][0-9]V') #****** no longer have nfiles.


	# v1 is wavelength, v2 is spectrum
	v1, v2 = np.loadtxt(templateFiles[0], unpack='True')

	# Using same keywords as fits headers
	CRVAL_temp = v1[0]		# starting wavelength
	NAXIS_temp = np.shape(v2)[0]   # Number of entries
	# wavelength increments (resolution?)
	CDELT_temp = (v1[NAXIS_temp-1]-v1[0])/(NAXIS_temp-1)

	lamRange_template = CRVAL_temp + [0, CDELT_temp*(NAXIS_temp-1)]

	log_temp_template, logLam_template, velscale = \
		util.log_rebin(lamRange_template, v1)

	FWHM_tem = 2.5 # Miles library has FWHM of 2.5A.
	FWHM_dif = np.sqrt(abs(FWHM_gal**2 - FWHM_tem**2))

	## Which templates to use are given in use_templates.pro. This is
	## transfered to the array templatesToUse.
	templatesToUse = use_templates(galaxy, cc.device=='glamdring')
	nfiles = len(templatesToUse)
	templates = np.zeros((len(log_temp_template), nfiles))

	## Reading the contents of the files into the array templates. 
	## Including rebinning them.
	for i in range(nfiles):
		v1, v2 = np.loadtxt(templateFiles[templatesToUse[i]], unpack='True')
		if FWHM_tem < FWHM_gal:
			sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels
			v2 = ndimage.gaussian_filter1d(v2,sigma)
		## Rebinning templates logarthmically
		log_temp_template, logLam_template, _ = \
			util.log_rebin(lamRange_template, v2, velscale=velscale)
	## ****************************************************************
	## ^^^ this has changed from the last time we called this: we called v1 before...
		templates[:,i] = log_temp_template

	#templates /= np.median(log_temp_template)
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = glob.glob(dir+"cubes/%s.cube.combined.fits" % (galaxy)) 
		
	galaxy_data, header = fits.getdata(dataCubeDirectory[0], 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory[0], 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory[0], 3)

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

	n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])
## ----------============= Spatially Binning ===============---------
	gal_spec = np.nansum(np.nansum(galaxy_data, axis=2),axis=1)
	bin_lin_noise = np.sqrt(np.nansum(np.nansum(galaxy_noise**2, axis=2),axis=1))
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	bin_lin, lam, cut = remove_anomalies(gal_spec, window=201, repeats=3, 
		lam=lam, set_range=set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
	bin_lin_noise = bin_lin_noise[cut]

	## smooth spectrum to fit with templates resolution
	if FWHM_gal < FWHM_tem:
		sigma = FWHM_dif/2.355/CDELT_spec # Sigma difference in pixels
		bin_lin = ndimage.gaussian_filter1d(bin_lin, sigma)
		bin_lin_noise = np.sqrt(ndimage.gaussian_filter1d(bin_lin_noise**2, sigma))

	## rebin spectrum logarthmically
	bin_log, logLam_bin, _ = util.log_rebin(lamRange, bin_lin, velscale=velscale)
	bin_log_noise, logLam_bin, _ = util.log_rebin(lamRange, 
		bin_lin_noise**2, velscale=velscale)
	bin_log_noise = np.sqrt(bin_log_noise)

	noise = bin_log_noise+0.0000000000001



	dv = (logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodpixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z) 
	lambdaq = np.exp(logLam_bin)
	start = [vel, sig] # starting guess

	goodpixels = determine_goodpixels(logLam_bin,lamRange_template,
		vel, z, gas=True)


	noise = np.abs(noise)
	bin_log_sav = bin_log
	noise_sav = noise

	return templates, bin_log, noise, velscale, start, goodpixels, moments, \
		degree, dv, lambdaq, not quiet, quiet







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

	templates, bin_log, noise, velscale, start, goodpixels,	moments, degree, dv, \
		lambdaq, plot, quiet = 	setup(galaxy, z=z, vel=vel, sig=sig, discard=discard, 
		set_range=set_range)

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


		pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodpixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet)

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

	dir = '%s/Data/vimos/'  %(cc.base_dir)

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
	galaxy_gals = np.loadtxt(data_file, usecols=(0,), unpack=True, 
		dtype=str, skiprows=1)
	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_kin_gals, SN_pop_gals = np.loadtxt(
		data_file, skiprows=1, usecols=(1,2,3,4,5,6,7), unpack=True,
		dtype='float,float,float,int,int,float,float')

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

	temp = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}{6:8}{7:8}\n"
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "z", "velocity", "vel dispersion", "x", "y", 
			"Kin SN", "Pop SN"))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(z_gals[i],7)), 
				str(round(vel_gals[i],4)), str(round(sig_gals[i],4)), 
				str(int(x_gals[i])), str(int(y_gals[i])), str(round(SN_kin_gals[i],2)),
				str(round(SN_pop_gals[i],2))))