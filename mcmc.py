# ==================================================================
# 		MCMC to find systematic v and sigma
# ==================================================================
# warrenj 20160210 Changed to fit the whole galaxy spectrum and have
# combined mcmc and mcmc_fit_bin into one routine. 
import numpy as np
import glob
from ppxf import ppxf
import ppxf_util as util
from checkcomp import checkcomp
from errors2 import determine_goodpixels
cc=checkcomp()

quiet = True
plot = True
c = 299792.458

repeats = 100
threshold =0.85
local_step = 30

def setup(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000]):
# ----------===============================================---------
# ----------============= Input parameters ================---------
# ----------===============================================---------
	FWHM_gal = 4*0.571 # The fibre FWHM on VIMOS is
					   # about 4px with a dispersion of
					   # 0.571A/px. (From: http://www.eso.org
					   # /sci/facilities/paranal/instruments
					   # /vimos/inst/ifu.html)
	moments = 4 # number of componants to calc with ppxf (see 
				# keyword moments in ppxf.pro for more details)
	degree = 4 # order of addative Legendre polynomial used to 
				# correct the template continuum shape during the fit 	
	dir = '%s/Data/vimos/'  %(cc.base_dir)
	dir2 = '%s/Data/idl_libraries' % (cc.base_dir)

# ----------===============================================---------
# ----------=============== Run analysis  =================---------
# ----------===============================================---------


## ----------=============== Miles library =================---------
	# Finding the template files
	templatesDirectory = dir2 + 'ppxf/MILES_library/'
	templateFiles = glob.glob(templatesDirectory + \
		'm0[0-9][0-9][0-9]V') #****** no longer have nfiles.


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


	## ****************************************************************
	## NB: shouldn't this be 0.9A as this is resolution?
	FWHM_tem = 2.5 # Miles spectra have a resolution
							   # FWHM of 2.5A.


	## Which templates to use are given in use_templates.pro. This is
	## transfered to the array templatesToUse.
	templatesToUse = use_templates(galaxy, glamdring)
	nfiles = len(templatesToUse)
	templates = np.zeros((len(log_temp_template), nfiles))
	FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
	sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels


	## Reading the contents of the files into the array templates. 
	## Including rebinning them.
	for i in range(nfiles):
		v1, v2 = np.loadtxt(templateFiles[templatesToUse[i]], unpack='True')
		v2 = ndimage.gaussian_filter1d(v2,sigma)
		## Rebinning templates logarthmically
		log_temp_template, logLam_template, velscale = \
			util.log_rebin(lamRange_template, v2, velscale=velscale)
	## ****************************************************************
	## ^^^ this has changed from the last time we called this: we called v1 before...
		templates[:,i] = log_temp_template

	#templates /= np.median(log_temp_template)
## ----------========= Reading Tessellation  ===============---------

	## Reads the txt file containing the output of the binning_spaxels
	## routine. 
	x,y,bin_num = np.loadtxt(tessellation_File, usecols=(0,1,2), \
		unpack=True, skiprows=1)

	n_bins = max(bin_num) + 1
	## Contains the order of the bin numbers in terms of index number.
	order = np.sort(bin_num)
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = glob.glob(dir+"cubes/%s.cube.combined.fits" % (galaxy)) 
		
	galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)
	galaxy_noise = pyfits.getdata(dataCubeDirectory[0], 1)

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

	n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])
## ----------============= Spatially Binning ===============---------
	gal_spec = np.nansum(np.nansum(galaxy_data, axis=2),axis=1)
	gal_noise = np.nansum(np.nansum(galaxy_noise**2, axis=2),axis=1)
	gal_noise = np.sqrt(gal_noise)
	

## ----------======== Finding limits of the spectrum =======---------
	## limits are the cuts in pixel units, while lamRange is the cuts in
	## wavelength unis.
	gap=12
	ignore = int((5581 - CRVAL_spec)/CDELT_spec) + np.arange(-gap+1,gap)  
	ignore2 = int((5199 - CRVAL_spec)/CDELT_spec) + np.arange(-gap+1,gap) 

	## h is the spectrum with the peak enclosed by 'ignore' removed.
	h = np.delete(gal_spec, ignore)
	h = np.delete(h,ignore2)


	half = s[0]/2
	a = np.delete(h,np.arange(-4,0)+len(h),None)/np.median(h[np.nonzero(h)]) - \
		h[4:]/np.median(h[np.nonzero(h)])
	a = np.where(np.isfinite(a), a, 0)

	if any(np.abs(a[:0.5*half]) > 0.2):
		lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
	else: 
		lower_limit = -1
	
	# lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
	if any(np.abs(a[1.5*half:]) > 0.2):
		upper_limit = min(np.where(np.abs(a[1.5*half:]) > 0.2)[0])+int(1.5*half)
	else:
		upper_limit = -1
		
	if upper_limit > ignore2[0]: upper_limit+=gap 
	if upper_limit > ignore[0]: upper_limit+=gap

	if lower_limit < 0:
		lower_limit = min(np.where(a[:half] != 0)[0]) + 5
		if lower_limit < 0: lower_limit = 0

	else:
		lower_limit +=5

	if upper_limit > s[0]-1 or upper_limit < half:
		upper_limit = s[0]-6 
	else:
		upper_limit += -5
## ----------========== Using set_range variable ===========---------
	if set_range is not None:
	## Change to pixel units
		set_range = ((set_range - CRVAL_spec)/CDELT_spec).astype(int)
		if set_range[0] > lower_limit: lower_limit = set_range[0] 
		if set_range[1] < upper_limit: upper_limit = set_range[1]

	lamRange = np.array([lower_limit, upper_limit])*CDELT_spec + \
		CRVAL_spec
## ----------=========== Writing the spectrum  =============---------
	bin_lin = gal_spec[lower_limit:upper_limit]
	bin_lin_noise = gal_noise[lower_limit:upper_limit]
## ----------========= Calibrating the spectrum  ===========---------
	## For calibrating the resolutions between templates and observations
	## using the gauss_smooth command
	# FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_gal**2)
	# sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels

	## smooth spectrum to fit with templates resolution
	# bin_lin = ndimage.gaussian_filter1d(bin_lin, sigma)
	# bin_lin_noise = ndimage.gaussian_filter1d(bin_lin_noise, sigma)
	
	lamRange = lamRange/(1+z)
	## rebin spectrum logarthmically
	bin_log, logLam_bin, velscale = util.log_rebin(lamRange, bin_lin, 
		velscale=velscale)
	bin_log_noise, logLam_bin, velscale = util.log_rebin(lamRange, 
		np.power(bin_lin_noise,2), velscale=velscale)
	bin_log_noise = np.sqrt(bin_log_noise)

	noise = bin_log_noise+0.0000000000001



	dv = (logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodPixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z) 
	lambdaq = np.exp(logLam_bin)
	start = [vel, sig] # starting guess

	goodPixels = determine_goodpixels(logLam_bin,lamRange_template,
		vel, z, gas=True)


	noise = np.abs(noise)
	bin_log_sav = bin_log
	noise_sav = noise

	return templates, bin_log, noise, velscale, start, 
		goodpixels, moments, degree, dv, lambdaq, not quiet, quiet







## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
def mcmc(galaxy, z=0.01, vel=0.0, sig=200.0, discard=2, set_range=[4200,10000]):
	
	results = np.zeros((2,repeats))

	for i in range(2):
		pp = ppxf(setup(galaxy, z=z, vel=vel, sig=sig, discard=discard, 
			set_range=set_range))
		z = (z + 1)*sqrt((1 + pp.sol[0]/c)/(1 - pp.sol[0]/c)) - 1 

	templates, bin_log, noise, velscale, start, goodpixels,	moments, degree, dv, \
		lambdaq, plot, quiet = 	setup(galaxy, z=z, vel=vel, sig=sig, discard=discard, 
		set_range=set_range)

	chi = 1000000
	for i in range(repeats):
		print 'i=',i
		v_sav = vel 
		sigma_sav = sig
		chi_sav = chi
		print "input v: ",v_sav
		print "input sig: ",sigma_sav
		print "input chi2: ",chi_sav

		start = [vel,sig]


		pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet)

		vel = pp.sol[0]
		sig = pp.sol[1]
		chi = pp.sol[6]

		results[:,i] = [vel,sig]
		if abs(chi-1) > abs(chi_sav-1) or np.random.uniform > threshold:
			vel += np.random.uniform(low=-1, high=1)*step_range 
			sig += np.random.uniform(low=-1, high=1)*step_range 
			chi = chi_sav
	vel = np.mean(results[0,:])
	sig = np.mean(results[1,:])

	print "Mean vel: :",vel
	print "MEAN vel dispersion: ",sig

# ----------===============================================---------
# ----------================= Save Result =================---------
# ----------===============================================---------

	dir = '%s/Data/vimos/'  %(cc.base_dir)

	fig, ax = plt.figure()
	ax.scatter(results[0,:],results[1,:])
	ax.scatter(vel,sig, marker='*')
	ax.xlabel("velocity")
	ax.ylabel("velocity dispersion")
	ax.title("MCMC for initial conditions")

	plt.show()
	f.savefig('%s/analysis/%s/MCMC_initial_fit.png' %(dir,galaxy))





	data_file = "%s/analysis/galaxies.txt" % (dir)
	galaxy_gals = np.loadtxt(data_file, usecols=(0,), unpack=True, 
		dtype=str, skiprows=1)
	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_used_gals = np.loadtxt(
		data_file, skiprows=1, usecols=(1,2,3,4,5,6), 
		dtype='float,float,float,int,int,float')

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

	np.savetxt(data_file, np.column_stack([galaxy_gals, z_gals, vel_gals[i],
		sig_gals[i], x_gals[i], y_gals[i], SN_used_gals]), 
		fmt='%s %10.9f %10.9f %10.5f %2i %2i %2.2f', 
		header ="Galaxy	  z	 velocity	 velocity dispersion	x	 y	 Target SN")
