## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import math # for sine functions
import matplotlib.pyplot as plt # used for plotting
import matplotlib.axes as ax # for adding text onto images
from scipy.optimize import curve_fit # for fitting a gaussian
from checkcomp import checkcomp
cc = checkcomp()
import cPickle as pickle

#---------------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, discard=0, wav_range="", plots=False):

	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""
		wav_range = ""
	galaxiesFile =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	galaxiesFile2 =  "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)
	galaxiesFile_Re =  "%s/Data/vimos/analysis/galaxies_R_e.txt" % (cc.base_dir)
	output = '%s/Data/vimos/analysis/%s/results/%s' % (cc.base_dir, galaxy, 
		wav_range_dir)
	pickleFile = open('%s/pickled/dataObj_%s.pkl' % (output, wav_range))
	D = pickle.load(pickleFile)
	pickleFile.close()


	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_gals = np.loadtxt(galaxiesFile, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5,6))
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	ellip_gals, star_mis, OIII_mis, Hbeta_mis, Hdelta_mis, Hgamma_mis = np.loadtxt(
		galaxiesFile2, unpack=True, skiprows=1, usecols=(1,2,3,4,5,6))
	gas = {'[OIII]5007d':OIII_mis, 'Hbeta':Hbeta_mis, 'Hdelta':Hdelta_mis, 'Hgamma':Hgamma_mis}


	log_R_e_RC3_gals, R_e_2MASS_gals = np.loadtxt(galaxiesFile_Re, unpack=True, 
		skiprows=1, usecols=(1,2))
	R_e_RC3 = 6*10**log_R_e_RC3_gals[i_gal]/2 # convert to arcsec
	R_e_2MASS = R_e_2MASS_gals[i_gal]

	R_e = np.nanmean([R_e_RC3,R_e_2MASS])
# ------------=============== Photometry =================----------
# ------------============= Fit photometry ===============----------
	save_to = "%s/Data/vimos/analysis/%s/results/" % (cc.base_dir, 
		galaxy) + "%splots/photometry_%s.png" % (wav_range_dir, wav_range)
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed ?
	ellip_gals[i_gal] = f.eps
	print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))


# ------------================ Kinemetry =================----------
# ------------============== Fit kinemetry ===============----------
	save_to = "%s/Data/vimos/analysis/%s/results/" % (cc.base_dir, 
		galaxy) + "%s/plots/stellar_kinematics_%s.png" % (wav_range_dir, 
		wav_range)
 	xBar = np.array(D.xBar)-f.xmed
 	yBar = np.array(D.yBar)-f.ymed


	k = fit_kinematic_pa(xBar, yBar, np.array(D.components['stellar'].plot['vel']), 
		quiet=True, plot=plots, sav_fig=save_to)

	print "PA_kin: " + str(k[0]) + "+/-" + str(k[1]/3)


# ------------============== Misalignment ================----------
	phot = math.radians(90-f.theta)
	kine = math.radians(k[0])


	mis = math.asin(abs(math.sin(phot-kine)))
	mis = math.degrees(mis)
	star_mis[i_gal] = mis
	print "Stars misalignment: " + str(mis)

# ------------================= Lambda_R ================----------
# R is distance from axis of rotation NOT distance from center of galaxy.	

	# distance from center
	R = np.sqrt(np.square(xBar)+np.square(yBar))
	# Angle of r from axis of rotation
	#ang = abs(math.asin(math.sin(math.radians(f.theta))) - np.abs(
	#	np.arctan((xBar)/(yBar))))
	#R = r * np.sin(ang)

	order = np.argsort(R)

	# NB: lam is ordered in terms of increasing R.
	lam_num = D.flux[order]*R[order]*np.abs(np.array(D.components['stellar'].plot['vel'])[order]) # numerator
	lam_den = D.flux[order]*R[order]*np.sqrt(np.square(
		np.array(D.components['stellar'].plot['vel'])[order]) + np.square(
		np.array(D.components['stellar'].plot['sigma'])[order])) # denominator

	# cumulative summation
	lam_num = np.cumsum(lam_num)
	lam_den = np.cumsum(lam_den)

	lam = lam_num/lam_den
	plt.figure()
	plt.title(r"Radial $\lambda_R$ profile")
	plt.ylabel(r"$\lambda_R$")
	plt.xlabel("Radius (R_e)")
	x = spxToRe(R[order], R_e)
	plt.plot(x, lam)
	ax =plt.gca()
	plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), verticalalignment='top',
		transform=ax.transAxes)
	plt.savefig("%s/Data/vimos/analysis/%s/results/" % (cc.base_dir, 
		galaxy) + "%s/plots/lambda_R_%s.png" % (wav_range_dir, wav_range), \
		bbox_inches="tight")
	if plots: 
		plt.show()
# ------------================= Gas ================----------
	for c in D.e_components:
		print c
		save_to = "%s/Data/vimos/analysis/%s/results/" % (cc.base_dir, 
			galaxy) + "%s/plots/%s_kinematics_%s.png" % (wav_range_dir, c, wav_range)
		gas_k = fit_kinematic_pa(xBar, yBar, np.array(D.components[c].plot['vel']), 
			quiet=True, plot=plots, sav_fig=save_to)

		print "%s_PA_kin: " % (c) + str(gas_k[0]) + "+/-" + str(gas_k[1]/3)

# ------------============== Misalignment ================----------
		gas_kine = math.radians(gas_k[0])

		gas_mis = math.asin(abs(math.sin(gas_kine-kine)))
		gas_mis = math.degrees(gas_mis)

		gas[c][i_gal] = gas_mis

		#gas_psi_gals[i_gal] = gas_mis
		print "Mis-alignment: " + str(gas_mis)


# ------------============= Save outputs =============----------
	f = open(galaxiesFile, 'wb')
	f.write('Galaxy    z     velocity     velocity dispersion      x    y' +
		'    Target SN \n')
	f2 = open(galaxiesFile2, 'wb')
	f2.write('Galaxy    Ellipticity    Misalignment: Stellar    OIII    Hbeta' +
		'    Hdelta   Hgamma\n')

	for i in range(len(galaxy_gals)):
		f.write(galaxy_gals[i] + '   ' + str(z_gals[i]) + '   ' + \
			str(vel_gals[i]) + '   ' + str(sig_gals[i]) + '   ' + \
			str(int(x_gals[i])) + '   ' + str(int(y_gals[i])) + '   ' + \
			str(SN_gals[i]) + '\n')

		f2.write(galaxy_gals[i] + '   ' + str(ellip_gals[i]) + '   ' + \
			str(star_mis[i]) + '   ' + str(OIII_mis[i]) + '   ' + \
			str(Hbeta_mis[i]) + '   ' + str(Hdelta_mis[i]) + '   ' + \
			str(Hgamma_mis[i]) + '\n')

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc3557'
	discard = 2
	wav_range = '4200-'

	kinematics(galaxy, discard=discard, wav_range=wav_range, plots=True)
