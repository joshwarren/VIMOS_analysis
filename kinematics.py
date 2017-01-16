## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import os
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import math # for sine functions
import matplotlib.pyplot as plt # used for plotting
import matplotlib.axes as ax # for adding text onto images
from scipy.optimize import curve_fit # for fitting a gaussian
from checkcomp import checkcomp
from classify import get_R_e
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


def kinematics(galaxy, discard=0, wav_range="", plots=False, D=None):

	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)

	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""
		wav_range = ""
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	galaxiesFile2 = "%s/galaxies2.txt" % (analysis_dir)
	output = '%s/%s/results/%s' % (analysis_dir, galaxy, wav_range_dir)
	if D is None:
		pickleFile = open('%s/pickled/dataObj_%s.pkl' % (output, wav_range))
		D = pickle.load(pickleFile)
		pickleFile.close()


	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_kin_gals, SN_pop_gals = np.loadtxt(
		galaxiesFile, unpack=True, skiprows=1, usecols=(1,2,3,4,5,6,7))
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	if os.path.exists(galaxiesFile2):
		lambda_R, ellip_gals, star_mis, OIII_mis, Hbeta_mis, Hdelta_mis, \
			Hgamma_mis = np.loadtxt(galaxiesFile2, unpack=True, skiprows=1, 
			usecols=(1,2,3,4,5,6,7))
		galaxy_gals2 = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
		i_gal2 = np.where(galaxy_gals2 == galaxy)[0]
		if len(i_gal2) == 0:
			i_gal2 = -1
			galaxy_gals2 = np.append(galaxy_gals2, galaxy)
			lambda_R = np.append(lambda_R, np.nan)
			ellip_gals = np.append(ellip_gals, np.nan)
			star_mis = np.append(star_mis, np.nan)
			OIII_mis = np.append(OIII_mis, np.nan)
			Hbeta_mis = np.append(Hbeta_mis, np.nan)
			Hdelta_mis = np.append(Hdelta, np.nan)
			Hgamma_mis = np.append(Hgamma_mis, np.nan)
	else:
		galaxy_gals2 = np.array([galaxy])
		lambda_R = np.array([np.nan])
		ellip_gals = np.array([np.nan])
		star_mis = np.array([np.nan])
		OIII_mis = np.array([np.nan])
		Hbeta_mis = np.array([np.nan])
		Hdelta_mis = np.array([np.nan])
		Hgamma_mis = np.array([np.nan])
		i_gal2 = 0

	gas = {'[OIII]5007d':OIII_mis, 'Hbeta':Hbeta_mis, 'Hdelta':Hdelta_mis, 
		'Hgamma':Hgamma_mis}

	R_e = get_R_e(galaxy)
# ------------=============== Photometry =================----------
	save_to = "%s/%s/results/" % (analysis_dir,	galaxy
		) + "%splots/photometry_%s.png" % (wav_range_dir, wav_range)
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed?
	ellip_gals[i_gal2] = f.eps
	print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))
# ------------================ Kinemetry =================----------
	save_to = "%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
		) + "/stellar_kinematics_%s.png" % (wav_range)
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
	star_mis[i_gal2] = mis
	print "Stars misalignment: " + str(mis)
# ------------================ Lambda_R ==================----------
	# distance from center
	R = np.sqrt(np.square(xBar)+np.square(yBar))
	# Angle of r from axis of rotation
	#ang = abs(math.asin(math.sin(math.radians(f.theta))) - np.abs(
	#	np.arctan((xBar)/(yBar))))
	#R = r * np.sin(ang)

	order = np.argsort(R)

	xBar_r = -xBar*math.cos(f.theta)-yBar*math.sin(f.theta)
	yBar_r = xBar*math.sin(f.theta) - yBar*math.cos(f.theta)

	R_m = np.sqrt(np.square(xBar_r)*(1-f.eps) + np.square(yBar_r)/(1+f.eps))
	h = 0.1 # discretization accuracy
		# discretizing the unit square
	x, y = np.mgrid[-f.xmed:40-2*discard-f.xmed:h, -f.ymed:40-2*discard-f.ymed:h]
	for i, r in enumerate(R_m):
			# set all points of ellipse that are inside 
		el = (np.square(-x*math.cos(f.theta)-y*math.sin(f.theta))*(1-f.eps) + \
			np.square(x*math.sin(f.theta) - y*math.cos(f.theta))/(1+f.eps))/r**2 <= 1
		
		if el.any():
			A_s = np.sum(el) * h * h
			A_ellipse = math.pi * r**2 * math.sqrt((1+f.eps**2)/(1-f.eps))
			if 0.85 * A_ellipse < A_s:
				R_m[i] = math.sqrt(A_s/math.pi)
			else:
				R_m[i] = np.nan

	order_m = np.argsort(R_m)

	# NB: lam is ordered in terms of increasing R.
	lam_num = D.flux[order_m]*R[order_m]*np.abs(np.array(
		D.components['stellar'].plot['vel'][order_m])) # numerator
	lam_den = D.flux[order_m]*R[order_m]*np.sqrt(np.square(
		np.array(D.components['stellar'].plot['vel'])[order_m]) + np.square(
		np.array(D.components['stellar'].plot['sigma'])[order_m])) # denominator

	# cumulative summation with mask from R_m
	lam_num = np.cumsum(lam_num[~np.isnan(R_m)])
	lam_den = np.cumsum(lam_den[~np.isnan(R_m)])

	lam = lam_num/lam_den
	plt.figure()
	plt.title(r"Radial $\lambda_R$ profile")
	plt.ylabel(r"$\lambda_R$")
	plt.xlabel("Radius (R_e)")
	x = spxToRe(R[order_m], R_e)[~np.isnan(R_m[order_m])] # Plotted as a fucntion of R not R_m
	order = np.argsort(x)
	lambda_R[i_gal2] = lam[order][-1]
	print 'lambda_R_MAX: ', lambda_R[i_gal2]
	plt.plot(x[order[5:]], lam[order[5:]])
	ax =plt.gca()
	plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), verticalalignment='top',
		transform=ax.transAxes)
	plt.savefig("%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
		) + "lambda_R_%s.png" % (wav_range), bbox_inches="tight")
	if plots: 
		plt.show()
# ------------=================== Gas ====================----------
	for c in D.e_components:
		print c
		save_to = "%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
			) + "%s_kinematics_%s.png" % (c, wav_range)
		gas_k = fit_kinematic_pa(xBar, yBar, np.array(D.components[c].plot['vel']), 
			quiet=True, plot=plots, sav_fig=save_to)

		print "%s_PA_kin: " % (c) + str(gas_k[0]) + "+/-" + str(gas_k[1]/3)
# ------------============== Misalignment ================----------
		gas_kine = math.radians(gas_k[0])

		gas_mis = math.asin(abs(math.sin(gas_kine-kine)))
		gas_mis = math.degrees(gas_mis)

		if 'NI' not in c: gas[c][i_gal2] = gas_mis

		print "Mis-alignment: " + str(gas_mis)
# ------------============== Save outputs ================----------
	template = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}{6:8}{7:8}\n"
	template2 = "{0:13}{1:9}{2:13}{3:15}{4:8}{5:8}{6:8}{7:8}\n" 

	f = open(galaxiesFile, 'wb')
	f.write(template.format("Galaxy", "z", "velocity", "vel dispersion", "x", "y", 
		"Kin SN", "Pop SN"))
	f2 = open(galaxiesFile2, 'wb')
	f2.write(template2.format('Galaxy', 'Lambda_R', 'Ellipticity', 
		'Misa: Stellar', 'OIII', 'Hbeta', 'Hdelta', 'Hgamma'))

	for i in range(len(galaxy_gals)):
		f.write(template.format(galaxy_gals[i], str(z_gals[i]),	str(vel_gals[i]), 
			str(sig_gals[i]), str(int(x_gals[i])), str(int(y_gals[i])), 
			str(SN_kin_gals[i]), str(SN_pop_gals[i])))

	for i in range(len(galaxy_gals2)):
		f2.write(template2.format(galaxy_gals2[i], str(round(lambda_R[i],3)), 
			str(round(ellip_gals[i],3)), str(round(star_mis[i],3)), 
			str(OIII_mis[i]), str(Hbeta_mis[i]), str(Hdelta_mis[i]) ,
			str(Hgamma_mis[i])))

	f.close()
	f2.close()

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc3557'
	discard = 2
	wav_range = '4200-'

	kinematics(galaxy, discard=discard, wav_range=wav_range, plots=False)
