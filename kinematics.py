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
from scipy.interpolate import interp1d
from checkcomp import checkcomp
from classify import get_R_e
cc = checkcomp()
import cPickle as pickle
from plot_results import set_lims

#---------------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, discard=0, opt="kin", plots=False, D=None):
	print '  kinematics'

	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)

	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	galaxiesFile2 = "%s/galaxies2.txt" % (analysis_dir)
	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)
	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()


	d = np.loadtxt(galaxiesFile, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	z_gals, vel_gals, sig_gals = d[1][1:].astype(float), d[2][1:].astype(float), \
		d[3][1:].astype(float),
	x_gals, y_gals = d[4][1:].astype(int), d[5][1:].astype(int)
	SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(6,len(d))}
	i_gal = np.where(galaxy_gals==galaxy)[0][0]


	if os.path.exists(galaxiesFile2):
		lambda_Re_gals, ellip_gals, pa_gals, star_mis, OIII_mis, Hbeta_mis, \
			Hdelta_mis, Hgamma_mis = np.loadtxt(galaxiesFile2, unpack=True, 
			skiprows=1, usecols=(1,2,3,4,5,6,7,8))
		galaxy_gals2 = np.loadtxt(galaxiesFile2, skiprows=1, usecols=(0,),
			dtype=str)
		i_gal2 = np.where(galaxy_gals2 == galaxy)[0]
		if len(i_gal2) == 0:
			i_gal2 = -1
			galaxy_gals2 = np.append(galaxy_gals2, galaxy)
			lambda_Re_gals = np.append(lambda_Re_gals, np.nan)
			ellip_gals = np.append(ellip_gals, np.nan)
			pa_gals = np.append(pa_gals, np.nan)
			star_mis = np.append(star_mis, np.nan)
			OIII_mis = np.append(OIII_mis, np.nan)
			Hbeta_mis = np.append(Hbeta_mis, np.nan)
			Hdelta_mis = np.append(Hdelta_mis, np.nan)
			Hgamma_mis = np.append(Hgamma_mis, np.nan)
		else: i_gal2 = i_gal2[0]
		if type(lambda_Re_gals) is np.float64:
			i_gal2 = 0
			galaxy_gals2 = np.array([galaxy_gals2])
			lambda_Re_gals = np.array([lambda_Re_gals])
			ellip_gals = np.array([ellip_gals])
			pa_gals = np.array([pa_gals])
			star_mis = np.array([star_mis])
			OIII_mis = np.array([OIII_mis])
			Hbeta_mis = np.array([Hbeta_mis])
			Hdelta_mis = np.array([Hdelta_mis])
			Hgamma_mis = np.array([Hgamma_mis])
	else:
		galaxy_gals2 = np.array([galaxy])
		lambda_Re_gals = np.array([np.nan])
		ellip_gals = np.array([np.nan])
		pa_gals = np.array([np.nan])
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
	# save_to = "%s/plots/photometry.png" % (output)
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed?
	ellip_gals[i_gal2] = f.eps
	pa_gals[i_gal2] = f.theta

	print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))
# ------------================ Kinemetry =================----------
# 	save_to = "%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
# 		) + "/stellar_kinematics_%s.png" % (wav_range)
#  	xBar = np.array(D.xBar)-f.xmed
#  	yBar = np.array(D.yBar)-f.ymed


# 	k = fit_kinematic_pa(xBar, yBar, np.array(D.components['stellar'].plot['vel']), 
# 		quiet=True, plot=plots, sav_fig=save_to)

# 	print "PA_kin: " + str(k[0]) + "+/-" + str(k[1]/3)
# ------------============== Misalignment ================----------
# 	phot = math.radians(90-f.theta)
# 	kine = math.radians(k[0])


# 	mis = math.asin(abs(math.sin(phot-kine)))
# 	mis = math.degrees(mis)
# 	star_mis[i_gal2] = mis
# 	print "Stars misalignment: " + str(mis)
# ------------================ Lambda_R ==================----------

	beta = np.tan(D.yBar/D.xBar) # Angle between point and center of galaxy and RA 
									# axis.
	pa = 90 - f.theta # Position angle
	res = 0.67 # arcsec. Spatial resolution of MUSE

	# Semi-major axis
	a = res* np.sqrt(((D.xBar - f.xpeak)**2 + (D.yBar - f.ypeak)**2) * \
		(np.sin(beta + pa)**2 + np.cos(beta + pa)**2/(1 - f.eps)**2))

	R_m = a * np.sqrt(1 - f.eps)
	R_m_sort = R_m.argsort() # argsort of R_m_sort will be the inverse opperation 
	
	# Area of ellipse
	A_ellipse = np.pi * a**2 * (1 - f.eps)
	# Area sampled
	A_s = np.cumsum(D.n_spaxels_in_bin[R_m_sort])[R_m_sort.argsort()] * res**2

	R_m[A_ellipse > A_s] = np.sqrt(A_s[A_ellipse>A_s]/np.pi)
	# R_max occurs when 0.85*A_ellipse = A_s
	# R_m[0.85*A_ellipse > A_s] = np.nan

	vel = D.components['stellar'].plot['vel']
	vel_lim = set_lims(vel, symmetric=True)
	mask = (vel < vel_lim[0]) + (vel > vel_lim[1])

	sigma = D.components['stellar'].plot['sigma']
	sigma_lim = set_lims(sigma, positive=True)
	mask += (sigma < sigma_lim[0]) + (sigma > sigma_lim[1])
	
	vel[mask], sigma[mask] = np.nan, np.nan


	# NB: numerator and denominator are in R_m order
	numerator = np.nancumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * np.abs(vel[R_m_sort]))

	denominator = np.nancumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * np.sqrt(vel**2 + sigma**2)[R_m_sort])

	lambda_R = numerator[R_m_sort.argsort()]/denominator[R_m_sort.argsort()]
	lambda_R[np.isnan(R_m)] = np.nan

	lambda_Re = interp1d(R_m[~np.isnan(R_m)], lambda_R[~np.isnan(R_m)],
		bounds_error=False, fill_value=(0, 
			lambda_R[R_m_sort][~np.isnan(lambda_R[R_m_sort])][-1]))(R_e)

	print 'lambda_Re: ', lambda_Re
	lambda_Re_gals[i_gal2] = lambda_Re

	fig, ax = plt.subplots()
	ax.set_title(r"Radial $\lambda_R$ profile")
	ax.set_ylabel(r"$\lambda_R$")
	ax.set_xlabel(r"Radius ($R_m$/$R_e$)")
	ax.plot(R_m[R_m_sort][10:]/R_e, lambda_R[R_m_sort][10:])
	# ax.text(0.02,0.98, "Galaxy: " + galaxy.upper())#, verticalalignment='top',
		# transform=ax.transAxes)
	plt.savefig("%s/plots/lambda_R.png" % (output), bbox_inches="tight")
	if plots: 
		plt.show()
# ------------=================== Gas ====================----------
# 	for c in D.e_components:
# 		print c
# 		save_to = "%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
# 			) + "%s_kinematics_%s.png" % (c, wav_range)
# 		gas_k = fit_kinematic_pa(xBar, yBar, np.array(D.components[c].plot['vel']), 
# 			quiet=True, plot=plots, sav_fig=save_to)

# 		print "%s_PA_kin: " % (c) + str(gas_k[0]) + "+/-" + str(gas_k[1]/3)
# # ------------============== Misalignment ================----------
# 		gas_kine = math.radians(gas_k[0])

# 		gas_mis = math.asin(abs(math.sin(gas_kine-kine)))
# 		gas_mis = math.degrees(gas_mis)

# 		if 'NI' not in c: gas[c][i_gal2] = gas_mis

# 		print "Mis-alignment: " + str(gas_mis)
# ------------============== Save outputs ================----------
	template = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}" + ''.join(
		['{%i:%i}'%(i+6,len(t)+1) for i, t in enumerate(SN_gals.keys())]
		) + "\n"

	SN_titles = list(SN_gals.keys())
	with open(galaxiesFile, 'w') as f:
		f.write(template.format("Galaxy", "z", "velocity", "sigma", "x", "y", 
			*(s for s in SN_titles)))
		for i in range(len(galaxy_gals)):
			f.write(template.format(galaxy_gals[i], str(round(z_gals[i],7)), 
				str(round(vel_gals[i],4)), str(round(sig_gals[i],4)), 
				str(int(x_gals[i])), str(int(y_gals[i])), 
				*(str(round(SN_gals[s][i],2)) for s in SN_titles)))

	
	template2 = "{0:13}{1:9}{2:13}{3:13}{4:15}{5:8}{6:8}{7:8}{8:8}\n" 
	with open(galaxiesFile2, 'wb') as f2:
	# f2 = open('/Data/vimos/analysis/test.txt', 'wb')
		f2.write(template2.format('Galaxy', 'Lambda_R', 'Ellipticity', 'PA (deg)',
			'Misa: Stellar', 'OIII', 'Hbeta', 'Hdelta', 'Hgamma'))
		for i in range(len(galaxy_gals2)):
			f2.write(template2.format(galaxy_gals2[i], 
				str(round(lambda_Re_gals[i],3)), str(round(ellip_gals[i],3)), 
				str(round(pa_gals[i],3)), str(round(star_mis[i],3)), 
				str(OIII_mis[i]), str(Hbeta_mis[i]), str(Hdelta_mis[i]),
				str(Hgamma_mis[i])))

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'eso443-g024'
	discard = 0

	kinematics(galaxy, discard=discard, plots=False)
