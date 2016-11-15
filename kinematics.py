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
import re # for regex expressions

#---------------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, discard=0, wav_range="", plots=False):

	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)

	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""
		wav_range = ""
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	galaxiesFile2 = "%s/galaxies2.txt" % (analysis_dir)
	classify_file = "%s/galaxies_classify.txt" % (analysis_dir)
	galaxiesFile_Re =  "%s/galaxies_R_e.txt" % (analysis_dir)
	output = '%s/%s/results/%s' % (analysis_dir, galaxy, wav_range_dir)
	pickleFile = open('%s/pickled/dataObj_%s.pkl' % (output, wav_range))
	D = pickle.load(pickleFile)
	pickleFile.close()


	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_gals = np.loadtxt(galaxiesFile, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5,6))
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	lambda_R, ellip_gals, star_mis, OIII_mis, Hbeta_mis, Hdelta_mis, \
		Hgamma_mis, k5k1 = np.loadtxt(galaxiesFile2, unpack=True, skiprows=1, 
		usecols=(1,2,3,4,5,6,7,8))
	gas = {'[OIII]5007d':OIII_mis, 'Hbeta':Hbeta_mis, 'Hdelta':Hdelta_mis, 'Hgamma':Hgamma_mis}


	log_R_e_RC3_gals, R_e_2MASS_gals = np.loadtxt(galaxiesFile_Re, unpack=True, 
		skiprows=1, usecols=(1,2))
	R_e_RC3 = 6*10**log_R_e_RC3_gals[i_gal]/2 # convert to arcsec
	R_e_2MASS = R_e_2MASS_gals[i_gal]

	R_e = np.nanmean([R_e_RC3,R_e_2MASS])
# ------------=============== Photometry =================----------
	save_to = "%s/%s/results/" % (analysis_dir,	galaxy
		) + "%splots/photometry_%s.png" % (wav_range_dir, wav_range)
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed ?
	ellip_gals[i_gal] = f.eps
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
	star_mis[i_gal] = mis
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
	lambda_R[i_gal] = lam[order][-1]
	print 'lambda_R_MAX: ', lambda_R[i_gal]
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

		if 'NI' not in c: gas[c][i_gal] = gas_mis

		#gas_psi_gals[i_gal] = gas_mis
		print "Mis-alignment: " + str(gas_mis)

# ------------================= RR/NRR ===================----------
	RR, NF, NR, KT, M2, KDC = np.loadtxt(classify_file, skiprows=1, unpack=True,
		usecols=(1,2,3,4,5,6), dtype='S3')

	file = '%s/%s/kinemetry_vel.txt' % (analysis_dir,galaxy)
	rad, pa, k1, k51 = np.loadtxt(file, usecols=(0,1,5,7), skiprows=1, unpack=True)
	rad *= 0.67 # Pix to arcsec


	# Finding smoothest pa by add or subtracting 360 deg	
	for j in range(1,len(pa)):
	    test = np.array([pa[j]-360, pa[j], pa[j]+360])
	    a = np.argmin(np.abs(test-pa[j-1]))
	    if a==0: 
	        pa[j:]-=360
	    elif a==2: 
	        pa[j:]+=360

	###### NEEDS TO BE CHANGED TO LUMINOSITY WEIGHTED MEAN
	k5k1[i_gal] = np.average(k51/k1, weights=np.ones(len(k1)))
	if k5k1[i_gal] < 0.04:
		RR[i_gal] = 'RR'
	else: RR[i_gal] = 'NRR'
# ------------=========== kinematic features =============----------
	i_rad = np.where(rad < R_e)[0][-1] # Largest radius below R_e
	feature = False

	# Non-rotators
	if k1[i_rad] < 5: 
		NR[i_gal] = 'NR'
		feature = True
	else:
		NR[i_gal] = '-'

	# kinematcs twist
	# smooth: are 80% of pionts increaded or decreased from their neighbour?
	difference = np.append(pa[:-1]-pa[1:], 0) 
	change = difference < 0
	smooth = 0.7 	# percentage of points to be moving the same direction to be 
					# considered smooth
	if np.sum(change) > smooth*len(pa):
		twist = abs(np.arcsin(np.sin(max(pa[change])-min(pa[change]))))
		if twist > 10: 
			KT[i_gal] = 'KT'
			feature= True
		else: KT[i_gal] = '-'
	elif np.sum(change) < (1-smooth)*len(pa):
		twist = abs(np.arcsin(np.sin(max(pa[~change])-min(pa[~change]))))
		if twist > 10: 
			KT[i_gal] = 'KT'
			feature = True
		else: KT[i_gal] = '-'
	else: KT[i_gal] = '-'

	# Kinematically distinct/decoupled core
	# CRC not included yet
	sharp = np.abs(difference) > 30
	kdc_location = np.logical_and(sharp, k1 < 0.15*max(k1))
	if any(kdc_location) and any(difference[0:np.median(np.where(kdc_location)[0])]<3):
		KDC[i_gal] = 'KDC'
		feature = True
	else: KDC[i_gal] = '-'

	# Double maximum (M2)
	k1_change = np.append(pa[1:-1]-pa[0:-2], 0) > 0
	k1_change = str(k1_change.astype(int))
	# at least 3 True's followed by at least 2 False' followed by at least 3 True's
	if re.search('1 1( 1)+ 0( 0)+ 1 1( 1)+', k1_change) is not None:
		M2[i_gal] = '2M'
		feature = True
	else: M2[i_gal] = '-'

	# No features
	if not feature: NF[i_gal] = 'NF' 
	else: NF[i_gal] = '-'
# ------------============== Save outputs ================----------
	template = "{0:13}{1:14}{2:15}{3:16}{4:4}{5:4}{6:11}\n" 
	template2 = "{0:13}{1:9}{2:13}{3:15}{4:8}{5:8}{6:8}{7:8}{8:8}\n" 
	template3 = "{0:13}{1:7}{2:5}{3:5}{4:5}{5:5}{6:5}\n"

	f = open(galaxiesFile, 'wb')
	f.write(template.format('Galaxy', 'z', 'velocity', 'vel dispersion', 'x', 
		'y', 'Target SN'))
	f2 = open(galaxiesFile2, 'wb')
	f2.write(template2.format('Galaxy', 'Lambda_R', 'Ellipticity', 
		'Misa: Stellar', 'OIII', 'Hbeta', 'Hdelta', 'Hgamma','k5/k1'))

	f3 = open(classify_file, 'wb')
	f3.write(template3.format('Galaxy', 'RR/NRR', 'NF', 'NR', 'KT', '2M', 'KDC'))

	for i in range(len(galaxy_gals)):
		f.write(template.format(galaxy_gals[i], str(z_gals[i]),	str(vel_gals[i]), 
			str(sig_gals[i]), str(int(x_gals[i])), str(int(y_gals[i])), 
			str(SN_gals[i])))

		f2.write(template2.format(galaxy_gals[i], str(round(lambda_R[i],3)), 
			str(round(ellip_gals[i],3)), str(round(star_mis[i],3)), 
			str(OIII_mis[i]), str(Hbeta_mis[i]), str(Hdelta_mis[i]) ,
			str(Hgamma_mis[i]), str(k5k1[i])))

		f3.write(template3.format(galaxy_gals[i], RR[i], NF[i], NR[i], KT[i], 
			M2[i], KDC[i]))
	f.close()
	f2.close()
	f3.close()

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc3557'
	discard = 2
	wav_range = '4200-'

	kinematics(galaxy, discard=discard, wav_range=wav_range, plots=True)
