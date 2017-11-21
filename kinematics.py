## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import os
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import matplotlib.pyplot as plt # used for plotting
# from scipy.optimize import curve_fit # for fitting a gaussian
from scipy.interpolate import interp1d
from checkcomp import checkcomp
from classify import get_R_e
cc = checkcomp()
import cPickle as pickle
from sav_for_kinemetry import find_mask

#----------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, discard=0, opt="kin", plots=False, D=None, 
	instrument='vimos'):

	print '  kinematics'

	analysis_dir = "%s/Data/%s/analysis" % (cc.base_dir, instrument)

	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	galaxiesFile2 = "%s/galaxies2.txt" % (analysis_dir)
	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)
	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()


	d = np.loadtxt(galaxiesFile, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	if instrument == 'vimos':
		x0, y0 = int(d[4][i_gal + 1]), int(d[5][i_gal + 1])
	elif instrument == 'muse':
		x0, y0 = int(d[1][i_gal + 1]), int(d[2][i_gal + 1])

	if os.path.exists(galaxiesFile2):
		lambda_Re_gals, ellip_gals, pa_gals, e_pa_gals, \
			star_kine_pa_gals, e_star_kine_pa_gals, \
			gas_kine_pa_gals, e_gas_kine_pa_gals = np.loadtxt(
			galaxiesFile2, unpack=True, skiprows=1, 
			usecols=(1,2,3,4,5,6,7,8))
		galaxy_gals2 = np.loadtxt(galaxiesFile2, skiprows=1, 
			usecols=(0,), dtype=str)
		i_gal2 = np.where(galaxy_gals2 == galaxy)[0]
		if len(i_gal2) == 0:
			i_gal2 = -1
			galaxy_gals2 = np.append(galaxy_gals2, galaxy)
			lambda_Re_gals = np.append(lambda_Re_gals, np.nan)
			ellip_gals = np.append(ellip_gals, np.nan)
			pa_gals = np.append(pa_gals, np.nan)
			e_pa_gals = np.append(e_pa_gals, np.nan)
			star_kine_pa_gals = np.append(star_kine_pa_gals, np.nan)
			e_star_kine_pa_gals = np.append(e_star_kine_pa_gals, np.nan)
			gas_kine_pa_gals = np.append(gas_kine_pa_gals, np.nan)
			e_gas_kine_pa_gals = np.append(e_gas_kine_pa_gals, np.nan)
		else: i_gal2 = i_gal2[0]
		if type(lambda_Re_gals) is np.float64:
			i_gal2 = 0
			galaxy_gals2 = np.array([galaxy_gals2])
			lambda_Re_gals = np.array([lambda_Re_gals])
			ellip_gals = np.array([ellip_gals])
			pa_gals = np.array([pa_gals])
			e_pa_gals = np.array([e_pa_gals])
			star_kine_pa_gals = np.array([star_kine_pa_gals])
			e_star_kine_pa_gals = np.array([e_star_kine_pa_gals])
			gas_kine_pa_gals = np.array([gas_kine_pa_gals])
			e_gas_kine_pa_gals = np.array([e_gas_kine_pa_gals])
	else:
		galaxy_gals2 = np.array([galaxy])
		lambda_Re_gals = np.array([np.nan])
		ellip_gals = np.array([np.nan])
		pa_gals = np.array([np.nan])
		e_pa_gals = np.array([np.nan])
		star_kine_pa_gals = np.array([np.nan])
		e_star_kine_pa_gals = np.array([np.nan])
		gas_kine_pa_gals = np.array([np.nan])
		e_gas_kine_pa_gals = np.array([np.nan])
		i_gal2 = 0

	R_e = get_R_e(galaxy)
# ------------=============== Photometry =================----------
	file = '%s/kinemetry/kinemetry_stellar_flux.txt' % (output)
	a, pa_r, e_pa_r, ellip, e_ellip, b4, e_b4 = np.loadtxt(file, 
		unpack=True, skiprows=1)

	A_ellipse = np.pi * a**2 * (1 - ellip)
	A_s = np.zeros(len(a))

	if instrument == 'vimos':
		res = 0.67 # arcsec/pix
		l = 40
	elif instrument == 'muse':
		res = 0.2 # arcsec/pix
		l = 150

	h = 0.1 # pix; discretization accuracy
	# discretizing FoV
	x, y = np.meshgrid(np.arange(0, l, h), np.arange(0, l, h))
	x -= x0
	y -= y0
	x *= res 
	y *= res

	for i in range(len(a)):
		# set all points of ellipse that are inside
		x_rot = x * np.cos(np.radians(90 + pa_r[i])) - \
			y * np.sin(np.radians(90 + pa_r[i]))
		y_rot = y * np.cos(np.radians(90 + pa_r[i])) + \
			x * np.sin(np.radians(90 + pa_r[i]))
		el = ((x_rot)/a[i])**2 + ((y_rot)/
			(a[i] * (1 - ellip[i])))**2 <= 1 
		# plt.imshow(el.astype(int))
		# plt.show()
		A_s[i] = np.sum(el) * h**2 * res**2

		if A_s[i] > 0.85 * A_ellipse[i]:
			R_max = a[i]

	R_m = np.sqrt(A_s/np.pi)
	m = np.array(R_m <= R_max)

	R_m = R_m[m]
	ellip, e_ellip, pa_r, e_pa_r, b4, e_b4 = ellip[m], e_ellip[m], \
		pa_r[m], e_pa_r[m], b4[m], e_b4[m]

	ellip_gals[i_gal2] = interp1d(R_m, ellip, bounds_error=False, 
		fill_value=(0, ellip[-1]))(R_e)
	pa_gals[i_gal2] = interp1d(R_m, pa_r, bounds_error=False, 
		fill_value=(0, pa_r[-1]))(R_e)
	e_pa_gals[i_gal2] = interp1d(R_m, e_pa_r, bounds_error=False, 
		fill_value=(0, e_pa_r[-1]))(R_e)
# ------------================ Lambda_R ==================----------
	m = find_mask(galaxy, instrument, D)

	vel = D.components['stellar'].plot['vel']
	sigma = D.components['stellar'].plot['sigma']

	vel[~m] = np.nan
	sigma[~m] = np.nan

	R = np.sqrt((D.xBar - x0)**2 + (D.yBar - y0)**2) * res

	lambda_R = np.zeros(len(R_m))

	for i in range(len(R_m)):
		x_rot = (D.xBar - x0) * res * np.cos(np.radians(90 + pa_r[i])) \
			- (D.yBar - y0) * res * np.sin(np.radians(90 + pa_r[i]))
		y_rot = (D.yBar - y0) * res * np.cos(np.radians(90 + pa_r[i])) \
			+ (D.xBar - x0) * res * np.sin(np.radians(90 + pa_r[i]))
		el = ((x_rot)/a[i])**2 + ((y_rot)/
			(a[i] * (1 - ellip[i])))**2 <= 1

		numerator = np.nansum(D.flux[el] * R[el] * np.abs(vel[el]))
		denominator = np.nansum(D.flux[el] * R[el] * 
			np.sqrt(vel[el]**2 + sigma[el]**2))

		lambda_R[i] = numerator/denominator

	lambda_Re = interp1d(R_m, lambda_R, bounds_error=False, 
		fill_value=(0, lambda_R[-1]))(R_e)

	print 'lambda_Re: ', lambda_Re

	lambda_Re_gals[i_gal2] = lambda_Re

	file = '%s/lambda_R.txt' % (output)

	with open(file, 'w') as f2:
		for i in range(len(R_m)):
			f2.write('%.4f   %.4f \n' %(R_m[i], lambda_R[i]))

	fig, ax = plt.subplots()
	ax.set_title(r"Radial $\lambda_R$ profile")
	ax.set_ylabel(r"$\lambda_R$")
	ax.set_xlabel(r"Radius ($R_m$/$R_e$)")
	ax.plot(R_m/R_e, lambda_R)

	if not os.path.exists('%s/plots' % (output)):
		os.makedirs('%s/plots' % (output))
	plt.savefig("%s/plots/lambda_R.png" % (output), bbox_inches="tight")
	if plots: 
		plt.show()
# # ------------========= Stellar Kinematics ===============----------
	save_to = "%s/plots/stellar_kinematics.png" % (output)
	star_kine_pa_gals[i_gal2], e_star_kine_pa_gals[i_gal2], vSyst = \
		fit_kinematic_pa(D.xBar-x0, D.yBar-y0, vel, quiet=True, 
		dvel=vel.uncert, plot=plots, sav_fig=save_to, vSyst=0)
# # ------------=========== Gas Kinematics =================----------
# # NB: this is not written for gas=2 or gas=3 options. 
# 	if D.gas == 1:
# 		save_to = "%s/plots/gas_kinematics.png" % (output)
# 		kgas = fit_kinematic_pa(D.xBar - f.xpeak, D.yBar - f.ypeak, 
# 			np.array(D.components['Hbeta'].plot['vel']), quiet=True, 
# 			plot=plots, sav_fig=save_to)
# 		gas_kine_pa_gals[i_gal2] = kgas[0]
# ------------============== Save outputs ================----------
	template2 = "{0:13}{1:9}{2:13}{3:9}{4:8}{5:17}{6:8}{7:8}{8:8}\n" 
	with open(galaxiesFile2, 'wb') as f2:
	# f2 = open('/Data/vimos/analysis/test.txt', 'wb')
		f2.write(template2.format('Galaxy', 'Lambda_R', 'Ellipticity', 
			'PA (deg)', 'e_PA', 'kine_pa: Stellar', 'e_star', 'Gas', 
			'e_gas'))
		for i in range(len(galaxy_gals2)):
			f2.write(template2.format(galaxy_gals2[i], 
				str(round(lambda_Re_gals[i],3)), 
				str(round(ellip_gals[i],3)), str(round(pa_gals[i],3)), 
				str(round(e_pa_gals[i],3)), 
				str(round(star_kine_pa_gals[i],3)), 
				str(round(e_star_kine_pa_gals[i],3)), 
				str(round(gas_kine_pa_gals[i],3)), 
				str(round(e_gas_kine_pa_gals[i],3))))

	return D

#####################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'pks0718-34'
	discard = 0

	kinematics(galaxy, discard=discard, plots=False)
