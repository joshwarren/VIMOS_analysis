## ==================================================================
## 		Save to text file for KINEMTRY IDL
## ==================================================================
from checkcomp import checkcomp
cc = checkcomp()
import cPickle as pickle
import os
import numpy as np

def find_mask(galaxy, instrument, D):
	data_file =  "%s/Data/%s/analysis/galaxies.txt" % (cc.base_dir,
		instrument)
	if instrument == 'vimos':
		cols = (4,5)
	elif instrument == 'muse':
		cols = (1,2)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=cols, dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),
		dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	center = np.array([x_cent_gals[i_gal], y_cent_gals[i_gal]])

	m = np.ones(D.number_of_bins).astype(bool)
	if instrument == 'muse' and galaxy == 'ngc1399':
		size = 1.7/0.2
		star_center = np.array([1.7, 17.5])/0.2
		m = np.sqrt((D.xBar - center[0] - star_center[0])**2 + (
			D.yBar - center[1] - star_center[1])**2) >= size

	if galaxy == 'pks0718-34':
		size = 3/0.67
		star_center = np.array([1, 10])/0.67
		m = np.sqrt((D.xBar - center[0] - star_center[0])**2 + (
			D.yBar - center[1] - star_center[1])**2) >= size

	return m

def sav_for_kinemetry(galaxy, opt='kin', D=None, instrument='vimos'):
	if 'kin' in opt:
		print "    Saving for KINEMETRY (IDL)"

		output = '%s/Data/%s/analysis/%s/%s' % (cc.base_dir, instrument,
			galaxy, opt)
		if not os.path.exists('%s/kinemetry' % (output)):
			os.makedirs('%s/kinemetry' % (output))

		if D is None:
			pickleFile = open("%s/pickled/dataObj.pkl" % (output), 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

		if D.norm_method != 'lws':
			D.norm_method = 'lws'
			D.find_restFrame()

		m = find_mask(galaxy, instrument, D)

		# Stellar velocity (for classifying)
		with open('%s/kinemetry/stellar_vel.dat' % (output), 'wb') as f:
			vel = D.components['stellar'].plot['vel']
			e_v = np.array(vel.uncert)
			v = np.array(vel)
			v[~m] = 9999

			for i in range(D.number_of_bins):
				f.write(str(v[i]) + '  ' + str(e_v[i]) 
					+ '\n')


		# Stellar flux (for properties at R_e)
		with open('%s/kinemetry/stellar_flux.dat' % (output), 'wb') as f:
			flux = D.flux
			e_f = np.array(flux.uncert)
			fl = np.array(flux)
			fl[~m] = 9999
			for i in range(D.number_of_bins):
				f.write(str(fl[i]) + '  ' + str(e_f[i]) 
					+ '\n')


		# Gas kinematics for deviation from circular
		# with open('%s/kinemetry/gas_flux.dat' % (output), 'wb') as f:
		# 	flux = D.gas_flux
		# 	flux[np.isnan(flux)] = 9999
		# 	for i in range(D.number_of_bins):
		# 		f.write(str(flux[i]) + '\n')

		# with open('%s/kinemetry/gas_vel.dat' % (output), 'wb') as f:
		# 	vel = D.components['[OIII]5007d'].plot['vel']
		# 	vel[np.isnan(vel)] = 9999
		# 	vel.uncert[np.isnan(vel.uncert)] = 9999
		# 	for i in range(D.number_of_bins):
		# 		f.write(str(vel[i]) + '  ' + str(vel.uncert[i]) + '\n')

		# with open('%s/kinemetry/gas_sigma.dat' % (output), 'wb') as f:
		# 	sigma = D.components['[OIII]5007d'].plot['sigma']
		# 	sigma[np.isnan(sigma)] = 9999
		# 	sigma.uncert[np.isnan(sigma.uncert)] = 9999
		# 	for i in range(D.number_of_bins):
		# 		f.write(str(sigma[i]) + '  ' + str(sigma.uncert[i]) + 
		# 		'\n')

		return D

if __name__=='__main__':
	gals = ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']
	gals = ['pks0718-34']
	for g in gals:
		sav_for_kinemetry(g, opt='kin')