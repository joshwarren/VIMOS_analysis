## ==================================================================
## 		Save to text file for KINEMTRY IDL
## ==================================================================
from checkcomp import checkcomp
cc = checkcomp()
import cPickle as pickle
import os
import numpy as np


def sav_for_kinemetry(galaxy, opt='kin', D=None):
	if 'kin' in opt:
		print "    Saving for KINEMETRY (IDL)"

		output = '%s/Data/vimos/analysis/%s/%s' % (cc.base_dir, galaxy, opt)
		if not os.path.exists('%s/kinemetry' % (output)):
			os.makedirs('%s/kinemetry' % (output))

		if D is None:
			pickleFile = open("%s/pickled/dataObj.pkl" % (output), 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()


		if D.norm_method != 'lws':
			D.norm_method = 'lws'
			D.find_restFrame()


		# Stellar velocity (for classifying)
		with open('%s/kinemetry/stellar_vel.dat' % (output), 'wb') as f:
			vel = D.components['stellar'].plot['vel']
			for i in range(D.number_of_bins):
				f.write(str(vel[i]) + '  ' + str(vel.uncert[i]) + '\n')


		# Gas kinematics for deviation from circular
		with open('%s/kinemetry/gas_flux.dat' % (output), 'wb') as f:
			flux = D.gas_flux
			flux[np.isnan(flux)] = 9999
			for i in range(D.number_of_bins):
				f.write(str(flux[i]) + '\n')

		with open('%s/kinemetry/gas_vel.dat' % (output), 'wb') as f:
			vel = D.components['Hbeta'].plot['vel']
			vel[np.isnan(vel)] = 9999
			vel.uncert[np.isnan(vel.uncert)] = 9999
			for i in range(D.number_of_bins):
				f.write(str(vel[i]) + '  ' + str(vel.uncert[i]) + '\n')

		with open('%s/kinemetry/gas_sigma.dat' % (output), 'wb') as f:
			sigma = D.components['Hbeta'].plot['sigma']
			sigma[np.isnan(sigma)] = 9999
			sigma.uncert[np.isnan(sigma.uncert)] = 9999
			for i in range(D.number_of_bins):
				f.write(str(sigma[i]) + '  ' + str(sigma.uncert[i]) + '\n')

		return D

if __name__=='__main__':
	gals = ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 
		'ngc3557', 'ngc7075', 'pks0718-34']
	for g in gals:
		sav_for_kinemetry(g, opt='kin')