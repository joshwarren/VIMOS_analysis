## ==================================================================
## 		Save to text file for KINEMTRY IDL
## ==================================================================
from checkcomp import checkcomp
cc = checkcomp()
import cPickle as pickle
import os


def sav_for_kinemetry(galaxy, opt='kin', D=None):
	if 'kin' in opt:

		output = '%s/Data/vimos/analysis/%s/%s' % (cc.base_dir, galaxy, opt)
		if not os.path.exists('%s/kinemetry' % (output)):
			os.makedirs('%s/kinemetry' % (output))

		if D is None:
			pickleFile = open("%s/pickled/dataObj.pkl" % (output), 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()


		print "    Saving flux for KINEMETRY (IDL)"
		with open('%s/kinemetry/flux.dat' % (output), 'wb') as f:
			flux  = D.flux
			for i in range(D.number_of_bins):
				f.write(str(flux[i]) + '\n')

		with open('%s/kinemetry/vel.dat' % (output), 'wb') as f:
			vel = D.components['stellar'].plot['vel']
			for i in range(D.number_of_bins):
				f.write(str(vel[i]) + '  ' + str(vel.uncert[i]) + '\n')

		with open('%s/kinemetry/sigma.dat' % (output), 'wb') as f:
			sig = D.components['stellar'].plot['sigma']
			for i in range(D.number_of_bins):
				f.write(str(sig[i]) + '  ' + str(sig.uncert[i]) + '\n')

		return D