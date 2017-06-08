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
			for i in range(D.number_of_bins):
				f.write(str(D.flux[i]) + '\n')

		with open('%s/kinemetry/vel.dat' % (output), 'wb') as f:
			for i in range(D.number_of_bins):
				f.write(str(D.components['stellar'].plot['vel'][i]) + '  ' + 
					str(D.components['stellar'].plot['vel'].uncert[i]) + '\n')

		with open('%s/kinemetry/sigma.dat' % (output), 'wb') as f:
			for i in range(D.number_of_bins):
				f.write(str(D.components['stellar'].plot['sigma'][i]) + '  ' + 
					str(D.components['stellar'].plot['sigma'].uncert[i]) + '\n')

		return D