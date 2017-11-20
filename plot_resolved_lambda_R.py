import numpy as np 
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
from errors2 import get_dataCubeDirectory
from classify import get_R_e
from prefig import Prefig
Prefig()

import cPickle as pickle

opt = 'kin'

def plot_lambda_R():
	fig, ax = plt.subplots()
	FR = ['ic1531', 'ngc0612', 'ngc1316', 'ngc3100', 'ngc3557', 
		'pks0718-34']
	for galaxy in ['eso443-g024',
				'ic1459',
				'ic1531', 
				'ic4296',
				'ngc0612',
				'ngc1399',
				'ngc3100',
				'ngc3557',
				'ngc7075',
				'pks0718-34']:
		lam_R_file = '%s/Data/vimos/analysis/%s/%s/lambda_R.txt' % (
			cc.base_dir, galaxy, opt)
		r, lam_R = np.loadtxt(lam_R_file, unpack=True)
		r /= get_R_e(galaxy)

		if galaxy in FR:
			ax.plot(r, lam_R, 'r')
		else:
			ax.plot(r, lam_R, 'r--')
		ax.text(r[-1], lam_R[-1], galaxy.upper(), color='r', 
			va='center', size='x-small')

	for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		lam_R_file = '%s/Data/muse/analysis/%s/%s/lambda_R.txt' % (
			cc.base_dir, galaxy, opt)
		r, lam_R = np.loadtxt(lam_R_file, unpack=True)
		r /= get_R_e(galaxy)

		if galaxy in FR:
			ax.plot(r, lam_R, 'b')
		else:
			ax.plot(r, lam_R, 'b--')
		ax.text(r[-1], lam_R[-1], galaxy.upper(), color='b', 
			va='center', size='x-small')

	ax.set_xlabel(r'Radius (R$_e$)')
	ax.set_ylabel(r'$\lambda_R$')

	fig.savefig('%s/Data/vimos/analysis/lambda_R.png' % (cc.base_dir))


if __name__=='__main__':
	plot_lambda_R()






















