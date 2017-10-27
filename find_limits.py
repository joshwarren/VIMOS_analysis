from astropy.io import fits # reads fits files (is from astropy)
import numpy as np # for array handling
import cPickle as pickle
from checkcomp import checkcomp
from plot_results import set_lims
cc = checkcomp()





def find_limits(galaxy, opt='kin', norm='fit_disk', D=None, 
	instrument='vimos', plots=None):

	if plots is None:
		if 'kin' in opt:
			plots = [
				'flux',
				"components['stellar'].plot['vel']",
				"components['stellar'].plot['sigma']",
				"components['stellar'].plot['vel'].uncert",
				"components['stellar'].plot['sigma'].uncert"
				]
		elif 'pop' in opt:
			plots = [
				"absorption_line('G4300')",
				"absorption_line('Fe4383')",
				"absorption_line('Ca4455')",
				"absorption_line('Fe4531')",
				"absorption_line('H_beta')",
				"absorption_line('Fe5015')",
				"absorption_line('Mg_b')",
				"absorption_line('G4300', uncert=True)[1]",
				"absorption_line('Fe4383', uncert=True)[1]",
				"absorption_line('Ca4455', uncert=True)[1]",
				"absorption_line('Fe4531', uncert=True)[1]",
				"absorption_line('H_beta', uncert=True)[1]",
				"absorption_line('Fe5015', uncert=True)[1]",
				"absorption_line('Mg_b', uncert=True)[1]"
				]

	if D is None:
		pickleFile = open("%s/Data/%s/analysis/%s/%s/pickled/dataObj.pkl" % (
			cc.base_dir, instrument, galaxy, opt), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()
	if D.norm_method != norm:
		D.norm_method = norm
		D.find_restFrame()

	output = '%s/Data/%s/analysis/lims.txt' % (cc.base_dir, instrument)

	lims = np.loadtxt(output, dtype=str)

	for p in plots:
		array = eval('D.' + p)
		symmetric = False
		positive = True
		if 'vel' in p:
			positive = False
			symmetric = True
		r = set_lims(array, symmetric=symmetric, positive=positive)

		try:
			row = np.where(lims[:,0]==p)[0][0] # attribute
			lims[row,1] = str(np.nanmin([round(r[0],3), float(lims[row,1])])) 
			lims[row,2] = str(np.nanmax([round(r[1],3), float(lims[row,2])]))
			if lims[row,1] == str(round(r[0],3)):
				lims[row,3] = galaxy # extreme galaxy
		except:
			lims = np.append(lims, [[p, r[0], r[1], galaxy]], axis=0)

	np.savetxt(output, lims, fmt='%s %s %s %s')


	
if __name__=='__main__':
	if 'home' in cc.device:
		instrument = 'vimos'
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

			plots = [
				'flux',
				"components['stellar'].plot['vel']",
				"components['stellar'].plot['sigma']",
				"components['stellar'].plot['vel'].uncert",
				"components['stellar'].plot['sigma'].uncert"
			]

			find_limits(galaxy, plots=plots, opt='kin', 
				instrument=instrument)

			plots = [
				"absorption_line('G4300')",
				"absorption_line('Fe4383')",
				"absorption_line('Ca4455')",
				"absorption_line('Fe4531')",
				"absorption_line('H_beta')",
				"absorption_line('Fe5015')",
				"absorption_line('Mg_b')",
				"absorption_line('G4300', uncert=True)[1]",
				"absorption_line('Fe4383', uncert=True)[1]",
				"absorption_line('Ca4455', uncert=True)[1]",
				"absorption_line('Fe4531', uncert=True)[1]",
				"absorption_line('H_beta', uncert=True)[1]",
				"absorption_line('Fe5015', uncert=True)[1]",
				"absorption_line('Mg_b', uncert=True)[1]"
			]
			find_limits(galaxy, plots=plots, opt='pop', 
				instrument=instrument)

	elif cc.device == 'uni':
		instrument = 'muse'
		for galaxy in [
				'ic1531', 
				'ic4296',
				'ngc1316',
				'ngc1399']:

			plots = [
				'flux',
				"components['stellar'].plot['vel']",
				"components['stellar'].plot['sigma']",
				"components['stellar'].plot['vel'].uncert",
				"components['stellar'].plot['sigma'].uncert",
			]

			find_limits(galaxy, plots=plots, opt='kin', 
				instrument=instrument)

			plots = [
				"absorption_line('H_beta')",
				"absorption_line('Fe5015')",
				"absorption_line('Mg_b')",
				"absorption_line('Fe5270')",
				"absorption_line('Fe5335')",
				"absorption_line('Fe5406')",
				"absorption_line('Fe5709')",
				"absorption_line('NaD')",
				"absorption_line('Ti01')",
				"absorption_line('Ti02')",
				"absorption_line('H_beta', uncert=True)[1]",
				"absorption_line('Fe5015', uncert=True)[1]",
				"absorption_line('Mg_b', uncert=True)[1]",
				"absorption_line('Fe5270', uncert=True)[1]",
				"absorption_line('Fe5335', uncert=True)[1]",
				"absorption_line('Fe5406', uncert=True)[1]",
				"absorption_line('Fe5709', uncert=True)[1]",
				"absorption_line('NaD', uncert=True)[1]",
				"absorption_line('Ti01', uncert=True)[1]",
				"absorption_line('Ti02', uncert=True)[1]"
			]
			find_limits(galaxy, plots=plots, opt='pop', 
				instrument=instrument)