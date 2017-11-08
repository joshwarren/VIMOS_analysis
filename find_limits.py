from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import matplotlib.pyplot as plt
from astropy.io import fits # reads fits files (is from astropy)
import numpy as np # for array handling
import cPickle as pickle
from plot_results import set_lims
from prefig import Prefig



def find_limits(galaxy, opt='kin', norm='fit_disk', D=None, 
	instrument='vimos', plots=None):

	if plots is None and instrument=='vimos':
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
				"absorption_line('G4300',uncert=True)[1]",
				"absorption_line('Fe4383',uncert=True)[1]",
				"absorption_line('Ca4455',uncert=True)[1]",
				"absorption_line('Fe4531',uncert=True)[1]",
				"absorption_line('H_beta',uncert=True)[1]",
				"absorption_line('Fe5015',uncert=True)[1]",
				"absorption_line('Mg_b',uncert=True)[1]"
				]
	if plots is None and instrument=='muse':
		if 'kin' in opt:
			plots = [
				'flux',
				"components['stellar'].plot['vel']",
				"components['stellar'].plot['sigma']",
				"components['stellar'].plot['vel'].uncert",
				"components['stellar'].plot['sigma'].uncert",
			]

		elif 'pop' in opt:
			plots = [
				"absorption_line('H_beta')",
				"absorption_line('Fe5015')",
				"absorption_line('Mg_b')",
				"absorption_line('Fe5270')",
				"absorption_line('Fe5335')",
				"absorption_line('Fe5406')",
				"absorption_line('Fe5709')",
				"absorption_line('Fe5782')",
				"absorption_line('NaD')",
				"absorption_line('TiO1')",
				"absorption_line('TiO2')",
				"absorption_line('H_beta',uncert=True)[1]",
				"absorption_line('Fe5015',uncert=True)[1]",
				"absorption_line('Mg_b',uncert=True)[1]",
				"absorption_line('Fe5270',uncert=True)[1]",
				"absorption_line('Fe5335',uncert=True)[1]",
				"absorption_line('Fe5406',uncert=True)[1]",
				"absorption_line('Fe5709',uncert=True)[1]",
				"absorption_line('Fe5782',uncert=True)[1]",
				"absorption_line('NaD',uncert=True)[1]",
				"absorption_line('TiO1',uncert=True)[1]",
				"absorption_line('TiO2',uncert=True)[1]"
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
		if not (p == "components['stellar'].plot['vel']" and \
			galaxy in ['ngc0612', 'pks0718-34']):
			
			array = eval('D.' + p)
			symmetric = False
			positive = True
			if 'vel' in p:
				positive = False
				symmetric = True
			if 'uncert' in p:
				positive = True
				symmetric = False
			print galaxy, p, symmetric, positive
			r = set_lims(array, symmetric=symmetric, positive=positive)

			try:
				row = np.where(lims[:,0]==p)[0][0] # attribute
				lims[row,1] = str(np.nanmin([round(r[0],3), float(lims[row,1])])) 
				lims[row,2] = str(np.nanmax([round(r[1],3), float(lims[row,2])]))
				if lims[row,2] == str(round(r[1],3)):
					lims[row,3] = galaxy # extreme galaxy
			except:
				lims = np.append(lims, [[p, r[0], r[1], galaxy]], axis=0)

	np.savetxt(output, lims, fmt='%s %s %s %s')





def find_ndec(): # Find ratio between [NII] and [NI]
	if 'home' in cc.device:
		raise ValueError('This routine is for MUSE data only')
	else:
		Prefig()
		fig, ax = plt.subplots()
		ax.set_ylabel('[NII] flux')
		ax.set_xlabel('[NI] flux')

		NII = []
		NI = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ic1459', 'ngc1316']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			if '[NII]6583d' in D.components.keys() and '[NI]d' in \
				D.components.keys():
				ax.scatter(D.components['[NI]d'].flux, 
					D.components['[NII]6583d'].flux, marker='.', 
					label=galaxy)

				NII.extend(D.components['[NII]6583d'].flux)
				NI.extend(D.components['[NI]d'].flux)


				print galaxy, np.sum(
					~np.isnan(D.components['[NII]6583d'].flux) * 
					~np.isnan(D.components['[NI]d'].flux)),'/',\
					D.number_of_bins

		NII = np.array(NII)
		NI = np.array(NI)

		m = ~np.isnan(NII) * ~np.isnan(NI)
		params, cov = np.polyfit(NI[m], NII[m], 1, cov=True)

		lims = np.array(ax.get_xlim())
		ax.plot(lims, np.poly1d(params)(lims), 'k')

		fig.text(0.14,0.84, r'Gadient Ratio: %.3f $\pm$ %.3f' % (
			params[0], np.sqrt(np.diag(cov))[0]))
		ax.legend(facecolor='w', loc=4)

		fig.savefig('%s/Data/muse/analysis/NII_NI_ratio.png' % (cc.base_dir))

		

























	
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
			print galaxy

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
				"absorption_line('G4300',uncert=True)[1]",
				"absorption_line('Fe4383',uncert=True)[1]",
				"absorption_line('Ca4455',uncert=True)[1]",
				"absorption_line('Fe4531',uncert=True)[1]",
				"absorption_line('H_beta',uncert=True)[1]",
				"absorption_line('Fe5015',uncert=True)[1]",
				"absorption_line('Mg_b',uncert=True)[1]"
			]
			find_limits(galaxy, plots=plots, opt='pop', 
				instrument=instrument)

	elif cc.device == 'uni':
		instrument = 'muse'
		find_ndec()
		# for galaxy in [
		# 		'ic1459', 
		# 		'ic4296',
		# 		'ngc1316',
		# 		'ngc1399']:

		# 	plots = [
		# 		'flux',
		# 		"components['stellar'].plot['vel']",
		# 		"components['stellar'].plot['sigma']",
		# 		"components['stellar'].plot['vel'].uncert",
		# 		"components['stellar'].plot['sigma'].uncert",
		# 	]

		# 	find_limits(galaxy, plots=plots, opt='kin', 
		# 		instrument=instrument)

		# 	plots = [
		# 		"absorption_line('H_beta')",
		# 		"absorption_line('Fe5015')",
		# 		"absorption_line('Mg_b')",
		# 		"absorption_line('Fe5270')",
		# 		"absorption_line('Fe5335')",
		# 		"absorption_line('Fe5406')",
		# 		"absorption_line('Fe5709')",
		# 		"absorption_line('Fe5782')",
		# 		"absorption_line('NaD')",
		# 		"absorption_line('TiO1')",
		# 		"absorption_line('TiO2')",
		# 		"absorption_line('H_beta',uncert=True)[1]",
		# 		"absorption_line('Fe5015',uncert=True)[1]",
		# 		"absorption_line('Mg_b',uncert=True)[1]",
		# 		"absorption_line('Fe5270',uncert=True)[1]",
		# 		"absorption_line('Fe5335',uncert=True)[1]",
		# 		"absorption_line('Fe5406',uncert=True)[1]",
		# 		"absorption_line('Fe5709',uncert=True)[1]",
		# 		"absorption_line('Fe5782',uncert=True)[1]",
		# 		"absorption_line('NaD',uncert=True)[1]",
		# 		"absorption_line('TiO1',uncert=True)[1]",
		# 		"absorption_line('TiO2',uncert=True)[1]"
		# 	]
		# 	find_limits(galaxy, plots=plots, opt='pop', 
		# 		instrument=instrument)

