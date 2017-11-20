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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset



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
		e_NII = []
		NI = []
		e_NI = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ic1459', 'ngc1316']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			if 'Hbeta' in D.components.keys() and 'Halpha' in \
				D.components.keys():
				bd = D.components['Halpha'].flux/D.components['Hbeta'].flux
				e_bd = np.sqrt((D.components['Halpha'].flux.uncert/
					D.components['Halpha'].flux)**2 + 
					(D.components['Hbeta'].flux.uncert/
					D.components['Hbeta'].flux))
				dust_grad = bd/2.86/(6563 - 4861) # Flux per A
				e_dust_grad = e_bd/2.86/(6563 - 4861)

				if '[NII]6583d' in D.components.keys() and '[NI]d' in \
					D.components.keys():
					NI_gal = D.components['[NI]d'].flux * dust_grad*\
						(6563 - 5200)
					e_NI_gal = (6563 - 5200) * np.sqrt(e_dust_grad**2 
						+ (D.components['[NI]d'].flux.uncert
						/ D.components['[NI]d'].flux)**2)
					ax.errorbar(NI_gal, D.components['[NII]6583d'].flux, 
						yerr = D.components['[NII]6583d'].flux.uncert,
						xerr = e_NI_gal, fmt='.', label=galaxy)

					NII.extend(D.components['[NII]6583d'].flux)
					e_NII.extend(D.components['[NII]6583d'].flux.uncert)
					NI.extend(NI_gal)
					e_NI.extend(e_NI_gal)


					print galaxy, np.sum(
						~np.isnan(NI_gal) * 
						~np.isnan(D.components['[NII]6583d'].flux)),'/',\
						D.number_of_bins

		NII = np.array(NII)
		e_NII = np.array(e_NII)
		NI = np.array(NI)
		e_NI = np.array(e_NI)

		m = ~np.isnan(NII) * ~np.isnan(NI)
		# params, cov = np.polyfit(NI[m], NII[m], 1, cov=True)

		from scipy import odr
		data = odr.RealData(NI[m], NII[m], sx=e_NI[m], sy=e_NII[m])
		myodr = odr.ODR(data, odr.unilinear, beta0=[1.,0.])
		output = myodr.run()

		lims = np.array(ax.get_xlim())
		ax.plot(lims, np.poly1d(output.beta)(lims), 'k')

		fig.text(0.14,0.84, r'[NII] = (%.3f $\pm$ %.3f) [NI] + (%.3f $\pm$ %.3f)' % (
			output.beta[0], output.sd_beta[0],
			output.beta[1], output.sd_beta[1]))
		ax.legend(facecolor='w', loc=4)

		fig.savefig('%s/Data/muse/analysis/NII_NI_ratio.png' % (cc.base_dir))



def find_odec(): # Find ratio between [OIII] and [OI]
	if 'home' in cc.device:
		raise ValueError('This routine is for MUSE data only')
	else:
		Prefig()
		fig, ax = plt.subplots()
		ax.set_xlabel('[OIII] flux')
		ax.set_ylabel('[OI] flux')
		axins = zoomed_inset_axes(ax, 6, loc=1) # zoom = 6


		OIII = []
		OI = []
		e_OIII = []
		e_OI = []
		n_bins = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ic1459', 'ngc1316']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			if 'Hbeta' in D.components.keys() and 'Halpha' in \
				D.components.keys():
				bd = D.components['Halpha'].flux/\
					D.components['Hbeta'].flux
				e_bd = np.sqrt((D.components['Halpha'].flux.uncert/
					D.components['Halpha'].flux)**2 + 
					(D.components['Hbeta'].flux.uncert/
					D.components['Hbeta'].flux))
				dust_grad = bd/2.86/(6563 - 4861) # Flux per A
				e_dust_grad = e_bd/2.86/(6563 - 4861)

				if '[OIII]5007d' in D.components.keys() and \
					'[OI]6300d' in D.components.keys():
					OIII_gal = D.components['[OIII]5007d'].flux * \
						dust_grad * (6300 - 5007)
					e_OIII_gal = (6300 - 5007) * np.sqrt(e_dust_grad**2 
						+ (D.components['[OIII]5007d'].flux.uncert
						/ D.components['[OIII]5007d'].flux)**2)
					ax.errorbar(OIII_gal, D.components['[OI]6300d'].flux,
						yerr=D.components['[OI]6300d'].flux.uncert,
						xerr=e_OIII_gal, fmt='.', label=galaxy)
					axins.errorbar(OIII_gal, D.components['[OI]6300d'].flux,
						yerr=D.components['[OI]6300d'].flux.uncert,
						xerr=e_OIII_gal, fmt='.')

					OIII.extend(OIII_gal)
					e_OIII.extend(e_OIII_gal)
					OI.extend(D.components['[OI]6300d'].flux)
					e_OI.extend(D.components['[OI]6300d'].flux.uncert)


					print galaxy, np.sum(~np.isnan(OIII_gal) * 
						~np.isnan(D.components['[OI]6300d'].flux)),'/',\
						D.number_of_bins

					n_bins.append(D.number_of_bins)

		OIII = np.array(OIII)
		e_OIII = np.array(e_OIII)
		OI = np.array(OI)
		e_OI = np.array(e_OI)

		axins.set_xlim(set_lims(OIII[n_bins[1]:], symmetric=False, positive=True))
		axins.set_ylim(set_lims(OI[n_bins[1]:], symmetric=False, positive=True))

		mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

		m = ~np.isnan(OIII) * ~np.isnan(OI)
		# weighting is assumed to be for y data, so following line is 
		# not strictly true
		# params, cov = np.polyfit(OIII[m], OI[m], 2, cov=True,
		# 	w=1/np.sqrt((e_OI/OI)**2 + (e_OIII/OIII)**2))


		from scipy import odr
		data = odr.RealData(OIII[m], OI[m], sx=e_OIII[m], sy=e_OI[m])
		myodr = odr.ODR(data, odr.quadratic, beta0=[0.,1.,0.])
		output = myodr.run()

		lims = np.array(ax.get_xlim())
		ax.plot(np.linspace(lims[0], lims[1], 100), 
			np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
			'k')
		axins.plot(np.linspace(lims[0], lims[1], 100), 
			np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
			'k')

		fig.text(0.14,0.84, r'[OI] = (%.3e $\pm$ %.3e) [OIII]$^2$ + (%.3f $\pm$ %.3f) [OIII] + (%.3f $\pm$ %.3f)' % (
			output.beta[0], output.sd_beta[0], 
			output.beta[1], output.sd_beta[1], 
			output.beta[2], output.sd_beta[2]))
		ax.legend(facecolor='w', loc=4)

		fig.savefig('%s/Data/muse/analysis/OIII_OI_ratio.png' % (
			cc.base_dir))
	

























	
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
		# find_ndec()
		find_odec()
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

