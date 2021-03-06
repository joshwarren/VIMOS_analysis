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
from scipy import odr
from lts_linefit import lts_linefit as lts
from Bin import myArray
from Bin2 import Data


c = 299792.458 # speed of light in km/s
H = 70.0 # Mpc/(km/s)

# in cm^2
def get_surface_area(galaxy):
	z_gals = np.loadtxt('%s/Data/vimos/analysis/galaxies.txt' % (cc.base_dir),
		usecols=(1,), skiprows=1, unpack=True)
	gals = np.loadtxt('%s/Data/vimos/analysis/galaxies.txt' % (cc.base_dir),
		usecols=(0,), skiprows=1, unpack=True, dtype=str)
	z = z_gals[np.where(gals==galaxy)[0][0]]

	return 4 * np.pi * (z*c/H * 3.0856775)**2 # 1e48 cm^2






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
				lims[row,1] = str(np.nanmin([round(r[0],3), 
					float(lims[row,1])])) 
				lims[row,2] = str(np.nanmax([round(r[1],3), 
					float(lims[row,2])]))
				if lims[row,2] == str(round(r[1],3)):
					lims[row,3] = galaxy # extreme galaxy
			except:
				lims = np.append(lims, [[p, r[0], r[1], galaxy]], axis=0)

	np.savetxt(output, lims, fmt='%s %s %s %s')



def find_stellardec(): # Find ratio between continuum at 6563A and 4861A
	if 'home' in cc.device:
		raise ValueError('This routine is for MUSE data only')
	else:
		Prefig()
		fig, ax = plt.subplots()
		ax.set_ylabel(r'$C_\mathrm{6563\AA}$')
		ax.set_xlabel(r'$C^\mathrm{corr}_\mathrm{4861\AA}$')
		# ax.set_xlabel(r'$C_\mathrm{4861\AA}$ (Corrected for dust)')

		Ha_cont = []
		e_Ha_cont = []
		Hb_cont = []
		e_Hb_cont = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ngc1316', 'ic1459']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			area = get_surface_area(galaxy)

			if 'Hbeta' in D.components.keys() and 'Halpha' in \
				D.components.keys():

				Ha = D.components['Halpha'].flux
				Hb = D.components['Hbeta'].flux
				
				bd = Ha/Hb	
				e_bd = bd * np.sqrt((Ha.uncert / Ha)**2 + (Hb.uncert / Hb)**2)

				dust_grad = bd/2.86/(6563 - 4861) # Flux per A
				e_dust_grad = e_bd/2.86/(6563 - 4861)


				Hb_cont_gal = myArray([bin.continuum[np.argmin(np.abs(bin.lam/
					(1 + bin.components['stellar'].vel/c) - 4861))] 
					for bin in D.bin]) * area

				Hb_cont_gal.uncert = np.array([
					bin.continuum.uncert[np.argmin(np.abs(bin.lam/
					(1 + bin.components['stellar'].vel/c) - 4861))] 
					for bin in D.bin]) * area

				Hb_cont_gal.uncert = np.sqrt((e_dust_grad / dust_grad)**2 
					+ (Hb_cont_gal.uncert / Hb_cont_gal)**2)

				Hb_cont_gal *= dust_grad * (6583 - 4861)
				Hb_cont_gal.uncert *= Hb_cont_gal

				Ha_cont_gal = myArray([bin.continuum[np.argmin(np.abs(bin.lam/
					(1 + bin.components['stellar'].vel/c) - 6563))] 
					for bin in D.bin])  * area
				Ha_cont_gal.uncert = np.array([
					bin.continuum.uncert[np.argmin(np.abs(bin.lam/
					(1 + bin.components['stellar'].vel/c) - 6563))] 
					for bin in D.bin]) * area
				ax.errorbar(Hb_cont_gal, Ha_cont_gal, 
					yerr = Ha_cont_gal.uncert, xerr = Hb_cont_gal.uncert, 
					fmt='.', label=galaxy)

				Ha_cont.extend(Ha_cont_gal)
				e_Ha_cont.extend(Ha_cont_gal.uncert)
				Hb_cont.extend(Hb_cont_gal)
				e_Hb_cont.extend(Hb_cont_gal.uncert)

				if galaxy == 'ngc1316':
					xlim = ax.get_xlim()
					ylim = ax.get_ylim()
					axins = zoomed_inset_axes(ax, 1.e9/2.5/np.ptp(ylim), 
						loc=2)

				axins.errorbar(Hb_cont_gal, Ha_cont_gal, 
					yerr = Ha_cont_gal.uncert, xerr = Hb_cont_gal.uncert, 
					fmt='.')


					# m = ~np.isnan(NII_gal) * ~np.isnan(NI_gal)

					# pivot=np.nanmean(NI_gal[m])
					# p = lts(NI_gal[m], NII_gal[m], e_NI_gal[m], 
					# 	NII_gal.uncert[m], pivot=pivot)

					# print 'NII to NI (lts)'
					# print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], 
					# 	p.ab_err[1], p.ab[0] - p.ab[1]*pivot, 
					# 	np.sqrt(p.ab_err[0]**2 + (p.ab_err[1]*pivot)**2 - 
					# 		2*p.ab_cov[0,1]*pivot))

					# lims = np.array(ax.get_xlim())

					# ax.plot(lims, np.poly1d(p.ab[::-1])(lims) 
					# 	- p.ab[1]*pivot, ':')

					# axins.plot(lims, np.poly1d(p.ab[::-1])(lims) 
					# 	- p.ab[1]*pivot, ':')


		Ha_cont = np.array(Ha_cont)
		e_Ha_cont = np.array(e_Ha_cont)
		Hb_cont = np.array(Hb_cont)
		e_Hb_cont = np.array(e_Hb_cont)

		lims = np.array(ax.get_xlim())
		ax_ylim = ax.get_ylim()

		m = ~np.isnan(Ha_cont) * ~np.isnan(Hb_cont)

		# pivot=np.nanmean(Hb_cont[m])
		# p = lts(Hb_cont[m], Ha_cont[m], e_Hb_cont[m], e_Ha_cont[m], 
		# 	pivot=pivot)

		# print 'C_6563 to C_4861 (lts)'
		# print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], p.ab_err[1],
		# 	p.ab[0] - p.ab[1]*pivot, np.sqrt(p.ab_err[0]**2 
		# 	+ (p.ab_err[1]*pivot)**2 - 2*p.ab_cov[0,1]*pivot))
		# print 'Variance matrix:'
		# for a in p.ab_cov[::-1]:
		# 	print '%.4g    %.4g' % (a[1], a[0])

		# ax.plot(lims, np.poly1d(p.ab[::-1])(lims) 
		# 	- p.ab[1]*pivot, 'k')

		# axins.plot(lims, np.poly1d(p.ab[::-1])(lims) 
		# 	- p.ab[1]*pivot, 'k')


		data = odr.RealData(Hb_cont[m], Ha_cont[m], sx=e_Hb_cont[m], 
			sy=e_Ha_cont[m])
		myodr = odr.ODR(data, odr.unilinear, beta0=[1.,0.])
		output = myodr.run()

		ax.plot(np.linspace(lims[0], lims[1], 100), 
			np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
			'b')
		axins.plot(np.linspace(xlim[0], xlim[1], 100), 
			np.poly1d(output.beta)(np.linspace(xlim[0], xlim[1], 100)), 
			'b')

		print 'C_6563 to C_4861 (ODR)'
		print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (output.beta[0], 
			output.sd_beta[0], output.beta[1], output.sd_beta[1])

		print 'Variance matrix:'
		for a in output.cov_beta:#[::-1]:
			print '%.4g    %.4g' % (a[0], a[1])


		ax.legend(facecolor='w', loc=4)

		ax.set_xlim(lims)
		ax.set_ylim(ax_ylim)

		axins.set_xlim(xlim)
		axins.set_ylim(ylim)
		axins.yaxis.set_tick_params(labelright='on', labelleft='off')

		mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

		fig.savefig('%s/Data/muse/analysis/stellar_ratio.png' % (cc.base_dir),
			dpi=240)


# Find relationship between [NII]/Hb and [NI]/Ha
def find_ratio_dec(thesis=False, paper=False): 	
	Prefig(size=(8,7))
	y, e_y = {}, {}
	x, e_x = {}, {}
	x2, e_x2 = {}, {}
	y2, e_y2 = {}, {}
	for galaxy in ['ic1459', 'ngc1316']:
		D = Data(galaxy, instrument='muse', opt='pop')

		NI = D.components['[NI]d'].flux
		e_NI = NI.uncert
		NII = D.components['[NII]6583d'].flux
		e_NII = NII.uncert

		NII /= 1.34
		e_NII /= 1.34

		alp = D.components['Halpha'].flux
		e_alp = alp.uncert
		bet = D.components['Hbeta'].flux
		e_bet = bet.uncert

		m = ~np.isnan(NI) * ~np.isnan(NII)
		x[galaxy] = NI[m]/bet[m]
		e_x[galaxy] = x[galaxy] * np.sqrt((e_NI[m]/NI[m])**2 + (e_bet[m]/bet[m])**2)
		y[galaxy] = NII[m]/alp[m]
		e_y[galaxy] = y[galaxy] * np.sqrt((e_NII[m]/NII[m])**2 
			+ (e_alp[m]/alp[m])**2)

		m = np.isfinite(D.components['Halpha'].equiv_width) * np.isfinite(D.components['Hbeta'].equiv_width)
		x2[galaxy] = D.components['Hbeta'].equiv_width[m]
		e_x2[galaxy] = x2[galaxy].uncert[m]
		y2[galaxy] = D.components['Halpha'].equiv_width[m]
		e_y2[galaxy] = y2[galaxy].uncert[m]

	def f(ab, pivot, x):
		return ab[1]*(x - pivot) + ab[0]
		
	pivot = np.nanmean(np.append(*x.values()))
	p = lts(np.append(*x.values()), np.append(*y.values()), 
		np.append(*e_x.values()), np.append(*e_y.values()), pivot=pivot, plot=True,
		text=False, corr=False)
	fig3 = plt.gcf()
	if not thesis and not paper:
		fig3.savefig('%s/Data/muse/analysis/ratio.png' % (cc.base_dir),
			dpi=240)
	fig, ax = plt.subplots(2,1, sharex=True,  gridspec_kw={'height_ratios':[3,1]})
	print 'NII/Ha to NI/Hb (lts)'
	print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], p.ab_err[1],
		p.ab[0] - p.ab[1]*pivot, np.sqrt(p.ab_err[0]**2 
		+ (p.ab_err[1]*pivot)**2 - 2*p.ab_cov[0,1]*pivot))
	print 'Variance matrix:'
	for a in p.ab_cov[::-1]:
		print '%.4g    %.4g' % (a[1], a[0])
	print ''

	for g in x.iterkeys():
		predicted = f(p.ab, pivot, x[g])
		fRes = (y[g] - predicted)/y[g]
		# yerr = res/x[g] * np.sqrt((e_x[g]/x[g])**2 + (e_y[g]**2 
		# 	+ pivot**2*p.ab_err[1]**2 + p.ab_err[0]**2 + p.ab[1]**2*e_x[g]**2 
		# 	+ x[g]**2*p.ab_err[1]**2)/res**2)
		fRes_err = fRes * np.sqrt(
			(((x[g] - pivot)**2 * p.ab_err[1]**2 + p.ab[1]**2 * e_x[g]**2 + p.ab_err[0]**2)
			/ predicted**2) + y[g]**2 * e_y[g]**2)

		ax[0].errorbar(x[g], y[g], xerr=e_x[g], yerr=e_y[g], fmt='.', label=g.upper())
		ax[1].errorbar(x[g], fRes, fmt='.', xerr=e_x[g], yerr=fRes_err)

	xlim = ax[0].get_xlim()
	ylim = ax[0].get_ylim()
	ax[0].plot(xlim, f(p.ab, pivot, xlim), 'k')
	ax[0].plot(xlim, f(p.ab, pivot, xlim) + p.rms, 'k--')
	ax[0].plot(xlim, f(p.ab, pivot, xlim) - p.rms, 'k--')
	ax[0].plot(xlim, f(p.ab, pivot, xlim) + 2.6 * p.rms, 'k:')
	ax[0].plot(xlim, f(p.ab, pivot, xlim) - 2.6 * p.rms, 'k:')
	ax[0].set_ylim(ylim)
	ax[0].set_xlim(xlim)

	ax[1].axhline(0, c='k')
	rms = np.std(np.append(*y.values()) - f(p.ab, pivot, np.append(*x.values())), 
		ddof=2)
	ax[1].axhline(+rms, ls='--', c='k')
	ax[1].axhline(-rms, ls='--', c='k')
	ax[1].axhline(+2.6*rms, ls=':', c='k')
	ax[1].axhline(-2.6*rms, ls=':', c='k')
	ax[1].set_ylim(-3.5*np.nanmax(rms), 3.5*np.nanmax(rms))

	ax[0].legend(facecolor='w', loc=4)

	ax[1].set_xlabel(r'$\mathrm{\frac{[NI]\lambda\lambda5197,5200}{H\beta}}$')
	# ax[1].set_xlabel(r'$\mathrm{\frac{[NI]\lambda5200}{H\beta}}$')
	ax[0].set_ylabel(r'$\mathrm{\frac{[NII]\lambda6584}{H\alpha}}$')
	ax[1].set_ylabel('Fractional\nResiduals')

	fig.subplots_adjust(hspace=0)
	if thesis:
		fig.savefig('%s/Documents/thesis/chapter5/ratio_fit.png' % (cc.home_dir),
			dpi=240, bbox_inches='tight')
	elif paper:
		fig.savefig('%s/Documents/paper/ratio_fit.png' % (cc.home_dir),
			dpi=240, bbox_inches='tight')
	plt.close('all')


	pivot = np.nanmean(np.append(*x2.values()))
	p = lts(np.append(*x2.values()), np.append(*y2.values()), 
		np.append(*e_x2.values()), np.append(*e_y2.values()), pivot=pivot, plot=True,
		text=False, corr=False)
	fig4 = plt.gcf()
	if not thesis and not paper:
		fig4.savefig('%s/Data/muse/analysis/EWrelationship.png' % (cc.base_dir),
			dpi=240)

	fig2, ax2 = plt.subplots(2,1, sharex=True,  gridspec_kw={'height_ratios':[3,1]})
	print 'EW(Ha) to EW(Hb) (lts)'
	print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], p.ab_err[1],
		p.ab[0] - p.ab[1]*pivot, np.sqrt(p.ab_err[0]**2 
		+ (p.ab_err[1]*pivot)**2 - 2*p.ab_cov[0,1]*pivot))
	print 'Variance matrix:'
	for a in p.ab_cov[::-1]:
		print '%.4g    %.4g' % (a[1], a[0])

	for g in x2.iterkeys():
		# res = y2[g] - f(p.ab, pivot, x2[g])
		# yerr = res/x2[g] * np.sqrt((e_x2[g]/x2[g])**2 + (e_y2[g]**2 
		# 	+ pivot**2*p.ab_err[1]**2 + p.ab_err[0]**2 + p.ab[1]**2*e_x2[g]**2 
		# 	+ x2[g]**2*p.ab_err[1]**2)/res**2)

		predicted = f(p.ab, pivot, x2[g])
		fRes = (y2[g] - predicted)/y2[g]
		fRes_err = fRes * np.sqrt(
			(((x2[g] - pivot)**2 * p.ab_err[1]**2 + p.ab[1]**2 * e_x2[g]**2 + p.ab_err[0]**2)
			/ predicted**2) + y2[g]**2 * e_y2[g]**2)

		ax2[0].errorbar(x2[g], y2[g], xerr=e_x2[g], yerr=e_y2[g], fmt='.', 
			label=g.upper())
		ax2[1].errorbar(x2[g], fRes, fmt='.', xerr=e_x2[g], yerr=fRes_err)

	xlim = ax2[0].get_xlim()
	ylim = ax2[0].get_ylim()
	ax2[0].plot(xlim, f(p.ab, pivot, xlim), 'k')
	ax2[0].plot(xlim, f(p.ab, pivot, xlim) + p.rms, 'k--')
	ax2[0].plot(xlim, f(p.ab, pivot, xlim) - p.rms, 'k--')
	ax2[0].plot(xlim, f(p.ab, pivot, xlim) + 2.6 * p.rms, 'k:')
	ax2[0].plot(xlim, f(p.ab, pivot, xlim) - 2.6 * p.rms, 'k:')
	ax2[0].set_ylim(ylim)
	ax2[0].set_xlim(xlim)

	ax2[1].axhline(0, c='k')
	rms = np.std(np.append(*y2.values()) - f(p.ab, pivot, np.append(*x2.values())), 
		ddof=2)
	ax2[1].axhline(+rms, ls='--', c='k')
	ax2[1].axhline(-rms, ls='--', c='k')
	ax2[1].axhline(+2.6*rms, ls=':', c='k')
	ax2[1].axhline(-2.6*rms, ls=':', c='k')
	ax2[1].set_ylim(-3.5*np.nanmax(rms), 3.5*np.nanmax(rms))

	ax2[0].legend(facecolor='w', loc=4)

	ax2[1].set_xlabel(r'EW(H$\beta$)')
	ax2[0].set_ylabel(r'EW(H$\alpha$)')
	ax2[1].set_ylabel('Fractional\nResiduals')

	fig2.subplots_adjust(hspace=0)
	
	if thesis:
		fig2.savefig('%s/Documents/thesis/chapter5/EqW_fit.png' % (cc.home_dir),
			dpi=240, bbox_inches='tight')
	elif paper:
		fig2.savefig('%s/Documents/paper/EqW_fit.png' % (cc.home_dir),
			dpi=240, bbox_inches='tight')
	else:
		fig2.savefig('%s/Data/muse/analysis/EqW_fit.png' % (cc.base_dir),
			dpi=240, bbox_inches='tight')




def find_ndec(): # Find ratio between [NII] and [NI]
	if 'home' in cc.device:
		raise ValueError('This routine is for MUSE data only')
	else:
		Prefig()
		fig, ax = plt.subplots()
		ax.set_ylabel(r'[NII]$\lambda\lambda$6548, 6584')
		ax.set_xlabel(r'[NI]$\lambda\lambda$5197, 5200 (Corrected for dust)')

		NII = []
		e_NII = []
		NI = []
		e_NI = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ngc1316', 'ic1459']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			area = get_surface_area(galaxy)

			if 'Hbeta' in D.components.keys() and 'Halpha' in \
				D.components.keys():

				Ha = D.components['Halpha'].flux
				Hb = D.components['Hbeta'].flux
				
				bd = Ha/Hb	
				e_bd = bd * np.sqrt((Ha.uncert / Ha)**2 + (Hb.uncert / Hb)**2)

				dust_grad = bd/2.86/(6563 - 4861) # Flux per A
				e_dust_grad = e_bd/2.86/(6563 - 4861)

				if '[NII]6583d' in D.components.keys() and '[NI]d' in \
					D.components.keys():
					NI_gal = D.components['[NI]d'].flux * area

					e_NI_gal = NI_gal.uncert * area
					e_NI_gal = np.sqrt((e_dust_grad / dust_grad)**2 
						+ (NI_gal.uncert / NI_gal)**2)

					NI_gal *= dust_grad * (6583 - 5200)
					e_NI_gal *= NI_gal

					NII_gal = D.components['[NII]6583d'].flux * area
					NII_gal.uncert *= area
					ax.errorbar(NI_gal, NII_gal, yerr = NII_gal.uncert,
						xerr = e_NI_gal, fmt='.', label=galaxy)

					NII.extend(NII_gal)
					e_NII.extend(NII_gal.uncert)
					NI.extend(NI_gal)
					e_NI.extend(e_NI_gal)

					if galaxy == 'ngc1316':
						xlim = ax.get_xlim()
						ylim = ax.get_ylim()
						# xlim = set_lims(NI_gal, symmetric=False, 
						# 	positive=True, n_std=2.5)
						# ylim = set_lims(NII_gal, symmetric=False, 
						# 	positive=True)
						axins = zoomed_inset_axes(ax, 3.e10/3.3/np.ptp(ylim), 
							loc=2)

					axins.errorbar(NI_gal, NII_gal, yerr=NII_gal.uncert,
						xerr=e_NI_gal, fmt='.')


					print galaxy, np.sum(~np.isnan(NI_gal) * ~np.isnan(NII_gal)
						), '/', D.number_of_bins


					# m = ~np.isnan(NII_gal) * ~np.isnan(NI_gal)

					# pivot=np.nanmean(NI_gal[m])
					# p = lts(NI_gal[m], NII_gal[m], e_NI_gal[m], 
					# 	NII_gal.uncert[m], pivot=pivot)

					# print 'NII to NI (lts)'
					# print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], 
					# 	p.ab_err[1], p.ab[0] - p.ab[1]*pivot, 
					# 	np.sqrt(p.ab_err[0]**2 + (p.ab_err[1]*pivot)**2 - 
					# 		2*p.ab_cov[0,1]*pivot))

					# lims = np.array(ax.get_xlim())

					# ax.plot(lims, np.poly1d(p.ab[::-1])(lims) 
					# 	- p.ab[1]*pivot, ':')

					# axins.plot(lims, np.poly1d(p.ab[::-1])(lims) 
					# 	- p.ab[1]*pivot, ':')


		NII = np.array(NII)
		e_NII = np.array(e_NII)
		NI = np.array(NI)
		e_NI = np.array(e_NI)

		m = ~np.isnan(NII) * ~np.isnan(NI)
		# params, cov = np.polyfit(NI[m], NII[m], 1, cov=True)

		lims = np.array(ax.get_xlim())
		ax_ylim = ax.get_ylim()

		# pivot=np.nanmean(NI[m])
		# p = lts(NI[m], NII[m], e_NI[m], e_NII[m], 
		# 	pivot=pivot)

		# print 'NII to NI (lts)'
		# print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (p.ab[1], p.ab_err[1],
		# 	p.ab[0] - p.ab[1]*pivot, np.sqrt(p.ab_err[0]**2 
		# 	+ (p.ab_err[1]*pivot)**2 - 2*p.ab_cov[0,1]*pivot))
		# print 'Variance matrix:'
		# for a in p.ab_cov[::-1]:
		# 	print '%.4g    %.4g' % (a[1], a[0])

		# ax.plot(lims, np.poly1d(p.ab[::-1])(lims) 
		# 	- p.ab[1]*pivot, 'r')

		# axins.plot(lims, np.poly1d(p.ab[::-1])(lims) 
		# 	- p.ab[1]*pivot, 'r')



		data = odr.RealData(NI[m], NII[m], sx=e_NI[m], sy=e_NII[m])
		myodr = odr.ODR(data, odr.unilinear, beta0=[1.,0.])
		output = myodr.run()

		ax.plot(np.linspace(lims[0], lims[1], 100), 
			np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
			'b')
		axins.plot(np.linspace(xlim[0], xlim[1], 100), 
			np.poly1d(output.beta)(np.linspace(xlim[0], xlim[1], 100)), 
			'b')

		print 'NII to NI (odr linear)'
		print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g' % (output.beta[0], 
			output.sd_beta[0], output.beta[1], output.sd_beta[1])

		print 'Variance matrix:'
		for a in output.cov_beta:#[::-1]:
			print '%.4g    %.4g' % (a[0], a[1])



		# data = odr.RealData(NI[m], NII[m], sx=e_NI[m], sy=e_NII[m])
		# myodr = odr.ODR(data, odr.quadratic, beta0=[0.,1.,0.])
		# output = myodr.run()

		# ax.plot(np.linspace(lims[0], lims[1], 100), 
		# 	np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
		# 	'k')
		# axins.plot(np.linspace(xlim[0], xlim[1], 100), 
		# 	np.poly1d(output.beta)(np.linspace(xlim[0], xlim[1], 100)), 
		# 	'k')

		# print 'NII to NI (odr quadratic)'
		# print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g, c = %.4g+/-%.4g' % (
		# 	output.beta[0], output.sd_beta[0], output.beta[1], output.sd_beta[1], 
		#	output.beta[2], output.sd_beta[2])

		# print 'Variance matrix:'
		# for a in output.cov_beta:#[::-1]:
		# 	print '%.4g    %.4g    %.4g' % (a[0], a[1], a[2])



		# fig.text(0.14,0.84,r'[NII] = (%.3f $\pm$ %.3f) [NI] + (%.3f $\pm$ %.3f)'%(
		# 	output.beta[0], output.sd_beta[0],
		# 	output.beta[1], output.sd_beta[1]))
		ax.legend(facecolor='w', loc=4)

		ax.set_xlim(lims)
		ax.set_ylim(ax_ylim)

		# axins.set_xlim(xlim)
		axins.set_xlim([0, 14])
		axins.set_ylim(ylim)
		axins.yaxis.set_tick_params(labelright='on', labelleft='off')

		mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

		fig.savefig('%s/Data/muse/analysis/NII_NI_ratio.png' % (cc.base_dir),
			dpi=240)



def find_odec(): # Find ratio between [OIII] and [OI]
	if 'home' in cc.device:
		raise ValueError('This routine is for MUSE data only')
	else:
		Prefig()
		fig, ax = plt.subplots()
		ax.set_xlabel(r'[OIII]$\lambda\lambda$5007, 4959')
		ax.set_ylabel(r'[OI]$\lambda\lambda$6300, 6364')

		OIII = []
		OI = []
		e_OIII = []
		e_OI = []
		n_bins = []
		# for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		for galaxy in ['ngc1316', 'ic1459']:
			pickleFile = open("%s/Data/muse/analysis/%s/%s/pickled" % (
				cc.base_dir, galaxy, 'pop')+"/dataObj.pkl", 'rb')
			D = pickle.load(pickleFile)
			pickleFile.close()

			area = get_surface_area(galaxy) # find intrinsic flux

			if 'Hbeta' in D.components.keys() and 'Halpha' in \
				D.components.keys():
				Ha = D.components['Halpha'].flux
				Hb = D.components['Hbeta'].flux
				
				bd = Ha/Hb	
				e_bd = bd * np.sqrt((Ha.uncert / Ha)**2 + (Hb.uncert / Hb)**2)

				dust_grad = bd/2.86/(6563 - 4861) # Flux per A
				e_dust_grad = e_bd/2.86/(6563 - 4861)

				if '[OIII]5007d' in D.components.keys() and \
					'[OI]6300d' in D.components.keys():
					OIII_gal = D.components['[OIII]5007d'].flux * area
					OIII_gal.uncert *= area
					e_OIII_gal = np.sqrt((e_dust_grad / dust_grad)**2 
						+ (OIII_gal.uncert / OIII_gal)**2)

					OIII_gal *= dust_grad * (6300 - 5007)
					e_OIII_gal *= OIII_gal

					OI_gal = D.components['[OI]6300d'].flux * area
					OI_gal.uncert *= area

					ax.errorbar(OIII_gal, OI_gal, yerr=OI_gal.uncert, 
						xerr=e_OIII_gal, fmt='.', label=galaxy)
					

					if galaxy == 'ngc1316':
						xlim = ax.get_xlim()
						ylim = ax.get_ylim()
						# xlim = set_lims(OIII_gal, symmetric=False, positive=True)
						# ylim = set_lims(OI_gal, symmetric=False, positive=True)
						axins = zoomed_inset_axes(ax, 5.e9/2.5/np.ptp(ylim), 
							loc=2)
						fig2, ax2 = plt.subplots()
						ax2.errorbar(OIII_gal, OI_gal, yerr=OI_gal.uncert,
							xerr=e_OIII_gal, fmt='.', label=galaxy)
						fig2.savefig('ngc1316_OIII_OI.png')

					axins.errorbar(OIII_gal, OI_gal, yerr=OI_gal.uncert,
						xerr=e_OIII_gal, fmt='.')


					OIII.extend(OIII_gal)
					e_OIII.extend(e_OIII_gal)
					OI.extend(OI_gal)
					e_OI.extend(OI_gal.uncert)


					print galaxy, np.sum(~np.isnan(OIII_gal) * ~np.isnan(OI_gal)
						), '/', D.number_of_bins

					n_bins.append(D.number_of_bins)

		OIII = np.array(OIII)
		e_OIII = np.array(e_OIII)
		OI = np.array(OI)
		e_OI = np.array(e_OI)

		# axins.set_xlim(set_lims(OIII[n_bins[1]:], symmetric=False, positive=True))
		# axins.set_ylim(set_lims(OI[n_bins[1]:], symmetric=False, positive=True))

		m = ~np.isnan(OIII) * ~np.isnan(OI)

		data = odr.RealData(OIII[m], OI[m], sx=e_OIII[m], sy=e_OI[m])
		myodr = odr.ODR(data, odr.quadratic, beta0=[0.,1.,0.])
		output = myodr.run()

		print 'OI to OIII'
		print 'a = %.4g+/-%.4g, b = %.4g+/-%.4g, c = %.4g+/-%.4g' % (output.beta[0], 
			output.sd_beta[0], output.beta[1], output.sd_beta[1], output.beta[2],
			output.sd_beta[2])

		print 'Variance matrix:'
		for a in output.cov_beta:#[::-1]:
			print '%.4g    %.4g    %.4g' % (a[0], a[1], a[2])

		lims = np.array(ax.get_xlim())
		ax_ylim = ax.get_ylim()
		ax.plot(np.linspace(lims[0], lims[1], 100), 
			np.poly1d(output.beta)(np.linspace(lims[0], lims[1], 100)), 
			'k')
		axins.plot(np.linspace(xlim[0], xlim[1], 100), 
			np.poly1d(output.beta)(np.linspace(xlim[0], xlim[1], 100)), 
			'k')

		# fig.text(0.14,0.84, r'[OI] = (%.3e $\pm$ %.3e) [OIII]$^2$ +' % (
		# 	output.beta[0], output.sd_beta[0]
		# 	) + r' (%.3f $\pm$ %.3f) [OIII] + (%.3f $\pm$ %.3f)' % (
		# 	output.beta[1], output.sd_beta[1], 
		# 	output.beta[2], output.sd_beta[2]))
		ax.legend(facecolor='w', loc=4)

		ax.set_xlim(lims)
		ax.set_ylim(ax_ylim)

		axins.set_xlim(xlim)
		axins.set_ylim(ylim)

		axins.yaxis.set_tick_params(labelright='on', labelleft='off')


		mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

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
		# find_stellardec()
		# find_ndec()
		find_ratio_dec(save_values=False)
		# find_odec()
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

