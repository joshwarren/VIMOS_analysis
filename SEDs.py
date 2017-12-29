# Routine to pull SEDs from .vot files in /Data/SED/[galaxy]/.
# Also to plot and find spectral index.

import numpy as np
import os
from astropy.io.votable import parse
from checkcomp import checkcomp
cc=checkcomp()
import matplotlib.pyplot as plt
from prefig import Prefig 
Prefig()

c = 299792458 # speed of light in m/s

def convert_ref(ref):
	convert = {
	'IX/10A/cor_nvs': 'Otrupcek et al. 1990',
	'IX/10A/cor_ver': 'Otrupcek et al. 1990',
	'J/A+A/386/97/tablea1': 'Jackson et al. 2002',
	'J/A+A/537/A99/lqac2': 'Souchay el.al. 2012',#
	'J/A+A/544/A18/master': 'van Velzen et al. 2012',#
	'J/A+A/544/A18/matches': 'van Velzen et al. 2012',#
	'J/A+A/583/A75/lqac3': 'Souchay el.al. 2015',#
	'J/A+A/598/A78/table3': 'Intema et al. 2017',#
	'J/AJ/126/2562/table3': 'Fomalont et al. 2003',
	'J/AJ/89/53/table3': 'Sadler 1984',#
	'J/ApJ/731/L41/table1': 'Brown et al. 2011',#
	'J/ApJ/737/45/table1': 'Ofek et al. 2011',#
	'J/ApJ/784/159/sample': 'Graham and Tingay 2014',#
	# 'J/ApJ/818/182/galaxies': '',
	'J/MNRAS/251/330/table2': 'Wright et al. 1991', #
	'J/MNRAS/402/2403/at20gcat': 'Murphy et al. 2010', #
	'J/MNRAS/415/1597/tables1': 'Coutois et al. 2011',#
	'J/MNRAS/417/2651/catalog': 'Mahony et al. 2011',#
	'J/MNRAS/434/956/table2': 'Chhetri et al. 2013',#
	'J/MNRAS/436/2915/table2': 'Massardi et al. 2013',#
	'J/MNRAS/440/696/tableb1': 'Allison et al. 2014',#
	'J/MNRAS/455/3249/paco': 'Massardi et al. 2016',#
	'VII/258/vv10': 'Allison et al. 2014',
	'VIII/100/gleamegc': 'Wayth et al. 2015',#
	'VIII/15/pkscat90': 'Otrupcek et al. 1990',#
	'VIII/81B/sumss212': 'Murphy et al. 2007',#
	'VIII/85A/spectra': 'Vollmer et al. 2009',#
	'VIII/85A/waste': 'Vollmer et al. 2009',#
	}
	if hasattr(ref, '__iter__'):
		return np.array([convert[r] if r in convert.keys() 
			else r for r in ref])
	else:
		if ref in convert.keys():
			return convert[ref]
		else:
			return ref



def SED(galaxy):

	SEDfile = '%s/Data/SED/%s/SED.txt' % (cc.base_dir, galaxy)
	if not os.path.exists(SEDfile):
		tab = parse('%s/Data/SED/%s/vizier_votable.vot' % (cc.base_dir, galaxy)
			).get_first_table()
		freq = np.array(tab.array['sed_freq'])
		flux = np.array(tab.array['sed_flux'])
		flux_err = np.array(tab.array['sed_eflux'])
		ref = np.array(tab.array['_tabname'])


		with open(SEDfile, 'w') as f:
			f.write('freq      flux       flux_err     ref\n')
			f.write('GHz       Jy         Jy           -\n')

			for i in range(len(freq)):
				f.write(str(round(freq[i],3))+'  '+str(flux[i])
					+'  '+str(flux_err[i])+'  '+ref[i]+'\n')
	else:
		freq, flux, flux_err = np.loadtxt(SEDfile, unpack=True, usecols=(0,1,2), 
			skiprows=2)
		ref = np.loadtxt(SEDfile, usecols=(3,), skiprows=2, dtype=str)

	def m_to_GHz(m):
		return c/m/10**9
	
	flux_24 = np.nanmean(flux[(freq > m_to_GHz(26*10**-6)) * (
		freq < m_to_GHz(22*10**-6))]) # flux at 24 micro-meters
	flux_24_err = np.nanstd(flux[(freq > m_to_GHz(26*10**-6)) * (
		freq < m_to_GHz(22*10**-6))]) 

	flux_14 = np.nanmean(flux[(freq < 1.5) * (freq > 1.3)]) # flux at 1.4 GHz
	flux_14_err = np.nanstd(flux[(freq < 1.5) * (freq > 1.3)]) 

	q24 = np.log10(flux_24/flux_14)
	q24_err = np.sqrt((flux_24_err/flux_24)**2 + (flux_14_err/flux_14)**2)

	# Only use radio flux
	flux = flux[freq < 15]
	flux_err = flux_err[freq < 15]
	ref = ref[freq < 15]
	freq = freq[freq < 15]

	order = np.argsort(freq)
	flux=flux[order]
	flux_err=flux_err[order]
	ref=ref[order]
	freq=freq[order]


	m = np.isfinite(np.log10(flux)) #* np.isfinite(1/np.log10(flux_err))
	flux = flux[m]
	flux_err = flux_err[m]
	# flux_err[np.isnan(flux_err)] = 1
	freq = freq[m]
	ref	= ref[m]

	line_x = np.array([freq[0]*0.8, freq[-1]*1.2])

	# Plot SED
	fig, ax = plt.subplots()
	for r in np.unique(ref):
		ax.errorbar(freq[ref==r], flux[ref==r], yerr=flux_err[ref==r], fmt='o')

	if galaxy == 'ic1459':
		Prefig(size=(14,8))
		fig2, ax2 = plt.subplots()
		a = np.array([])
		b = np.array([])
		ref = convert_ref(ref)
		# for r in np.unique(ref):
		for r in ['Sadler 1984',
			'Otrupcek et al. 1990',
			'Brown et al. 2011',
			'Mahony et al. 2011',
			'Ofek et al. 2011',
			'van Velzen et al. 2012',
			'Chhetri et al. 2013',
			'Graham and Tingay 2014',
			'Wayth et al. 2015',
			'Massardi et al. 2016']:

			if any((freq[ref==r] < 15) * (flux[ref==r] > 0.1) 
				* (flux[ref==r] < 2)):
				ax2.errorbar(freq[ref==r], flux[ref==r], 
					yerr=flux_err[ref==r], fmt='o', ms=5, 
					label=r)
				a = np.append(a, freq[ref==r])
				b = np.append(b, flux[ref==r])
			# if any((freq[ref==r] < 15) * (flux[ref==r] > 0.1) 
			# 	* (flux[ref==r] < 2)) and not include(r):
			# 	ax2.errorbar(freq[ref==r], flux[ref==r], 
			# 		yerr=flux_err[ref==r], fmt='x', ms=5, 
			# 		label=r)
		ax2.set_ylim([0.1,2])
		ax2.set_xscale('log')
		ax2.set_xlabel('Freq (GHz)')
		ax2.set_yscale('log')
		ax2.set_ylabel('Flux Density (Jy)')

		a = a[(b>0.1) * (b<2)]
		b = b[(b>0.1) * (b<2)]

		# Approximate function
		def func(x, e,f,g):
			return e*x**(f*np.log10(g*x))

		lims = ax2.get_xlim()
		x = np.logspace(np.log10(lims[0]), np.log10(lims[1]), num=100)
		from scipy.optimize import curve_fit
		e,f,g = curve_fit(func, a, b)[0]
		ax2.plot(x, func(x, e,f,g), 'k')
		ax2.set_xlim(lims)

		ax2.tick_params(which='major', direction='in', length=10, 
				width=2)
		ax2.tick_params(which='minor', direction='in', length=6, width=2)
		box = ax2.get_position()
		ax2.set_position([box.x0, box.y0, box.width * 0.65, box.height])	
		ax2.set_position([box.x0, box.y0, box.width * 0.65, box.height])
		ax2.legend(bbox_to_anchor=(1, 1), loc='upper left')
		
		# https://stackoverflow.com/questions/21920233/
		# 	matplotlib-log-scale-tick-label-number-formatting
		from matplotlib import ticker
		for axis in [ax2.xaxis, ax2.yaxis]:
			axis.set_major_formatter(
				ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))

		# ax2.set_title('Radio spectrum for %s' % (galaxy))
		fig2.savefig('%s/Data/SED/%s/index_detailed.png' % (
			cc.base_dir, galaxy))
		fig2.savefig('%s/Documents/thesis/appendix/appendix2/%s.png' %
			(cc.home_dir, galaxy), bbox_inches='tight', dpi=240)

	params, cov = np.polyfit(np.log10(freq), np.log10(flux), 1, cov=True)
	alpha, alpha_err = params[0], np.sqrt(np.diag(cov))[0]
	fit = 10**(np.poly1d(params)(np.log10(line_x)))
	ax.plot(line_x, fit, c='b', label='All points')

	fig.text(0.13, 0.16, 'Spectral index (all points): %.2f +/- %.2f' % (params[0], 
		np.sqrt(np.diag(cov))[0]), color='b')

	m = np.isfinite(np.log10(flux)) * np.isfinite(1/np.log10(flux_err))
	flux = flux[m]
	flux_err = flux_err[m]
	# flux_err[np.isnan(flux_err)] = 1
	freq = freq[m]
	params, cov = np.polyfit(np.log10(freq), np.log10(flux), 1, 
		w=1/np.log10(flux_err), cov=True)
	fit = 10**(np.poly1d(params)(np.log10(line_x)))
	ax.plot(line_x, fit, c='k', label='Only points with errors')
	
	fig.text(0.13, 0.12, 'Spectral index (only points with errors): %.2f +/- %.2f' % (
		params[0], np.sqrt(np.diag(cov))[0]), color='k')

	ax.legend(facecolor='w')

	ax.set_xlim(line_x[0]/0.9, line_x[1]/1.1)
	ax.set_xscale('log')
	ax.set_xlabel('Freq (GHz)')
	ax.set_yscale('log')
	ax.set_ylabel('Flux Density (Jy)')

	ax.set_title('Radio spectrum for %s' % (galaxy))

	fig.savefig('%s/Data/SED/%s/index.png' % (cc.base_dir, galaxy))
	plt.close('all')

	return q24, q24_err, alpha, alpha_err


if __name__=='__main__':
	SED('ic1459')
	kldsjflskd
	galaxies = ['eso443-g024','ic1459','ic1531','ic4296','ngc0612','ngc1316', 
		'ngc1399','ngc3100','ngc3557','ngc7075','pks0718-34']
	q24 = np.array([])
	q24_err = np.array([])
	alpha = np.array([])
	alpha_err = np.array([])
	for galaxy in galaxies: 
		q, qe, a, ae = SED(galaxy)

		q24 = np.append(q24, q)
		q24_err = np.append(q24_err, qe)
		alpha = np.append(alpha, a)
		alpha_err = np.append(alpha_err, ae)


	with open('%s/Data/galaxies_properties2.txt' % (cc.base_dir), 'w') as f:
		f.write('Galaxy      q24     q24_err     spec_index    spec_index_err\n')
		for i in range(len(galaxies)):
			f.write('%s   %.3f     %.3f     %.3f     %.3f \n' % (galaxies[i], q24[i],
				q24_err[i], alpha[i], alpha_err[i]))
