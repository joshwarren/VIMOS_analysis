## ==================================================================
## 		Plot the absorption indices
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from astropy.io import fits
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from plot_results import set_lims, add_
from checkcomp import checkcomp
cc = checkcomp()
from errors2 import get_dataCubeDirectory, apply_range, set_params, run_ppxf
from pop import get_absorption
from prefig import Prefig
Prefig()

def plot_absorption(galaxy, D=None, uncert=True, opt='pop', overplot={}, 
	gradient=True):
	# Find lines:
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
		#'Mg_1', 'Mg_2', 
		'Mg_b']
	limits = {#'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
		'H_beta':[1.0,2.9], 'Fe5015':[3.5,5.9], 'Mg_b':[3.1,4.7]}

	print 'Absorption lines'


	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots/absorption" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None and gradient != 'only':
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header
	if not gradient: f.close()

	data_file =  "%s/galaxies.txt" % (out_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]
	z_gals, x_cent_gals, y_cent_gals, SN_target_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,4,5,col), dtype='float,int,int,float')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	SN_target=SN_target_gals[i_gal]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	figs = {}
	axs = {}
	fig,ax =plt.subplots()
	figs['Mg_sigma'] = fig
	axs['Mg_sigma'] = ax
	for line in lines:
		fig,ax =plt.subplots()
		figs[line] = fig
		axs[line] = ax
	if gradient !='only':
		for i, line in enumerate(lines):
			print "    " + line

			if uncert:
				ab_line, ab_uncert = D.absorption_line(line, uncert=True)
			else:
				ab_line = D.absorption_line(line)

			abmin, abmax = set_lims(ab_line, positive=True)

			# if line in limits.keys():
			# 	abmin = limits[line][0]
			# 	abmax = limits[line][1]
			

			saveTo = '%s/%s.png' % (out_plots, line)
			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
				D.xBar, D.yBar, ab_line, header, vmin=abmin, 
				vmax=abmax, nodots=True, colorbar=True, 
				label='Index strength ('+r'$\AA$'+')', title=line, 
				cmap='gnuplot2', redshift=z, flux_unbinned=D.unbinned_flux, 
				signal_noise=D.SNRatio, signal_noise_target=SN_target, 
				center=center, save=saveTo)
			ax.saveTo = saveTo
			count = 0
			if overplot:
				for o, color in overplot.iteritems():
					count +=1
					add_(o, color, ax, galaxy, close = count==len(overplot))
			
			if uncert:
				abmin, abmax = set_lims(ab_uncert)

				ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, ab_uncert, 
					header, vmin=abmin, vmax=abmax, nodots=True, colorbar=True, 
					label='Index strength ('+r'$\AA$'+')', title=line, cmap='gnuplot2', 
					flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
					signal_noise_target=SN_target, center=center, 
					save='%s/%s_uncert.png' % (out_plots, line), close=True)

			if gradient:
				r = np.sqrt((D.xBar - center[0])**2 + (D.yBar - center[1])**2)
				if uncert:
					axs[line].errorbar(r, ab_line, yerr=ab_uncert, fmt='.', 
						c='k')
				else:
					axs[line].scatter(r, ab_line, marker = 'x', c='k')

				if line == 'Mg_b':
					x = np.log10(D.components['stellar'].plot['sigma'])
					xerr = D.components['stellar'].plot['sigma'].uncert/\
						D.components['stellar'].plot['sigma']/np.log(10)

					axs['Mg_sigma'].errorbar(x, ab_line, 
						xerr=xerr, yerr=ab_uncert, fmt='.', c='k')

					if np.isfinite(ab_line[0]):
						params, cov = np.polyfit(x,ab_line, 1, 
							w=np.sqrt(1/(xerr)**2 + 1/ab_uncert**2), cov=True)

						lims = np.array(axs['Mg_sigma'].get_xlim())
						axs['Mg_sigma'].plot(lims, np.poly1d(params)(lims), 
							'k')



	if gradient:
		print '    Gradients'
		index = np.zeros((40,40,2))
		for i in range(40):
			for j in range(40):
				index[i,j,:] = np.array([i,j]) - center

		step_size = 2
		annuli = np.arange(2, 26, step_size)

		mg = np.array([])
		e_mg = np.array([])
		sigma = np.array([])
		e_sigma = np.array([])

		for a in annuli:

			params = set_params(reps=10, opt='pop', gas=1, produce_plot=False)

			mask = (np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) < a) * (
				np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) > a - step_size)

			spec = np.nansum(f[0].data[:,mask], axis=1)
			noise = np.sqrt(np.nansum(f[1].data[:,mask]**2, axis=1))

			lam = np.arange(len(spec))*header['CDELT3'] + header['CRVAL3']
			spec, lam, cut = apply_range(spec, lam=lam, 
				set_range=params.set_range, return_cuts=True)
			lamRange = np.array([lam[0],lam[-1]])
			noise = noise[cut]

			pp = run_ppxf(galaxy, spec, noise, lamRange, header['CDELT3'], 
				params)

			absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos')

			for line in lines:
				axs[line].errorbar(a, absorp[line], yerr=uncert[line], c='r',
					fmt='.')
			axs['Mg_sigma'].errorbar(np.log10(pp.sol[0][1]), absorp['Mg_b'], 
				xerr=np.std(pp.MCstellar_kin[:,1])/pp.sol[0][1]/np.log(10), 
				yerr=uncert['Mg_b'], c='r', fmt='.')

			mg = np.append(mg, absorp['Mg_b'])
			e_mg = np.append(e_mg, uncert['Mg_b'])
			sigma = np.append(sigma, np.log10(pp.sol[0][1]))
			e_sigma = np.append(e_sigma, 
				np.std(pp.MCstellar_kin[:,1])/pp.sol[0][1]/np.log(10))

		for line in lines:
			axs[line].set_xlabel('Radius')
			axs[line].set_ylabel(r'Index \AA')
			figs[line].savefig('%s/%s_grad.png' %(out_plots, line))
			plt.close(figs[line])

		if np.isfinite(mg[0]):
			mask = np.isfinite(mg)
			params, cov = np.polyfit(sigma[mask], mg[mask], 1, 
				w=np.sqrt(1/(e_sigma)**2 + 1/e_mg**2)[mask], cov=True)

			lims = np.array(axs['Mg_sigma'].get_xlim())
			axs['Mg_sigma'].plot(lims, np.poly1d(params)(lims), 
				'r')
			grad = params[0]
			e_grad = np.sqrt(np.diag(cov))[0]
		else:
			grad = np.nan 
			e_grad = np.nan


		axs['Mg_sigma'].set_xlabel(r'log $\sigma$ [km s$^{-1}$]')
		axs['Mg_sigma'].set_ylabel(r'Mg$_b \, \AA$')
		if galaxy not in ['pks0718-34', 'ngc0612']:
			axs['Mg_sigma'].set_ylim([0.9*min(mg), 1.1*max(mg)])
		axs['Mg_sigma'].set_xlim(np.log10(200), np.log10(350))
		figs['Mg_sigma'].savefig('%s/Mg_sigma.png' % (out_plots))
		plt.close(figs['Mg_sigma'])


		grad_file = '%s/galaxies_resolved_mg_sigma.txt' % (out_dir)
		galaxy_gals, grad_gals, e_grad_gals = np.loadtxt(grad_file, 
			unpack=True, dtype=str)
		i_gal = np.where(galaxy_gals == galaxy)[0][0]

		grad_gals[i_gal] = str(round(grad, 3))
		e_grad_gals[i_gal] = str(round(e_grad, 3))

		with open(grad_file, 'w') as f:
			temp = "{0:12}{1:15}{2:15}\n"
			for i in range(len(galaxy_gals)):
				f.write(temp.format(galaxy_gals[i], grad_gals[i], 
					e_grad_gals[i]))

	return D







##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ngc3100', overplot={'CO':'c', 'radio':'r'}, gradient='only')

