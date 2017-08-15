## ==================================================================
## 		Plot stellar population maps
## ==================================================================
## warrenj 20170331 Routine to plot the stellar populations found by pop.py
## on Glamdring. 


import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from astropy.io import fits
from errors2 import get_dataCubeDirectory
# from plot_results import set_lims
from checkcomp import checkcomp
cc = checkcomp()
from prefig import Prefig
Prefig()#subplots=(3,2))

def plot_stellar_pop(galaxy, method='median', D=None, opt='pop'):
	print 'Plotting stellar population'

	if cc.device == 'glamdring': vin_dir = '%s/analysis/%s/%s/pop' % (cc.base_dir,
		galaxy, opt)
	else: vin_dir = '%s/Data/vimos/analysis/%s/%s/pop' % (cc.base_dir, galaxy, opt)

	# Load pickle file from pickler.py
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None:
		pickle_file = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (pickle_file), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header
	f.close()

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

	age = np.zeros(D.number_of_bins)
	met = np.zeros(D.number_of_bins)
	alp = np.zeros(D.number_of_bins)
	unc_age = np.zeros(D.number_of_bins)
	unc_met = np.zeros(D.number_of_bins)
	unc_alp = np.zeros(D.number_of_bins)

	f, ax_array = plt.subplots(3, 2)

	if method == 'median':
		for i in xrange(D.number_of_bins):
			ag, me, al = np.loadtxt('%s/%i.dat' % (vin_dir, i), unpack=True)
			
			age[i] = ag[0]
			unc_age[i] = ag[1]
			met[i] = me[0]
			unc_met[i] = me[1]
			alp[i] = al[0]
			unc_alp[i] = al[1]

		f.suptitle('%s median and standard deviation' %(galaxy.upper()))


		

	elif method == 'mostlikely':
		# from peakdetect import peakdetect

		age1 = np.zeros(D.number_of_bins)
		met1 = np.zeros(D.number_of_bins)
		alp1 = np.zeros(D.number_of_bins)

		age2 = np.zeros(D.number_of_bins)
		met2 = np.zeros(D.number_of_bins)
		alp2 = np.zeros(D.number_of_bins)
		for i in xrange(D.number_of_bins):
			ag, me, al = np.loadtxt('%s/distribution/%i.dat' % (vin_dir, i), 
				unpack=True)

			for plot, unc_plot, pop in zip([age,met,alp],[unc_age,unc_met,unc_alp],
				[ag,me,al]):
				hist = np.histogram(pop, bins=40)
				x = (hist[1][0:-1]+hist[1][1:])/2
				hist = hist[0]
				# peaks = np.array(peakdetect(hist, x_axis=x, lookahead=4)[0])
				# plot[i] = peaks[np.argmax(peaks[:,1]), 0]
				plot[i] = x[np.argmax(hist)]

				gt_fwhm = hist >= np.max(hist)/2
				unc_plot[i] = np.max(x[gt_fwhm]) - np.min(x[gt_fwhm])

			f.suptitle('%s mostlikely and fwhm' % (galaxy.upper()))


	ax_array[0,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
		vmin=0, vmax=15, title='Age', ax=ax_array[0,0], cmap='gnuplot2', 
		flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
		signal_noise_target=SN_target, center=center)

	ax_array[1,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, met, header, nodots=True, colorbar=True, 
		label='Metalicity [Z/H]', vmin=-2.25, vmax=0.67, title='Metalicity', 
		ax=ax_array[1,0], cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
		signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)

	ax_array[2,0] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, alp, header, nodots=True, colorbar=True, 
		label='Element Ratio [alpha/Fe]', 
		vmin=-0.3, vmax=0.5, title='Alpha Enhancement', ax=ax_array[2,0], 
		cmap='gnuplot2', flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
		signal_noise_target=SN_target, center=center)


	ax_array[0,1] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, unc_age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
		vmin=0, vmax=15, title='Age Uncertainty', ax=ax_array[0,1], 
		cmap='gnuplot2', flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
		signal_noise_target=SN_target, center=center)

	ax_array[1,1] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, unc_met, header, nodots=True, colorbar=True, label='Metalicity', 
		vmin=0, vmax=0.67+2.25, title='Metalicity Uncertainty [Z/H]', 
		ax=ax_array[1,1], cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
		signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)

	ax_array[2,1] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
		D.xBar, D.yBar, unc_alp, header, nodots=True, colorbar=True, 
		label='Element Ratio [alpha/Fe]', 
		vmin=0, vmax=0.5+0.3, title='Alpha Enhancement Uncertainty', 
		ax=ax_array[2,1], cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
		signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)


	# f.set_size_inches(8.5,3*1.8)

	print 'Saving plot'

	saveTo = "%s/population.pdf" % (out_plots)
	# f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[1,0].set_xlabel('')
	ax_array[1,1].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	ax_array[2,1].set_ylabel('')
	f.savefig(saveTo)#, bbox_inches="tight")

	return D








##############################################################################

# Use of plot_stellar_pop.py

if __name__ == '__main__':
	plot_stellar_pop('eso443-g024')

