## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.
## warrenj 20151216 Added section to plot residuals.
## warrenj 20160111 Add section to plot histgrams of the fields.
## warrenj 20160405 Added keyword CO to overlay CO maps from ALMA if avaible.
## This supersedes plot_results_CO.py
## warrenj 20160726 Now plots in a more object orientated way and creates a
## grid of plots too. This supersedes plot_results2.py

## *************************** KEYWORDS ************************* ##
# galaxy 		Name of the galaxy being plotted: used to find 
#				correct files and to print onto the plot.
# discard	0	Interger giving the number of rows and columns 
#				to be removed from the plot to remove edge 
#				effects.
# wav_range 	null	Imposed wavelength range on top of the automated 
#				limits.	
# vLimit 	2 	Integer giving the number of lowest and highest 
#				results in the plot to be discarded. Defualt 
#				ignores 2 highest and 2 lowest bins.
# norm		"lwv"	Normalisation methods for velocity fields:
#				lwv: luminosity weighted mean of the whole 
#				field is set to 0.
#				lum: velocity of the brightest spaxel is set 
#				to 0.
#				sig: Noralised to the mean velocity of 5 bins with the
#				highest LOSVD.
# plots 	False   Boolean to show plots as routine runs.
# nointerp 	False 	Boolean to use interpolation between bins in 
#				plots or not.
# residual 	False	Method to measure the residuals:
#			mean: use the mean of the residuals in each 
#				bin.
#			median: use the median of the residuals in 
#				each bin.
#			max: use the maximum of the residuals in 
#				each bin.
#			False: do not calculate and produce plot of 
#				residuals.
# CO	   False	Boolean to show ALMA CO plots overlaied (if they exist)
# D 		None Option to pass in the Data object instead of loading it.
## ************************************************************** ##

#import matplotlib # 20160202 JP to stop lack-of X-windows error
#matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
#from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from plot_histogram import plot_histogram
#import ppxf_util as util
#from numpy.polynomial import legendre
import os
#import colormaps as cm
from sauron_colormap2 import sauron2 as sauron
#from Bin import Data
import cPickle as pickle
from checkcomp import checkcomp
cc = checkcomp()

# Give axes a saveTo property
plt.axes.saveTo = property(lambda self:str())
# Give axes an x and y on figure grid property
plt.axes.figx = property(lambda self:int())
plt.axes.figy = property(lambda self:int())
# give axes a property to hold a colorbar axes
plt.axes.cax = property(lambda self:plt.axes())
# give axes a property to hold 2 additional axes for showing other axis
plt.axes.ax2 = property(lambda self:plt.axes())
plt.axes.ax3 = property(lambda self:plt.axes())

vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
vin_dir_cube = '%s/Data/vimos/cubes' % (cc.base_dir)
ain_dir = '%s/Data/alma' % (cc.base_dir)
out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)

#-----------------------------------------------------------------------------
def set_lims(galaxy, v, vLimit, plot_species, plot_type, mean_centered=False,
	positive=False, symmetric=False):

	if 'uncert' in plot_type:
		uncert = True
		positive = True
	else:
		uncert = False

	if 'stellar' in plot_species:
		plot_species = 'stellar'
	elif 'Hbeta' in plot_species:
		plot_species = 'Hbeta'
	elif 'Hdelta' in plot_species:
		plot_species = 'Hdelta'
	elif 'Hgamma' in plot_species:
		plot_species = 'Hgamma'
	elif 'OIII' in plot_species:
		plot_species = 'OIII'

	if 'vel' in plot_type or 'velocity' in plot_type:
		plot_type = 'velocity'
		mean_centered = False
		symmetric = True
	elif 'sigma' in plot_type:
		plot_type = 'sigma'
	elif 'h3' in plot_type:
		plot_type = 'h3'
		symmetric = True
	elif 'h4' in plot_type:
		plot_type = 'h4'
	elif 'image' in plot_type or 'Flux' in plot_type:
		plot_type = 'image'
		postive = True
	elif 'Equivalent Width' in plot_type:
		plot_type = 'equivalent_width'
		positive = True
	if 'line ratio' in plot_type and 'OIII' in plot_type:
		plot_type = 'lrOIII'
		positive = True
	elif 'line ratio' in plot_type and 'Hbeta' in plot_type:
		plot_type = 'lrHbeta'
		positive = True

	# Read Limits file
	f = '%s/%s/limits.dat' % (vin_dir, galaxy)
	species, types = np.loadtxt(f, unpack=True, dtype=str, usecols=(0,1),
								skiprows=1)
	if uncert:
		c = (4,5)
	else:
		c = (2,3)
	mins, maxs = np.genfromtxt(f, unpack=True, usecols=c, skip_header=1, 
		missing_values='nan', filling_values=np.nan)

	# Is the plot in the limits file?
	if plot_species in species:
		s = np.where(species == plot_species)[0]
		if plot_type in types[s]:
			n = s[np.where(types[s] == plot_type)[0]]
			vmin, vmax = mins[n], maxs[n]
		else:
			vmin, vmax = np.nan, np.nan
	else: 
		vmin, vmax = np.nan, np.nan

	if np.isnan(vmin) or np.isnan(vmax):
		v_sorted = np.sort(np.unique(v[~np.isnan(v)]))
		# Is v_sorted too short for vLimit - normally means mostly nan.
		if len(v_sorted) < 2*vLimit:
			#v_sorted = np.sort(v[~np.isnan(v)])
			vmin = np.nanmin(v)
			vmax = np.nanmax(v)
		else:
			if np.isnan(vmin):
				vmin = v_sorted[vLimit]
			if np.isnan(vmax):
				vmax = v_sorted[-vLimit-1]
	
	# min must be less than max!
	if vmin > vmax:
		raise ValueError("min must be less than max limit.")
		vmin, vmax = vmax, vmin

	# center the range on the mean.
	if mean_centered:
		mean = np.mean(v)
		d = 2*mean - vmax
		vmin = mean-d

	# Uncertinty should always positive.
	if positive and vmin < 0:
		vmin = 0

	# Velocity plots should be symmetric about 0
	if symmetric:
		 if abs(vmin)>vmax:
			 # Check if both are < or > 0.
			 vmin=-abs(vmax)
		 else:
			 vmax=abs(vmin)

	return vmin, vmax
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def set_ax_y(plt_title):
	if "gas" in plt_title:
		ax_y=2
	elif "SF" in plt_title:
		ax_y=2
	elif "Shocks" in plt_title:
		ax_y=4
	elif 'Hbeta' in plt_title:
		ax_y=4
	elif 'Hgamma' in plt_title:
		ax_y=6
	elif 'OIII' in plt_title:
		ax_y=2
	elif 'stellar' in plt_title:
	   ax_y=0
	elif 'Hdelta' in plt_title or 'NI' in plt_title:
		ax_y = 8
	   
	return ax_y
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_CO(ax, galaxy, header, close=False):
	CO_image_dir="%s/%s-mom0.fits" % (ain_dir, galaxy)
	if os.path.exists(CO_image_dir):
		CO_image, CO_header = pyfits.getdata(CO_image_dir, 0, header=True)

		#remove random extra dimenisons.
		CO_image = np.sum(np.sum(CO_image,axis=0), axis=0)

		CO_x = np.arange(CO_header['NAXIS1'])*CO_header['CDELT1']*60*60
		CO_y = np.arange(CO_header['NAXIS2'])*CO_header['CDELT2']*60*60

		#x += max(ax.get_xlim())
		#y -= max(ax.get_ylim())

		# Coordinates of VIMOS pointing
		vc = SkyCoord(header['HIERARCH CCD1 ESO INS IFU RA'], 
			header['HIERARCH CCD1 ESO INS IFU DEC'], 
			unit=(u.deg, u.deg))

		# Coordinates of ALMA pointing
		ac = SkyCoord(CO_header['CRVAL1'], CO_header['CRVAL2'],
			unit=(u.deg, u.deg))

		# Offset between the two pointings
		CO_x -= ((vc.ra.degree - header['CRPIX1']*header['CDELT1']/(60*60)) -
			(ac.ra.degree +
			CO_header['CRPIX1']*CO_header['CDELT1']/(60*60)))*60*60
				
		CO_y += ((vc.dec.degree - header['CRPIX2']*header['CDELT2']/(60*60)) -
			(ac.dec.degree +
			CO_header['CRPIX2']*CO_header['CDELT2']/(60*60)))*60*60
			
		cs = ax.contour(CO_x,CO_y,CO_image, colors='k')

		saveTo = os.path.dirname(ax.saveTo)+"/withCO/" + \
			os.path.basename(ax.saveTo)
		if not os.path.exists(os.path.dirname(saveTo)):
			os.makedirs(os.path.dirname(saveTo))
		plt.savefig(saveTo, bbox_inches="tight")

		if close:
			plt.close()
		else:
			# Make lines thinner for pdf by finding the line objects
			for o in ax.get_children():
				if type(o) is LineCollection:
					o.set_linewidth(0.3)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def plot_results(galaxy, discard=0, wav_range="", vLimit=2, norm="lwv", 
	plots=False, nointerp=False, residual=False, CO=False, show_bin_num=False,
	D=None, **kwargs):	

	data_file =  "%s/galaxies.txt" % (vin_dir)
	# different data types need to be read separetly
	z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,4,5))
	#x_gals, y_gals = x_gals.astype(int), y_gals.astype(int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]


	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""

	dataCubeDirectory = "%s/%s.cube.combined.fits" % (vin_dir_cube, galaxy)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_plots = "%splots" % (output)
	out_nointerp = "%s/notinterpolated" % (out_plots)
	vin_dir_gasMC = "%s/%s/gas_MC" % (vin_dir, galaxy)
	out_pickle = '%s/pickled' % (output)

	# lists the files produced by man_errors[2].py
	outputs = glob.glob(output+'gal_*.dat')
	#outputs = glob.glob(output+'gal_stellar*.dat')
	#outputs = []

	# Used for CO plotting
	cubeFile = pyfits.open(dataCubeDirectory)
	header = cubeFile[0].header
	cubeFile.close()
# ------------== Reading pickle file and create plot  ===----------

	# Load pickle file from pickler.py
	if D is None:
		pickleFile = open("%s/dataObj_%s.pkl" % (out_pickle, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	# Create figure and array for axes
	n_rows = 2+2*len(D.e_components) + int(np.ceil(len(D.e_components)*
		(len(D.e_components)-1)/6.0))
	f = plt.figure(frameon=False)
	ax_array = []
# ------------=============== Plot image ================----------
	
	print "    Image"
	
	title = "Total Flux"
	CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"

	ax = f.add_subplot(111, aspect='equal')
	saveTo = "%s/total_image_%s.png" % (out_nointerp, wav_range)
	ax.saveTo = saveTo
	ax.figx, ax.figy = 0, 0

	fmin, fmax = set_lims(galaxy, D.flux, vLimit, 'stellar', title)

	ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.flux, vmin=fmin, 
		vmax=fmax, nodots=True, show_bin_num=show_bin_num, colorbar=True, 
		label=CBLabel, title=title, cmap='gist_yarg', ax=ax)#"gist_yarg", ax=ax)
	ax_array.append(ax)
	f.delaxes(ax)
	f.delaxes(ax.cax)
	if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
	if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
	
	if plots:
		plt.show()
# ------------========= Plot intensity (& EW) ===========----------
	print "    gas map(s) and equivalent widths"
	
	for c in D.e_components:
		print "        " + c

		if 'OIII' in c:
			c_title = '[OIII]'
		elif 'Hbeta' in c:
			c_title = r'H$_\beta$'
		elif 'Hgamma' in c:
			c_title = r'H$_\gamma$'
		else:
			c_title = c

		f_title = "%s Flux" % (c_title)
		fh_title = "%s Flux Histogram" % (c_title)
		# from header
		fCBtitle = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
		f_min, f_max = set_lims(galaxy, D.e_line[c].flux, vLimit, c, f_title)

		saveTo = "%s/%s_flux_hist_%s.png" % (out_plots, c, wav_range)
		plot_histogram(D.e_line[c].flux, galaxy=galaxy.upper(), redshift=z,
			vmin=f_min,vmax=f_max, weights=D.n_spaxels_in_bin, title=fh_title,
			xaxis=fCBtitle, save=saveTo)
		
		ax_y = set_ax_y(c)

		ax = f.add_subplot(111, aspect='equal')
		saveTo = "%s/%s_img_%s.png" % (out_nointerp, c, wav_range)
		ax.saveTo = saveTo
		ax.figx, ax.figy = 0, ax_y
		
		if 'OIII' in c:
			np.savetxt('flux.txt', D.e_line[c].mask)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.e_line[c].flux, 
			vmin=f_min, vmax=f_max, colorbar=True, nodots=True, label=fCBtitle, 
			  title=f_title, cmap = 'gist_yarg', ax=ax)
		ax_array.append(ax)
		f.delaxes(ax)
		f.delaxes(ax.cax)
		if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
		if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
			
		if plots: plt.show()
		

		eq_title = "%s Equivalent Width" % (c_title)
		eqh_title = "%s Equivalent Width Histogram" % (c_title)
		eqCBtitle = r"Equivalent Width ($\AA$)"

		eq_min, eq_max = set_lims(galaxy, D.e_line[c].equiv_width, vLimit, c, 
			eq_title)

		saveTo = "%s/%s_eqWidth_hist_%s.png" % (out_plots, c, wav_range)
		plot_histogram(D.e_line[c].equiv_width, galaxy=galaxy.upper(), redshift=z,
			vmin=eq_min,vmax=eq_max, weights=D.n_spaxels_in_bin, title=eqh_title,
			xaxis=eqCBtitle, save=saveTo)
		
		ax = f.add_subplot(111, aspect='equal')
		saveTo = "%s/%s_equiv_width_%s.png" % (out_nointerp, c, wav_range)
		ax.saveTo = saveTo
		ax.figx, ax.figy = 0, ax_y+1


		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
			D.e_line[c].equiv_width, vmin=eq_min, vmax=eq_max, colorbar=True, 
			nodots=True, label=eqCBtitle, title=eq_title, ax=ax)
		ax_array.append(ax)
		f.delaxes(ax)
		f.delaxes(ax.cax)
		if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
		if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
# ------------============ Amplitude/Noise ==============----------
		amp_title = '%s Amplitude to Noise ratio' % (c_title)
		amp_min, amp_max = set_lims(galaxy, D.e_line[c].amp_noise, vLimit, c, 
			amp_title)
		saveTo = "%s/%s_amp_nosie_%s.png" % (out_nointerp, c, wav_range)

		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
			D.e_line[c].amp_noise, vmin=amp_min, vmax=amp_max, colorbar=True, 
			nodots=True, title=amp_title, save=saveTo, close=not CO)
		if CO:
			ax1.saveTo = saveTo
			add_CO(ax1, galaxy, header, close=True)
# ------------=========== Setting titles etc ============----------
	for plot in outputs:
		plot_title = plot.split('gal_')[-1].split('.')[0]
		print "    " + plot_title
		[c, k] = plot_title.split('_')
		im_type = c
		if im_type == "gas":
			im_type=""
			ax_y=set_ax_y(plot_title)
		elif im_type == "SF":
			im_type=" (Star Forming)"
			ax_y=set_ax_y(plot_title)
		elif im_type == "Shocks":
			im_type=" (Shocking)"
			ax_y=set_ax_y(plot_title)
		elif 'Hbeta' in im_type:
			im_type=" ("+r'H$_\beta$'+")"
			ax_y=set_ax_y(plot_title)
		elif 'Hgamma' in im_type:
			im_type=" ("+r'H$_\gamma$'+")"
			ax_y=set_ax_y(plot_title)
		elif 'OIII' in im_type:
			im_type=" (OIII)"
			ax_y=set_ax_y(plot_title)
		else:
			im_type=" (" + im_type + ")"
			ax_y=set_ax_y(plot_title)

			
		CBLabel = None
		if "vel" in plot_title:
			ax_x=1
			title = 'Velocity'
			CBLabel = "V (km s$^{-1}$)"
			cmap = sauron

		else:
			cmap = sauron#cm.blue
		if "sigma" in plot_title:
			ax_x=2
			title = 'Velocity Dispersion'
			CBLabel = r'$\mathrm{\sigma}$ (km s$^{-1}$)'

		if "h3" in plot_title:
			ax_x=1
			ax_y+=1
			title = 'h3'

		if "h4" in plot_title:
			ax_x=2
			ax_y+=1
			title = 'h4'		



		if "stellar" in plot_title:
			utitle = "Stellar Uncertainty " + title + " Map"
			htitle = "Stellar " + title + " Histogram"
			uhtitle = "Stellar Uncertainty " + title + " Histogram"
			title = "Stellar " + title + " Map"
		else:
			utitle = "Ionised" + im_type + " Gas Uncertainty " + title + " Map"
			htitle = "Ionised" + im_type + " Gas " + title + " Histogram"
			uhtitle = "Ionised" + im_type + " Gas Uncertainty " + title + \
				" Histogram"
			title = "Ionised" + im_type + " Gas\n" + title + " Map"
# ------------============ Setting v range ==============----------
		vmin, vmax = set_lims(galaxy, D.components[c].plot[k], vLimit, plot_title, 
			plot_title)
		v_uncert_min, v_uncert_max = set_lims(galaxy, D.components[c].plot[k].uncert, 
			vLimit, plot_title, utitle)
# # ------------============== Plot Histogram =============----------
# 		# Field histogram
# 		saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title, wav_range)
# 		plot_histogram(D.components[c].plot[k], galaxy=galaxy.upper(), redshift=z,
# 			vmin=vmin,vmax=vmax, weights=D.n_spaxels_in_bin, title=htitle,
# 			xaxis=CBLabel, save=saveTo)
# 		# Uncertainty histogram
# 		saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title+'_uncert', wav_range)
# 		plot_histogram(D.components[c].plot[k].uncert, galaxy=galaxy.upper(), redshift=z,
# 			vmin=v_uncert_min,vmax=v_uncert_max, weights=D.n_spaxels_in_bin,
# 			title=uhtitle, xaxis=CBLabel, save=saveTo)

# 		if plots:
# 			plt.show()
# ------------==== Plot velfield - no interperlation ====----------
		if nointerp:
			# Field plot
			ax = f.add_subplot(111, aspect='equal')
			saveTo = ("%s/%s_field_%s.png" % (out_nointerp, plot_title, wav_range))
			ax.saveTo = saveTo
			ax.figx, ax.figy = ax_x, ax_y
			if 'vel' in k and 'OIII' in c:
				np.savetxt('vel.txt', D.components[c].plot[k]) ## save mask
			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
				D.yBar, D.components[c].plot[k], vmin=vmin, vmax=vmax, #flux_type='notmag',
				nodots=True, show_bin_num=show_bin_num, colorbar=True, 
				label=CBLabel, #flux_unbinned=D.unbinned_flux, 
				title=title, cmap=cmap, ax=ax)
			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)

			# Uncertainty plot
			saveTo = "%s/%s_field_%s.png" % (out_nointerp,plot_title+'_uncert', 
				wav_range)
			ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				D.components[c].plot[k].uncert, vmin=v_uncert_min, vmax=v_uncert_max,
				flux_type='notmag', nodots=True, show_bin_num=show_bin_num,
				colorbar=True, label=CBLabel, galaxy = galaxy.upper(),
				redshift = z, title=utitle, save=saveTo, close=not CO)#, cmap=cm.blue)
			if CO:
				ax1.saveTo = saveTo
				add_CO(ax1, galaxy, header, close=True)
			

		if plots:
			plt.show()
		#if CO:
		#	D.unbinned_flux = D.unbinned_flux_sav
# ------------============= Plot residuals ==============----------
	if residual:
		print "    " + residual + " residuals"

		average_residuals = np.zeros(D.number_of_bins)
		for i in range(D.number_of_bins):
			bestfit = np.loadtxt('%s/bestfit/%d.dat' % (vin_dir_gasMC, i))
			spectrum = np.loadtxt('%s/input/%d.dat' % (vin_dir_gasMC, i))
			residuals = np.abs(spectrum - bestfit)
			# remove edge pixels
			residuals = np.delete(residuals, [np.arange(5), 
			len(residuals)+np.arange(-5,0)], axis=0)

			if residual=="mean":
				average_residuals[i] = np.mean(residuals)
			elif residual=="median":
				average_residuals[i] = np.median(residuals)
			elif residual=="max":
				average_residuals[i] = np.max(np.abs(residuals))
				
		minres, maxres = set_lims(galaxy, average_residuals, vLimit, 'stellar', 
			'residuals', positive=True) #mean_centered=True,
		
		CBLabel = "Residuals"
		title = str.capitalize(residual) + \
		" Residuals of Bestfit to Normalised Spectrum"
		saveTo = "%s/%s_residual_%s.png" % (out_nointerp, residual, wav_range)

		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
			average_residuals, vmin=minres, vmax=maxres, flux_type='notmag',
			nodots=True, show_bin_num=show_bin_num, colorbar=True, 
			label=CBLabel, flux_unbinned=D.unbinned_flux, 
			galaxy = galaxy.upper(), redshift = z, title=title, 
			save=saveTo, close=not CO)#, cmap = cm.blue)
		if plots:
			plt.show()
		if CO:
			ax1.saveTo = saveTo
			add_CO(ax1, galaxy, header, close=True)
# ------------=============== Plot Chi2/DOF =============----------
	print "    chi2"

	chi2 = np.zeros(D.number_of_bins)
	for i in range(D.number_of_bins):
		chi2[i] = np.loadtxt("%s/chi2/%d.dat" % (vin_dir_gasMC, i))

	minchi2, maxchi2 = set_lims(galaxy, chi2, vLimit, 'stellar', 'chi2',
		positive = True)
	
	CBLabel = "Chi2/DOF"
	title = "Chi2/DOF of the bestfit"
	saveTo = "%s/chi2_%s.png" % (out_nointerp, wav_range)

	ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, chi2, 
		vmin=minchi2, vmax=maxchi2, flux_type='notmag',
		nodots=True, show_bin_num=show_bin_num, colorbar=True, 
		label=CBLabel, flux_unbinned=D.unbinned_flux, 
		galaxy = galaxy.upper(), redshift = z, title=title, 
		save=saveTo, close=not CO)#, cmap=cm.blue)
	if plots:
		plt.show()
	if CO:
		ax1.saveTo = saveTo
		add_CO(ax1, galaxy, header, close=True)
# ------------============ Line ratio maps ==============----------
	if any('OIII' in o for o in outputs):
		print "    line ratios"

		t_num = (len(D.e_components)-1)*len(D.e_components)/2
		for n in range(t_num):
			i = 0
			m = t_num
			while m > n:
				i += 1
				m -= i

			cA = D.e_components[len(D.e_components)-i-1]
			cB = D.e_components[len(D.e_components)-i+n-m]

			line_ratio = np.log10(D.e_line[cB].flux/D.e_line[cA].flux)
			if 'OIII' in cA:
				cA_title = '[OIII]'
			elif 'Hbeta' in cA:
				cA_title = r'H$_\beta$'
			elif 'Hdelta' in cA:
				cA_title = r'H$_\delta$'
			elif 'Hgamma' in cA:
				cA_title = r'H$_\gamma$'
			else:
				cA_title = cA

			if 'OIII' in cB:
				cB_title = '[OIII]'
			elif 'Hbeta' in cB:
				cB_title = r'H$_\beta$'
			elif 'Hdelta' in cB:
				cB_title = r'H$_\delta$'
			elif 'Hgamma' in cB:
				cB_title = r'H$_\gamma$'
			else:
				cB_title = cB
				
			lr_title = "%s/%s Line Ratio" % (cB_title, cA_title)
			lrCBtitle = r"log$_{10}$ (%s/%s)" %(cB_title,cA_title)

			lr_min, lr_max = set_lims(galaxy, line_ratio, vLimit, cA, lr_title)


			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/lineratio/%s_%s_line_ratio_%s.png" % (out_nointerp, cB, cA, 
				wav_range)
			ax.saveTo = saveTo
			ax.figx, ax.figy = n, n_rows-int(np.ceil(t_num/3))


			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				line_ratio, vmin=lr_min, vmax=lr_max, colorbar=True,
				nodots=True, title=lr_title, label=lrCBtitle, ax=ax)

			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)	
# ------------============= Plot and save ===============----------

	print "    Plotting and saving"

	for i, a in enumerate(ax_array):
		f.add_axes(a)
		#a.axis('tight')
		f.add_axes(a.cax)
		if hasattr(a,'ax2'): f.add_axes(a.ax2)
		if hasattr(a,'ax2'): f.add_axes(a.ax3)
		if not os.path.exists(os.path.dirname(a.saveTo)):
			os.makedirs(os.path.dirname(a.saveTo))  
		f.savefig(a.saveTo, bbox_inches="tight")

		if CO:
			add_CO(a, galaxy, header)

		f.delaxes(a)
		f.delaxes(a.cax)
		a.change_geometry(n_rows, 3, a.figy*3+a.figx+1)
		if hasattr(a,'ax2'): f.delaxes(a.ax2)		
		if hasattr(a,'ax2'): f.delaxes(a.ax3)



	for a in ax_array:
		if not np.isnan(a.figy):
			f.add_axes(a)
			f.add_axes(a.cax)
			a.xaxis.set_visible(False)
			a.yaxis.set_visible(False)
			a.axis('off')
			a.autoscale(False)

	f.set_size_inches(8.5,n_rows*1.8)
	f.tight_layout(h_pad=0.5)#pad=0.4, w_pad=0.5, h_pad=1.0)
	f.subplots_adjust(top=0.94)
	f.suptitle(galaxy.upper())

	saveTo = "%s/grid_%s.pdf" % (out_plots, wav_range)
	f.savefig(saveTo, bbox_inches="tight",format='pdf')

	return D



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[5]

	wav_range="4200-"
	discard = 2 # rows of pixels to discard- must have been the same 
			#	for all routines 
	vLimit = 2 #

	plot_results(galaxy, discard=discard, vLimit=vLimit, 
		wav_range=wav_range, plots=False, nointerp = True, residual = "median")