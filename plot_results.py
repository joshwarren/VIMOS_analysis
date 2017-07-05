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
# overplot  {}	Dictionary with keys with plots to be overplotted on the map, such 
#				as CO or radio, with the associated values being the matplotlib colors
#				that the contours are to be plotted in.
# D 		None Option to pass in the Data object instead of loading it.
## ************************************************************** ##

import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp 
from plot_histogram import plot_histogram
from errors2 import get_dataCubeDirectory
import os
from sauron_colormap2 import sauron2 as sauron
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
ain_dir = '%s/Data/alma' % (cc.base_dir)
out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)

#-----------------------------------------------------------------------------
class mapping(object):
	def __init__(self):
		self.SNR = True
		self.image = True
		self.equivalent_width = True
		self.amp_noise = True
		self.kinematics = True
		self.plot_resid = True
		self.line_ratios = True
	@property
	def all(self):
		return all([self.image, self.equivalent_width, self.kinematics, self.line_ratios])
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def set_lims(v, positive=False, symmetric=False):
	if all(~np.isfinite(v)):
		return 0, 0

	v = v[np.isfinite(v)]

	for i in range(2):
		av = np.median(v)
		std = np.std(v)

		include = (v >= av - 3*std) * (v <= av + 3*std)
		v = v[include]

	vmin, vmax = min(v), max(v)

	if symmetric:
		vmax = np.mean([vmax, abs(vmin)])
		vmin = -vmax

	if positive:
		vmin = max(vmin, 0)
		vmax = max(vmax, 0)

	return vmin, vmax
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
class set_ax_y_object(object):
	def __init__(self):
		self.Hd = False
		self.NI = False
	def set_ax_y(self, plt_title):
		if "stellar" in plt_title:
			ax_y=0
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
			if not self.Hd or self.NI: # assessing if NI or Hd have been plotted before.
				ax_y = 8
				if 'Hdelta' in plt_title:
					self.Hd = True
				else: 
					self.NI = True
			else:
				if 'Hdelta' in plt_title:
					ax_y = 8 if self.Hd else 10
				else:
					ax_y = 8 if self.NI else 10		   
		return ax_y
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_R_e(ax, galaxy, discard=0):
	from classify import get_R_e
	from  matplotlib.patches import Ellipse
	R_e = get_R_e(galaxy)
	
	data_file =  "%s/galaxies.txt" % (vin_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(4,5), dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str, unpack=True)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent_pix = x_cent_gals[i_gal]
	y_cent_pix = y_cent_gals[i_gal]

	xlims = ax.get_xlim()
	ylims = ax.get_ylim()

	x_cent = xlims[0] + (xlims[1] - xlims[0])/(40-discard*2)*x_cent_pix
	y_cent = ylims[0] + (ylims[1] - ylims[0])/(40-discard*2)*y_cent_pix




	data_file =  "%s/galaxies2.txt" % (vin_dir)
	ellip_gals, pa_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(2,3), dtype=float)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str, unpack=True)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	ellip = ellip_gals[i_gal]
	pa = pa_gals[i_gal]

	if ax.RaDec:
		patch = Ellipse([x_cent, y_cent], R_e*(1-ellip)/60/60, R_e/60/60, angle=pa, 
			fill=False)
	else:
		patch = Ellipse([x_cent, y_cent], R_e*(1-ellip), R_e, angle=pa, fill=False)
	ax.add_patch(patch)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_(overplot, color, ax, galaxy, close=False):
	image_dir=getattr(get_dataCubeDirectory(galaxy), overplot)
	
	if os.path.exists(image_dir):
		f = fits.open(image_dir)[0]

		# ****** NB: NOTE THE -VE SIGN ON CDELT1 ******
		x = (np.arange(f.header['NAXIS1']) - f.header['CRPIX1']) *\
			-f.header['CDELT1'] + f.header['CRVAL1'] + image_dir.RAoffset/(60**2)
		y = (np.arange(f.header['NAXIS2'])-f.header['CRPIX2']) *\
			f.header['CDELT2'] + f.header['CRVAL2'] + image_dir.decoffset/(60**2)
	
		#remove random extra dimenisons.
		s = np.array(f.data.shape)
		if any(s==1):
			image = np.sum(f.data, axis=tuple(np.where(s==1)[0]))
		else:
			image = f.data
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()

		# Discard noise from outer parts of the galaxy -  for radio
		if overplot == 'radio':
			lim = np.nanmean(image) + np.nanstd(image)
			image[image < lim] = lim
		else:
			Warning('Are you sure the x axis has the correct sign?')
		image = np.log(image)

		# Plot
		cs = ax.contour(x, y, image, colors=color, linestyles='solid', linewidth=1)
		# cs = ax.contour(image, colors=color, linestyles='solid', linewidth=1)
		cs.collections[0].set_label(overplot)

		ax.set_xlim(xlim)
		ax.set_ylim(ylim)

		leg = ax.legend(facecolor='w')

		# Save
		if hasattr(ax, 'saveTo'):
			saveTo = os.path.dirname(ax.saveTo)+"/Overplot/" + \
				os.path.basename(ax.saveTo)
			if not os.path.exists(os.path.dirname(saveTo)):
				os.makedirs(os.path.dirname(saveTo))
			plt.savefig(saveTo, bbox_inches="tight")

		if close:
			plt.close()
		else:
			leg.remove()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# @profile
def plot_results(galaxy, discard=0, norm="lwv", plots=False, residual=False, overplot={}, 
	show_bin_num=False, D=None, mapping=None, opt='kin'):	

	s = set_ax_y_object()
	data_file =  "%s/galaxies.txt" % (vin_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]

	# different data types need to be read separetly
	z_gals, x_cent_gals, y_cent_gals, SN_target_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,4,5,col), dtype='float,int,int,float')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	SN_target=SN_target_gals[i_gal]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots" % (output)
	out_nointerp = "%s/notinterpolated" % (out_plots)
	vin_dir_gasMC = "%s/%s/%s/MC" % (vin_dir, galaxy, opt)
	out_pickle = '%s/pickled' % (output)

	cubeFile = fits.open(get_dataCubeDirectory(galaxy))
	header = cubeFile[0].header
	cubeFile.close()
# ------------== Reading pickle file and create plot  ===----------

	# Load pickle file from pickler.py
	if D is None:
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	# for bin in D.bin:
	# 	for l in bin.e_line.keys():
	# 		bin.components[l].__threshold__ = 3
	
	if D.norm_method != norm:
		D.norm_method = norm
		D.find_restFrame()

	# Create figure and array for axes
	n_rows = 2+2*len(D.e_components) + int(np.ceil(len(D.e_components)*
		(len(D.e_components)-1)/6.0))
	f = plt.figure()#frameon=False)
	ax_array = []


	if mapping.SNR or mapping is None:	
		saveTo = "%s/SNR.png" % (out_nointerp)
		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.SNRatio, 
			header, colorbar=True, nodots=True, title='SNR', save=saveTo, close=True, 
			flux_unbinned=D.unbinned_flux, center=center)
# ------------=============== Plot image ================----------
	if mapping.image or mapping is None:
		print "    Image"
		
		title = "Total Flux"
		CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"

		ax = f.add_subplot(111, aspect='equal')
		saveTo = "%s/total_image.png" % (out_nointerp)
		ax.saveTo = saveTo
		ax.figx, ax.figy = 0, 0

		fmin, fmax = set_lims(D.flux, positive=True)

		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.flux, header,
			vmin=fmin, vmax=fmax, nodots=True, show_bin_num=show_bin_num, colorbar=True, 
			label=CBLabel, title=title, cmap='gist_yarg', ax=ax, 
			flux_unbinned=D.unbinned_flux, center=center)

		if plots:
			plt.show()

		ax_array.append(ax)
		f.delaxes(ax)
		f.delaxes(ax.cax)
		if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
		if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
# ------------========= Plot intensity (& EW) ===========----------
	if mapping.equivalent_width or mapping is None:
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
			f_min, f_max = set_lims(D.e_line[c].flux, positive=True)

			saveTo = "%s/%s_flux_hist.png" % (out_plots, c)
			plot_histogram(D.e_line[c].flux, galaxy=galaxy.upper(), redshift=z,
				vmin=f_min,vmax=f_max, weights=D.n_spaxels_in_bin, title=fh_title,
				xaxis=fCBtitle, save=saveTo)
			
			ax_y = s.set_ax_y(c)

			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/%s_img.png" % (out_nointerp, c)
			ax.saveTo = saveTo
			ax.figx, ax.figy = 0, ax_y
			
			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].flux, header, vmin=f_min, vmax=f_max, colorbar=True, 
				nodots=True, label=fCBtitle, title=f_title, ax=ax,
				flux_unbinned=D.unbinned_flux, center=center, 
				signal_noise=D.e_line[c].amp_noise, signal_noise_target=5)
			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
				
			if plots: plt.show()

			# Uncertainy in flux
			fu_title = "%s Flux Uncertainty" % (c_title)
			f_uncert_min, f_uncert_max = set_lims(D.e_line[c].flux.uncert, 
				positive=True)

			saveTo = "%s/%s_img_uncert.png" % (out_nointerp, c)
			ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				D.e_line[c].flux.uncert, header, vmin=f_uncert_min, vmax=f_uncert_max, 
				flux_unbinned=D.unbinned_flux, nodots=True, 
				show_bin_num=show_bin_num, colorbar=True, label=fCBtitle, 
				galaxy = galaxy.upper(), redshift = z, title=fu_title, 
				save=saveTo, close=True, center=center)

			

			eq_title = "%s Equivalent Width" % (c_title)
			eqh_title = "%s Equivalent Width Histogram" % (c_title)
			eqCBtitle = r"Equivalent Width ($\AA$)"

			eq_min, eq_max = set_lims(D.e_line[c].equiv_width, positive=True)

			saveTo = "%s/%s_eqWidth_hist.png" % (out_plots, c)
			plot_histogram(D.e_line[c].equiv_width, galaxy=galaxy.upper(), redshift=z,
				vmin=eq_min,vmax=eq_max, weights=D.n_spaxels_in_bin, title=eqh_title,
				xaxis=eqCBtitle, save=saveTo)
			
			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/%s_equiv_width.png" % (out_nointerp, c)
			ax.saveTo = saveTo
			ax.figx, ax.figy = 0, ax_y+1


			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].equiv_width, header, vmin=eq_min, vmax=eq_max, 
				colorbar=True, nodots=True, label=eqCBtitle, title=eq_title, ax=ax, 
				flux_unbinned=D.unbinned_flux, center=center, 
				signal_noise=D.e_line[c].amp_noise, signal_noise_target=5)
			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)

			# Uncertainy in EW
			equ_title = "%s Equivalent Width Uncertainty" % (c_title)
			eq_uncert_min, eq_uncert_max = set_lims(D.e_line[c].equiv_width.uncert, 
				positive=True)

			saveTo = "%s/%s_equiv_width_uncert.png" % (out_nointerp, c)
			ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				D.e_line[c].equiv_width.uncert, header, vmin=eq_uncert_min, 
				vmax=eq_uncert_max, flux_unbinned=D.unbinned_flux, nodots=True, 
				show_bin_num=show_bin_num, colorbar=True, label=eqCBtitle, 
				galaxy = galaxy.upper(), redshift = z, title=equ_title, 
				save=saveTo, close=True, center=center)
# ------------============ Amplitude/Noise ==============----------
	if mapping.amp_noise or mapping is None:
		for c in D.e_components:
			amp_title = '%s Amplitude to Noise ratio' % (c_title)
			amp_min, amp_max = set_lims(D.e_line[c].amp_noise, positive=True)
			saveTo = "%s/%s_amp_noise.png" % (out_nointerp, c)

			ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].amp_noise, header, vmin=amp_min, vmax=amp_max, 
				colorbar=True, nodots=True, title=amp_title, save=saveTo, 
				close=True, flux_unbinned=D.unbinned_flux, center=center)
# ------------=========== Setting titles etc ============----------
	if mapping.kinematics or mapping is None:
		print '    Kinematics'
		# for c in ['stellar']: # For debugging
		for c in D.independent_components:
			print '        %s' % (c)
			im_type = c
			pl = c

			if im_type == "gas":
				im_type=""
				pl = 'Hbeta'
			elif im_type == "SF":
				im_type=" (Star Forming)"
				pl = '[OIII]5007d'
			elif im_type == "Shocks":
				im_type=" (Shocking)"
				pl = 'Hbeta'
			elif 'Hbeta' in im_type:
				im_type=" ("+r'H$_\beta$'+")"
			elif 'Hgamma' in im_type:
				im_type=" ("+r'H$_\gamma$'+")"
			elif 'OIII' in im_type:
				im_type=" (OIII)"
			else:
				im_type=" (" + im_type + ")"


			for k in D.components[pl].plot.keys():
				ax_y = s.set_ax_y(c)

				symmetric=False
				positive=False
					
				CBLabel = None
				if k == "vel":
					ax_x=1
					title = 'Velocity'
					CBLabel = "V (km s$^{-1}$)"
					symmetric=True

				if  k == "sigma":
					ax_x=2
					title = 'Velocity Dispersion'
					CBLabel = r'$\mathrm{\sigma}$ (km s$^{-1}$)'
					positive = True

				if k == "h3":
					ax_x=1
					ax_y+=1
					title = 'h3'
					symmetric = True

				if k == "h4":
					ax_x=2
					ax_y+=1
					title = 'h4'

				SNR = D.SNRatio
				if pl != 'stellar':
					SNR = D.gas_dynamics_SN
					SN_target = 5


				if c == "stellar":
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
				vmin, vmax = set_lims(D.components[pl].plot[k][D.SNRatio > 
					0.75*SN_target], positive=positive, symmetric=symmetric)
				v_uncert_min, v_uncert_max = set_lims(D.components[pl].plot[k].uncert, 
					positive=True)
# # ------------============== Plot Histogram =============----------
			# Field histogram
			# saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title, wav_range)
			# plot_histogram(D.components[c].plot[k], galaxy=galaxy.upper(), redshift=z,
			# 	vmin=vmin,vmax=vmax, weights=D.n_spaxels_in_bin, title=htitle,
			# 	xaxis=CBLabel, save=saveTo)
			# # Uncertainty histogram
			# saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title+'_uncert', wav_range)
			# plot_histogram(D.components[c].plot[k].uncert, galaxy=galaxy.upper(), 
			# 	redshift=z, vmin=v_uncert_min,vmax=v_uncert_max, 
			# 	weights=D.n_spaxels_in_bin, title=uhtitle, xaxis=CBLabel, save=saveTo)

			# if plots:
			# 	plt.show()
# ------------==== Plot velfield - no interperlation ====----------
				# Field plot
				ax = f.add_subplot(111, aspect='equal')
				saveTo = ("%s/%s_%s_field.png" % (out_nointerp, c, k))
				ax.saveTo = saveTo
				ax.figx, ax.figy = ax_x, ax_y
				ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
					D.yBar, D.components[pl].plot[k], header, vmin=vmin, vmax=vmax, 
					# flux_type='notmag',
					nodots=True, show_bin_num=show_bin_num, colorbar=True, 
					label=CBLabel,galaxy = galaxy.upper(), redshift = z,
					title=title, ax=ax, signal_noise=D.SNRatio,
					signal_noise_target=SN_target, flux_unbinned=D.unbinned_flux, 
					center=center)
				# plots=True
				if plots:
					plt.show()
				ax_array.append(ax)
				f.delaxes(ax)
				f.delaxes(ax.cax)
				if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
				if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
				
				# Uncertainty plot
				saveTo = "%s/%s_%s_uncert_field.png" % (out_nointerp, c, k)
				ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
					D.components[pl].plot[k].uncert, header, vmin=v_uncert_min, 
					vmax=v_uncert_max, nodots=True, show_bin_num=show_bin_num,
					colorbar=True, label=CBLabel, galaxy = galaxy.upper(),
					redshift = z, title=utitle, save=saveTo, close=True, 
					flux_unbinned=D.unbinned_flux, center=center)
					
				#plots=False
				if plots:
					plt.show()
# ------------============= Plot residuals ==============----------
	if residual and (mapping.plot_resid or mapping is None):
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
				
		minres, maxres = set_lims(average_residuals, positive=True) #mean_centered=True,
		
		CBLabel = "Residuals"
		title = str.capitalize(residual) + \
		" Residuals of Bestfit to Normalised Spectrum"
		saveTo = "%s/%s_residual.png" % (out_nointerp, residual)

		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
			average_residuals, header, vmin=minres, vmax=maxres, flux_type='notmag',
			nodots=True, show_bin_num=show_bin_num, colorbar=True, 
			label=CBLabel, flux_unbinned=D.unbinned_flux, 
			galaxy = galaxy.upper(), redshift = z, title=title, 
			save=saveTo, close=True, center=center)
		if plots:
			plt.show()
# ------------============ Line ratio maps ==============----------
	if len(D.list_components)>2 and (mapping.line_ratios or mapping is None):
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

			lr_min, lr_max = set_lims(line_ratio)


			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/lineratio/%s_%s_line_ratio.png" % (out_nointerp, cB, cA)
			ax.saveTo = saveTo
			ax.figx, ax.figy = n%3, n_rows - int(np.ceil(t_num/3.0)) + int(np.ceil(n/3))

			ANRatio = np.min([D.e_line[cB].amp_noise, D.e_line[cA].amp_noise], axis=0)

			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				line_ratio, header, vmin=lr_min, vmax=lr_max, colorbar=True,
				nodots=True, title=lr_title, label=lrCBtitle, ax=ax,
				show_bin_num=show_bin_num, galaxy = galaxy.upper(), redshift = z,
				flux_unbinned=D.unbinned_flux, center=center, signal_noise=ANRatio,
				signal_noise_target=5)

			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
# ------------============= Plot and save ===============----------

	print "    Plotting and saving"

	cbar_position = []
	for i, a in enumerate(ax_array):
		f.add_axes(a)
		#a.axis('tight')
		f.add_axes(a.cax)
		if hasattr(a,'ax2'): f.add_axes(a.ax2)
		if hasattr(a,'ax3'): f.add_axes(a.ax3)
		if not os.path.exists(os.path.dirname(a.saveTo)):
			os.makedirs(os.path.dirname(a.saveTo))
		f.savefig(a.saveTo)#, bbox_inches="tight")

		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, a, galaxy)

			# Make lines thinner for pdf by finding the line objects
			for o in a.get_children():
				if type(o) is LineCollection:
					o.set_linewidth(0.3)

		cbar_position.append(a.cax.get_position())
		f.delaxes(a)
		f.delaxes(a.cax)
		if hasattr(a,'ax2'): f.delaxes(a.ax2)
		if hasattr(a,'ax3'): f.delaxes(a.ax3)
		a.change_geometry(n_rows, 3, a.figy*3+a.figx+1)
	if mapping.all or mapping is None:
		for i, a in enumerate(ax_array):
			if not np.isnan(a.figy):
				a2 = f.add_axes(a)
				# cax2 = f.add_axes(a.cax, position=cbar_position[i])
				p = a2.get_position()
				c = f.add_axes([p.x1,p.y0, 0.01, p.height])
				# c = f.add_axes([p.x1,p.y0*1.06-0.004,0.005, p.height*1.05])
				cb = plt.colorbar(a.images[0], cax=c)
				cb.outline.set_visible(False)
				cb.ax.tick_params(labelsize='x-small') 
				a.set_title(a.get_title(), fontdict={'fontsize':'small'})


					

				a.xaxis.set_visible(False)
				a.yaxis.set_visible(False)
				a.axis('off')
				# a.autoscale(False)
				# c.autoscale(False)
				if hasattr(a,'gal_name'): a.gal_name.remove()
				if hasattr(a, 'gal_z'): a.gal_z.remove()

		f.set_size_inches(8.5,n_rows*1.8)
		# f.tight_layout(h_pad=0.5)#pad=0.4, w_pad=0.5, h_pad=1.0)
		# f.subplots_adjust(top=0.94)
		f.suptitle(galaxy.upper())

		saveTo = "%s/grid.pdf" % (out_plots)
		f.savefig(saveTo)#, bbox_inches="tight",format='pdf')

	return D



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[5]

	discard = 0 # rows of pixels to discard- must have been the same 
			#	for all routines 
	print galaxy

	plot_results(galaxy, discard=discard, plots=False, residual = "median")