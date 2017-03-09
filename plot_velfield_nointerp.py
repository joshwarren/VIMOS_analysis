## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20151015 Routine to plot velocity (or similar) fields using
## methods from Michele Cappellari's plot_velfield routine and his 
## voronoi binning routine (both contained in the voronoi binning package),
## with some of my own inventions.
# warrenj 20160209 Added show_bin_num keyword.
# warrenj 20160413 Added flux_type keyword
# warrenj 20160615 New routine to handle plot as just ax item so that it
#                  can be used in subplots.
# warrenj 20161020 Added to plot in RA and dec if header keyword is provided
## *************************** KEYWORDS ************************* ##
# x_pix                 (Int Array) x/y coord of pixel
# y_pix                 As above
# bin_num               (Int Array) of bin numbers in order of above
# xBar_pix              (Array) x/y center of each bin in bin number order in 
#                           pixel units
# yBar_pix              As xBar_pix
# vel                   (Array) Values to be plotted in bin number order?
# vmin          None    (Double) Limits on values of vel: values outside of 
#                           this will be set to vmin/vmax
# vmax          None    As vmin
# nodots        True	(Boolean) to show a dot a x_pix,y_pix for each bin
# label         None    (String) Colorbar label
# flux          None    (Array) contains binned flux for plotting isophotes. 
#                           (Not recommended - use flux_unbinned)
# flux_unbinned None    (Array) contains flux at every spaxel for plotting isophotes. 
# galaxy        None    (String) of galaxy name for printing on plot (top left)
# redshift      None    (String) of galaxy redshift for printing on plot (top left)
#                           and display secondary axis with length scales.
# nticks        4       (Int) Default number of ticks on colorbar
# ncolors       64      (Int) Number of color levels in colar chart
# title         None    (String) of title of the plot
# save          None    (String) of location to save the plot to
# show_bin_number False (Boolean) to show bin number at x_pix,y_pix for each bin
#                           This overrides nodots=True
# show_vel      False (Boolean) to show value of vel at x_pix,y_pix for each bin
#                           This overrides nodots=True
# flux_type     'mag'   'mag'   Plot isophots in magnitudes (log)
#                       else    Plot in flux (linear)
# ax            None    (matplotlib.axes.Axes) axes to create the plot on. New
#                           axes are created if this is not supplied.
# header 		None 	(Header dictionary) If supplied, plot will be given in 
#						terms of RA and Dec. Will be in arcsec if not.
# close         False   (Boolean) to close the figure.
## ************************************************************** ##




import numpy as np # for reading files
from scipy.spatial import distance
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
# from sauron_colormap2 import sauron2 as sauron
from sauron_colormap import sauron
import math
import matplotlib.pyplot as plt # used for plotting
import matplotlib as mpl
import os



def plot_velfield_nointerp(x_pix, y_pix, bin_num, xBar_pix, yBar_pix, vel, 
	vmin=None, vmax=None, nodots=True, colorbar=False, label=None, flux=None, 
	flux_unbinned=None, galaxy = None, redshift = None, nticks=4, 
	ncolors=64, title=None, save=None, show_bin_num=False, flux_type='mag',
	ax = None, close=False, show_vel=False, header=None, signal_noise=None, 
	signal_noise_target=None, **kwargs):

	kwg = {}
	kwg.update(kwargs)

	
	if len(vel) != max(bin_num)+1:
		print "Not enough bins provided to vel keyword"
		return

	if ax is None:
		fig, ax = plt.subplots(nrows=1,ncols=1)

	if title is not None:
		ax.set_title(title, fontdict={'fontsize':'small'})

	if vmin is None:
		vmin = np.nanmin(vel)

	if vmax is None:
		vmax = np.nanmax(vel)

	xBar_pix, yBar_pix, vel = map(np.ravel, [xBar_pix, yBar_pix, vel])
	levels = np.linspace(vmin, vmax, ncolors)

	# Plot in arcsecs
	if header is None:
		res = 0.67 #arcsec per pixel
		xBar = xBar_pix*res
		yBar = yBar_pix*res
		x = x_pix*res
		y = y_pix*res

		x -= max(x)/2
		y -= max(y)/2
		x_label = r'$\Delta$ RA (arcsec)'
		y_label = r'$\Delta Dec (arcsec)'

		# Tells add_CO method in plot_results coords were not supplied
		ax.RaDec = False

	# Plot in RA and Dec
	else:
		from astropy.coordinates import SkyCoord
		from astropy import units as u

		coords = SkyCoord(header['HIERARCH CCD1 ESO INS IFU RA'], 
			header['HIERARCH CCD1 ESO INS IFU DEC'], 
			unit=(u.deg, u.deg))

		res = header['CDELT1'] #arcsec per pixel
		xBar = (xBar_pix-header['CRPIX1']-1)*res/(60**2) + coords.ra.degree
		yBar = (yBar_pix-header['CRPIX2']-1)*res/(60**2) + coords.dec.degree
		x = (x_pix-header['CRPIX1']-1)*res/(60**2) + coords.ra.degree
		y = (y_pix-header['CRPIX2']-1)*res/(60**2) + coords.dec.degree
		
		x_label = "RA"
		y_label = "Dec"

		# Tells add_CO in plot_results.py that plot is in terms of RA and Dec
		ax.RaDec = True

	tick_formatter = ticker.ScalarFormatter(useOffset=False)
	ax.xaxis.set_major_formatter(tick_formatter)
	ax.yaxis.set_major_formatter(tick_formatter)

	bin_num = bin_num.astype(int)

	pixelSize = np.min(distance.pdist(np.column_stack([x, y])))
	xmin, xmax = np.min(x), np.max(x)
	ymin, ymax = np.min(y), np.max(y)
	nx = int(round((xmax - xmin)/pixelSize) + 1)
	ny = int(round((ymax - ymin)/pixelSize) + 1)
	img = np.full((nx, ny), np.nan)  # use nan for missing data
	j = np.round((x - xmin)/pixelSize).astype(int)
	k = np.round((y - ymin)/pixelSize).astype(int)
	# Clipping vel and creating 2D array of image (img)
	img[j, k] = vel[bin_num].clip(vmin, vmax)
	
	if 'cmap' not in kwargs:
		cmap=kwargs.get('cmap',sauron)
	else:
		cmap = kwg['cmap']
		if isinstance(cmap, str):
			cmap = plt.get_cmap(cmap)

	# Change to RGBA 
	pic = cmap((img-np.nanmin(img))/(np.nanmax(img)-np.nanmin(img)))
	if signal_noise is not None:
		if signal_noise_target is None:
			pic[j, k, 3] = (signal_noise[bin_num]/max(signal_noise))*0.5+0.5
		else:
			pic[j, k, 3] = ((signal_noise[bin_num] - signal_noise_target/2)/\
				(signal_noise_target/2)).clip(0.05,1)
	# Set bad pixels grey
	pic[np.isnan(img),:] = [0.5,0.5,0.5,1]
	ax.set_axis_bgcolor('grey')


	# cs = ax.imshow(np.rot90(img), interpolation='none', clim=[vmin, vmax],
	# 	cmap=cmap, extent=[xmin, xmax, ymin, ymax])

	cs = ax.imshow(np.rot90(pic), interpolation='none', extent=[xmin, xmax, ymin, ymax],
		clim=[vmin, vmax], cmap=cmap) # clim and cmap supplied for colorbar

	# RA increases right to left
	ax.invert_xaxis()


	ax.set_ylabel(y_label)
	ax.set_xlabel(x_label)
	ax.set_aspect('equal')
	ax.autoscale(tight=True)
	ax.minorticks_on()
	ax.tick_params(length=10, which='major')
	ax.tick_params(length=5, which='minor')


	xmin_sav, xmax_sav = ax.get_xlim()
	ymin_sav, ymax_sav = ax.get_ylim()
	xlim=np.array([xmin_sav, xmax_sav])
	ylim=np.array([ymin_sav, ymax_sav])

	if galaxy is not None:
		gal_name = plt.text(0.02,0.98, galaxy, color='black',
			verticalalignment='top',transform=ax.transAxes)
		ax.gal_name = gal_name
		if redshift is not None:
			gal_z = plt.text(0.02,0.93, "Redshift: " + str(round(redshift,3)), 
				color = 'black',verticalalignment='top',
				transform=ax.transAxes)
			ax.gal_z = gal_z


	if flux is not None:
		# 1 mag contours
		ax.contour(-2.5*np.log10(flux/np.max(flux)), levels=np.arange(20), colors='k',
			extent=[xmin, xmax, ymin, ymax]) 

	if flux_unbinned is not None:
		if flux_type == 'mag':
			contours = -2.5*np.log10(flux_unbinned/np.max(flux_unbinned))
			# 1 mag contours
			ax.contour(contours, levels=np.arange(20), colors='k', 
				extent=[xmin, xmax, ymin, ymax])

		else:
			ax.contour(np.reshape(x,np.shape(flux_unbinned)),
					   np.reshape(y,np.shape(flux_unbinned)), 
					   flux_unbinned, colors='k')

	if not nodots and not show_bin_num and not show_vel:
		ax.plot(xBar, yBar, '.k',
				markersize=kwargs.get("markersize", 3))

	if show_bin_num and not show_vel:
		for i in range(len(xBar)):
			ax.text(xBar[i], yBar[i], str(i), 
				color='grey', fontsize=5)
	if show_vel:
		for i in range(len(xBar)):
			if i%5==0:
				ax.text(xBar[i]+pixelSize/2, yBar[i]+pixelSize/2, str(vel[i]), 
					color='grey', fontsize=5)
		

	if redshift is not None:
		c = 299792 #km/s
		#H = 67.8 #(km/s)/Mpc # From Planck
		H = 70.0 # value used by Bolonga group.
		xlim = np.radians(xlim/(60.0*60.0)) * redshift*c/H
		ylim = np.radians(ylim/(60.0*60.0)) * redshift*c/H
		xmax = xlim[1]
		ymax = ylim[1]
		#xlim -= xmax/2
		#ylim -= ymax/2
		axis_label = "Distance (Mpc)"
		if max(xlim) < 1.0:
			xlim *= 1000
			ylim *= 1000
			axis_label = "Distance (kpc)"


		ax2 = ax.twinx()
		#ax2.set_aspect('equal')
		#ax2.autoscale(tight=True)
		ax2.minorticks_on()
		ax2.tick_params(length=10, which='major')
		ax2.tick_params(length=5, which='minor')
		ax2.set_ylim(ylim[0],ylim[1])
		ax2.set_ylabel(axis_label, rotation=270)

		ax3 = ax.twiny()
		#ax3.set_aspect('equal')
		#ax3.autoscale(tight=True)
		ax3.minorticks_on()
		ax3.tick_params(length=10, which='major')
		ax3.tick_params(length=5, which='minor')
		ax3.set_xlim(xlim[0],xlim[1])
		ax3.set_xlabel(axis_label)

		ax.ax2 = ax2
		ax.ax3 = ax3

	if colorbar:
		ticks = ticker.MaxNLocator(nbins=nticks)
		if hasattr(ax,'ax3'):
			cbar = plt.colorbar(cs, ax=[ax,ax2,ax3], ticks=ticks)
		else:
			cbar = plt.colorbar(cs, ax=ax, ticks=ticks)

		cbar.ax.tick_params(labelsize=6)
		
		if label:
			#cbar.set_label(label, rotation=270, fontsize='small')
			cbar.ax.text(4.0,0.5, label, rotation=270, fontsize=8,
				verticalalignment='center')
			
		ax.cax = cbar.ax
	
	if save is not None:
		ax.saveTo = save
		if not os.path.exists(os.path.dirname(save)):
			os.makedirs(os.path.dirname(save))  

		plt.savefig(save, bbox_inches="tight")

	if close:
		plt.close(fig)

	return ax