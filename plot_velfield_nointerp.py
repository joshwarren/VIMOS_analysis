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
# nodots        False   (Boolean) to show a dot a x_pix,y_pix for each bin
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
# close         False   (Boolean) to close the figure.
## ************************************************************** ##




import numpy as np # for reading files
from scipy.spatial import distance
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from sauron_colormap import sauron
import math
import matplotlib.pyplot as plt # used for plotting
import matplotlib as mpl
import os



def plot_velfield_nointerp(x_pix, y_pix, bin_num, xBar_pix, yBar_pix, vel, 
	vmin=None, vmax=None, nodots=False, colorbar=False, label=None, flux=None, 
	flux_unbinned=None, galaxy = None, redshift = None, nticks=4, 
	ncolors=64, title=None, save=None, show_bin_num=False, flux_type='mag',
	ax = None, close=False, show_vel=False, **kwargs):

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

	res = 0.67 #arcsec per pixel
	xBar = xBar_pix*res
	yBar = yBar_pix*res
	x = x_pix*res
	y = y_pix*res

	x -= max(x)/2
	y -= max(y)/2
	axis_label = "Position (arcsec)"


	#im_xBar = np.copy(xBar)
	#im_yBar = np.copy(yBar)
	bin_num = bin_num.astype(int)

	v = vel[bin_num].clip(vmin,vmax)
	pixelSize = np.min(distance.pdist(np.column_stack([x, y])))

	xmin, xmax = np.min(x), np.max(x)
	ymin, ymax = np.min(y), np.max(y)
	nx = round((xmax - xmin)/pixelSize) + 1
	ny = round((ymax - ymin)/pixelSize) + 1
	img = np.full((nx, ny), np.nan)  # use nan for missing data
	j = np.round((x - xmin)/pixelSize).astype(int)
	k = np.round((y - ymin)/pixelSize).astype(int)
	img[j, k] = v
	
	if 'cmap' not in kwargs:
		cmap=kwargs.get('cmap',sauron)
	else:
		cmap = kwg['cmap']
		if isinstance(cmap, str):
			cmap = plt.get_cmap(cmap)
	#cmap.set_bad('darkslategray',1.0)
	cmap.set_bad('grey',1.0)

	cs = ax.imshow(np.rot90(img[:,:]), interpolation='none', 
		cmap=cmap,extent=[xmin - pixelSize/2, 
		xmax + pixelSize/2, ymin - pixelSize/2, ymax + pixelSize/2],
		clim = (vmin,vmax))

	# RA increases right to left
	ax.invert_xaxis()


	

	ax.set_ylabel(axis_label)
	ax.set_xlabel(axis_label)
	ax.set_aspect('equal')
	ax.autoscale(tight=True)
	ax.minorticks_on()
	ax.tick_params(length=10, which='major')
	ax.tick_params(length=5, which='minor')



	xmin_sav, xmax_sav = ax.get_xlim()
	ymin_sav, ymax_sav = ax.get_ylim()
	xlim=np.array([xmin_sav, xmax_sav])
	ylim=np.array([ymin_sav, ymax_sav])


   #y1, y2 = ax.get_ylim()   

	if galaxy is not None:
		gal_name = plt.text(0.02,0.98, "Galaxy: " + galaxy, color='black',
			verticalalignment='top',transform=ax.transAxes)
		ax.gal_name = gal_name
		if redshift is not None:
			gal_z = plt.text(0.02,0.93, "Redshift: " + str(round(redshift,3)), 
				color = 'black',verticalalignment='top',
				transform=ax.transAxes)
			ax.gal_z = gal_z


	if flux is not None:
		ax.tricontour(x, y, -2.5*np.log10(flux.ravel()/np.max(flux)),
					  levels=np.arange(20), colors='k') # 1 mag contours

	# NB: have assumed a square image!!!!
	if flux_unbinned is not None:
		#flux_unbinned[477]=0.001
		if flux_type == 'mag':
			contours = -2.5*np.log10(flux_unbinned.ravel()/
									 np.max(flux_unbinned))
			# 1 mag contours
			ax.tricontour(x, y, contours, levels=np.arange(20),
					  colors='k')

		else:
			ax.contour(np.reshape(x,np.shape(flux_unbinned)),
					   np.reshape(y,np.shape(flux_unbinned)), 
					   flux_unbinned, colors='k')

	if not nodots and not show_bin_num and not show_vel:
		#**********************************################
		#xBar /= res # no idea why this needs removing...
		#yBar /= res
		#**********************************################
		ax.plot(xBar-max(xBar)/2, yBar-max(yBar)/2, '.k',
				markersize=kwargs.get("markersize", 3))

	if show_bin_num and not show_vel:
		for i in range(len(xBar)):
			ax.text(xBar[i]-max(xBar)/2+pixelSize/2, 
				yBar[i]-max(yBar)/2+pixelSize/2, str(i), 
				color='grey', fontsize=5)
	if show_vel:
		for i in range(len(xBar)):
			if i%3==0:
				ax.text(xBar[i]-max(xBar)/2+pixelSize/2, 
					yBar[i]-max(yBar)/2+pixelSize/2, str(vel[i]), 
					color='grey', fontsize=5)
		

	#ax.axis('equal')

	#ax2 = ax.twinx()
	#ax3 = ax.twiny()
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
		#divider = make_axes_locatable(ax)
		#cax = divider.append_axes("right", size="5%", pad=0.05)
		#fig.subplots_adjust(right=0.66, top=0.84)
		#cax = fig.add_axes([0.75, 0.1, 0.02, 0.74])
		ticks = MaxNLocator(nbins=nticks)
		#cbar = plt.colorbar(cs, cax=cax, ticks=ticks)
		if hasattr(ax,'ax3'):
			cbar = plt.colorbar(cs, ax=[ax,ax2,ax3], ticks=ticks)
		else:
			cbar = plt.colorbar(cs, ax=ax, ticks=ticks)

		cbar.ax.tick_params(labelsize=6)

		# Manual method - documentation advices not to do it this way
		#cax = mpl.colorbar.make_axes(ax, location='right', fraction=0.05, pad=0.05, shrink=0.95)
		#norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
		#cbar = mpl.colorbar.ColorbarBase(cax[0], cmap=cmap, norm=norm, orientation='vertical', ticks=ticks)
		#cbar.ax.tick_params(labelsize=8)         


		
		if label:
			#cbar.set_label(label, rotation=270, fontsize='small')
			cbar.ax.text(4.0,0.5, label, rotation=270, fontsize=6,
				verticalalignment='center')
			
		ax.cax = cbar.ax
		#ax.cbar = cbar

	
	if save is not None:
		if not os.path.exists(os.path.dirname(save)):
			os.makedirs(os.path.dirname(save))  

		plt.savefig(save, bbox_inches="tight")

	if close:
		plt.close(fig)

	return ax




# if __name__ == '__main__':
#     pass