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
# flux_unbinned None    (Array) contains flux at every spaxel for plotting 
# 							isophotes. 
# galaxy        None    (String) of galaxy name for printing on plot (top left)
# redshift      None    (String) of galaxy redshift for printing on plot (top 
# 							left) and display secondary axis with length scales.
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
from prefig import Prefig 

## Adding the origentation arrows to a plot. 
## PA is given in degrees
## NB: currently will give arrow different lengths if image is not square.
def add_orientation(pa=0):
	pa = -np.radians(pa)
	ax=plt.gca()

	xlim = ax.get_xlim()
	xr = xlim[1]-xlim[0]
	ylim = ax.get_ylim()
	yr = ylim[1]-ylim[0]

	# East pointing arrow
	ax.arrow(0.8*xr+xlim[0], 0.8*yr+ylim[0], -0.1*xr*np.cos(pa), 
		0.1*yr*np.sin(pa), length_includes_head=True, head_width=0.02*xr, 
		head_length=0.03*xr, fc='k', ec='k')

	# North pointing arrow
	ax.arrow(0.8*xr+xlim[0], 0.8*yr+ylim[0], 0.1*xr*np.sin(pa), 
		0.1*yr*np.cos(pa), length_includes_head=True, head_width=0.02*xr, 
		head_length=0.03*xr, fc='k', ec='k')


	ax.text(0.77*xr+xlim[0]-0.1*xr*np.cos(pa), 
		0.78*yr+ylim[0]+0.1*yr*np.sin(pa), 'E')
	ax.text(0.78*xr+xlim[0]+0.1*xr*np.sin(pa), 
		0.805*yr+ylim[0]+0.1*yr*np.cos(pa), 'N')


def render_numbers(vmin, vmax):
    order = int(np.floor(np.log10(abs(vmax))))
    if order > 0: # vmax >= 10
        out = str(int(round(vmin, 0))), str(int(round(vmax, 0)))
    elif order == 0: # 1 >= vmax > 10
        out = '%.1f' % (vmin), '%.1f' % (vmax)
    else: # vmax < 1
        out = str(round(vmin, abs(order)+2)), str(round(vmax, abs(order)+2))
    return out


def plot_velfield_nointerp(x_pix, y_pix, bin_num, xBar_pix, yBar_pix, vel, 
	header, vmin=None, vmax=None, nodots=True, colorbar=False, label=None, 
	flux=None, flux_unbinned=None, galaxy = None, redshift = None, nticks=4, 
	ncolors=64, title=None, save=None, show_bin_num=False, flux_type='mag',
	ax = None, close=False, show_vel=False, signal_noise=None, debug=False, 
	dms=False, signal_noise_target=None, pa=None, center=None, alpha=None, 
	galaxy_labelcolor=None, label_size=None, lim_labels=False, 
	lim_labels_units=None, **kwargs):
	Prefig()

	kwg = {}
	kwg.update(kwargs)
	x_pix = np.array(x_pix)
	y_pix = np.array(y_pix)
	bin_num = np.array(bin_num)
	xBar_pix = np.array(xBar_pix)
	yBar_pix = np.array(yBar_pix)
	vel = np.array(vel)
	
	if len(vel) != max(bin_num)+1:
		print "Not enough bins provided to vel keyword"
		return

	if ax is None:
		fig, ax = plt.subplots()
	else:
		fig = plt.gcf()
	ax.set_aspect('equal')
	ax.set_zorder(10)


	if vmin is None:
		vmin = np.nanmin(vel)

	if vmax is None:
		vmax = np.nanmax(vel)

	# Flattens all arrays - I think
	xBar_pix, yBar_pix, vel = map(np.ravel, [xBar_pix, yBar_pix, vel])
	# Steps in color scale
	levels = np.linspace(vmin, vmax, ncolors)


	try:
		# VIMOS parameters
		header['CRVAL1'] = header['HIERARCH CCD1 ESO INS IFU RA']
		header['CRVAL2'] = header['HIERARCH CCD1 ESO INS IFU DEC']
		header['CD1_1'] = -abs(header['CDELT1']/(60**2))
		header['CD2_2'] = header['CDELT2']/(60**2)
	except KeyError:
		header['CD1_1'] = -abs(header['CD1_1'])
		#pass # MUSE has the correctly labelled headers

	x = (x_pix - header['CRPIX1']) * header['CD1_1'] + header['CRVAL1']
	y = (y_pix - header['CRPIX2']) * header['CD2_2'] + header['CRVAL2']

	xBar = (xBar_pix - header['CRPIX1']) * header['CD1_1'] + header['CRVAL1']
	yBar = (yBar_pix - header['CRPIX2']) * header['CD2_2'] + header['CRVAL2']

	if center is None:
		center = (max(x_pix)/2, max(y_pix)/2)

	x_label = r'$\Delta$ RA (arcsec)'
	y_label = r'$\Delta$ Dec (arcsec)'

	# Create display axis
	
	if not debug:
		ax_dis = fig.add_axes(ax.get_position(), aspect='equal', facecolor=None)
		ax_dis.set_xlim((np.array([0, header['NAXIS1']]) - (
			max(x_pix) - center[0]) + 0.5)*header['CD1_1']*60*60)
		ax_dis.set_ylim((np.array([0, header['NAXIS2']]) - center[1] - 0.5)*
			header['CD2_2']*60*60)

		ax_dis.set_xlabel(x_label)
		ax_dis.set_ylabel(y_label)

		ax_dis.minorticks_on()
		ax_dis.tick_params(length=10, which='major')
		ax_dis.tick_params(length=5, which='minor')

		tick_formatter = ticker.ScalarFormatter(useOffset=False)
		ax_dis.xaxis.set_major_formatter(tick_formatter)
		ax_dis.yaxis.set_major_formatter(tick_formatter)

		ax.ax_dis = ax_dis
		ax_dis.set_zorder(10)

		# Hide ax
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
	else:
		ax.set_xlabel('RA')
		ax.set_ylabel('Dec')


	
	bin_num = bin_num.astype(int)

	pixelSize = np.min(distance.pdist(np.column_stack([x, y])))
	xmax = (0 - header['CRPIX1']) * header['CD1_1'] + header['CRVAL1']
	xmin = (header['NAXIS1'] - header['CRPIX1']) * header['CD1_1'] + \
		header['CRVAL1']
	ymin = (0 - header['CRPIX2']) * header['CD2_2'] + header['CRVAL2']
	ymax = (header['NAXIS2'] - header['CRPIX2']) * header['CD2_2'] + \
		header['CRVAL2']

	nx = int(round((xmax - xmin)/pixelSize))
	ny = int(round((ymax - ymin)/pixelSize))
	img = np.full((nx, ny), np.nan)  # use nan for missing data
	j = np.round((x - xmin)/pixelSize).astype(int)-1
	k = np.round((y - ymin)/pixelSize).astype(int)-1
	# Clipping vel and creating 2D array of image (img)
	img[j, k] = vel[bin_num].clip(vmin, vmax)
	
	if 'cmap' not in kwargs:
		cmap=kwargs.get('cmap',sauron)
	else:
		cmap = kwg['cmap']
		if isinstance(cmap, str):
			cmap = plt.get_cmap(cmap)

	if vmin > 0: # Assume even moment
		vmin -= min((vmax-vmin)*0.05, vmin)

	# Change to RGBA 
	pic = cmap((img-vmin)/(vmax-vmin))
	if signal_noise is not None:
		if signal_noise_target is None:
			pic[j, k, 3] = (signal_noise[bin_num]/max(signal_noise))*0.5+0.5
		else:
			pic[j, k, 3] = ((signal_noise[bin_num] - signal_noise_target/2)/\
				(signal_noise_target/2)).clip(0.05,1)
	# Set bad pixels grey
	pic[np.isnan(img),:] = [0.5,0.5,0.5,1]
	ax.set_facecolor('grey')

	if alpha is not None:
		pic[j, k, 3] = alpha[bin_num]

	# RA increases right to left, thus xmax, xmin,...
	cs = ax.imshow(np.rot90(pic), interpolation='none', 
		extent=[xmax, xmin, ymin, ymax], clim=[vmin, vmax], 
		cmap=cmap) # clim and cmap supplied for colorbar
	ax.cs = cs
	
	if galaxy is not None:
		if galaxy_labelcolor is None:
			galaxy_labelcolor = 'k'
		gal_name = ax.text(0.08,0.92, galaxy, color=galaxy_labelcolor,
			verticalalignment='top',transform=ax.transAxes, bbox=dict(
			edgecolor='none', facecolor='w', zorder=99, alpha=0.6), zorder=100)
		ax.gal_name = gal_name
		if redshift is not None:
			gal_z = ax.text(0.02,0.93, "Redshift: " + str(round(redshift,3)), 
				color = galaxy_labelcolor, verticalalignment='top',
				transform=ax.transAxes)
			ax.gal_z = gal_z


	if flux is not None:
		# 1 mag contours
		ax.contour(-2.5*np.log10(flux/np.max(flux)), levels=np.arange(20), 
			colors='k', extent=[xmin, xmax, ymin, ymax])#, linewidths=1)

	if flux_unbinned is not None:
		if flux_type == 'mag':
			contours = -2.5*np.log10(np.rot90(flux_unbinned[:,::-1])/
				np.max(flux_unbinned))
			# 1 mag contours
			cont = ax.contour(contours, levels=np.arange(20), colors='k', 
				extent=[xmin, xmax, ymin, ymax])#, linewidths=1)
			cont.collections[0].set_label('Flux (mag)')


		else:
			cont = ax.contour(np.rot90(flux_unbinned[:,::-1]), colors='k', 
				extent=[xmin, xmax, ymin, ymax])#, linewidths=1)
			cont.collections[0].set_label('Flux (linear)')


	if not nodots and not show_bin_num and not show_vel:
		ax.plot(xBar, yBar, '.k',
				markersize=kwargs.get("markersize", 3))

	if show_bin_num and not show_vel:
		# for i in range(0,len(xBar),max(1,int(np.rint(np.log(len(xBar))))-3)):
		# 	ax.text(xBar[i], yBar[i], str(i), color='grey', fontsize=4)
		for i in range(len(xBar)):
			# number 100 bins only.
			if i%int((max(bin_num)+1)/100.0)==0:
				ax.text(xBar[i]+pixelSize/2, yBar[i]+pixelSize/2, str(i), 
					color='grey', fontsize=5)

	if show_vel:
		for i in range(len(xBar)):
			# number 100 bins only.
			if i%int((max(bin_num)+1)/100.0)==0:
				ax.text(xBar[i]+pixelSize/2, yBar[i]+pixelSize/2, str(vel[i]), 
					color='grey', fontsize=5)

	if redshift is not None and not debug:
		c = 299792 #km/s
		#H = 67.8 #(km/s)/Mpc # From Planck
		H = 70.0 # value used by Bolonga group.
		xlim = np.asarray(ax_dis.get_xlim())
		ylim = np.asarray(ax_dis.get_ylim())
		xlim = np.radians(xlim/(60.0*60.0)) * redshift*c/H
		ylim = np.radians(ylim/(60.0*60.0)) * redshift*c/H
		xmax = xlim[1]
		ymax = ylim[1]
		axis_label = "Distance (Mpc)"
		if max(xlim) < 1.0:
			xlim *= 1000
			ylim *= 1000
			axis_label = "Distance (kpc)"

		ax2 = fig.add_axes(ax.get_position(), aspect='equal', facecolor=None)
		ax2.minorticks_on()
		ax2.tick_params(length=10, which='major')
		ax2.tick_params(length=5, which='minor')
		ax2.tick_params(which='both',bottom=0, top=1, left=0, right=1, 
			labelbottom=0, labeltop=1, labelleft=0, labelright=1)
		
		ax2.set_ylim(ylim[0],ylim[1])
		ax2.yaxis.set_label_position("right")
		ax2.set_ylabel(axis_label, rotation=270)
		ax2.get_yaxis().get_major_formatter().set_useOffset(False)

		ax2.set_xlim(xlim[0],xlim[1])
		ax2.xaxis.set_label_position("top")
		ax2.set_xlabel(axis_label)
		ax2.get_xaxis().get_major_formatter().set_useOffset(False)

		ax.ax2 = ax2

	if title is not None:
		if hasattr(ax,'ax2'):
			ax.set_title(title, y=1.1)
		else:
			ax.set_title(title, y=1.02)

	if colorbar and not debug:
		ticks = ticker.MaxNLocator(nbins=nticks)
		ax_loc = ax.get_position()
		if hasattr(ax,'ax2'):
			cax = fig.add_axes([1.0, ax_loc.y0, 0.03, ax_loc.height])
		else:
			cax = fig.add_axes([0.93, ax_loc.y0, 0.03, ax_loc.height])
		cbar = plt.colorbar(cs, cax=cax, ticks=ticks)
	elif colorbar:
		ticks = ticker.MaxNLocator(nbins=nticks)
		if hasattr(ax,'ax2'):
			cbar = plt.colorbar(cs, ax=[ax, ax2], ticks=ticks, pad=0.1, 
				use_gridspec=True)
		else:
			cbar = plt.colorbar(cs, ax=[ax], ticks=ticks, pad=0.1, 
				use_gridspec=True)
	elif lim_labels:
		print vmin
		if lim_labels_units is None:
			ax_dis.text(1, 0, '%s/%s' % render_numbers(vmin, vmax), 
				horizontalalignment='right', rotation=270, 
				verticalalignment='bottom', bbox=dict(edgecolor='k', 
				facecolor='w', zorder=99), transform=ax.transAxes, zorder=100)
		else:
			ax_dis.text(1, 0, '%s/%s ' % render_numbers(vmin,vmax)+lim_labels_units, 
				horizontalalignment='right', rotation=270, 
				verticalalignment='bottom', bbox=dict(edgecolor='k', 
				facecolor='w', zorder=99), transform=ax.transAxes, zorder=100)


	if colorbar:	
		if label:
			t = cbar.ax.text(4.0,0.5, label, rotation=270, 
				verticalalignment='center')
			if label_size is not None:
				t.set_size(label_size*t.get_size())
		ax.cax = cbar.ax
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	ax.get_xaxis().get_major_formatter().set_useOffset(False)

	if pa is not None:
		add_orientation(pa)

	
	if debug and dms:
		plt.draw()
		from astropy.coordinates import Angle
		import astropy.units as u

		labels = ax.get_xticklabels()
		first = True
		second = True
		for i, l in enumerate(labels):
			if l.get_text() != '':
				coord = Angle(l.get_text().strip('$'), unit=u.deg)
				if first:
					labels[i] = r"%.0f$\degree$ %.0f' %.2f" % (coord.dms[0], 
						abs(coord.dms[1]), abs(coord.dms[2])) + '"'
					first = False
				elif second:
					labels[i] = ''
					second = False
				else:
					labels[i] = '%.2f"' % (coord.dms[2])
		ax.set_xticklabels(labels)

		labels = ax.get_yticklabels()
		first = True
		for i, l in enumerate(labels):
			if l.get_text() != '':
				coord = Angle(l.get_text().strip('$'), unit=u.deg)
				if i == len(labels)-2:
					labels[i] = r"%.0f$\degree$ %.0f' %.2f" %(coord.dms[0], 
						abs(coord.dms[1]), abs(coord.dms[2])) + '"'
					first = False
				else:
					labels[i] = '%.2f"' % (coord.dms[2])
		ax.set_yticklabels(labels)
	
	if save is not None:
		ax.saveTo = save
		if not os.path.exists(os.path.dirname(save)):
			os.makedirs(os.path.dirname(save))  

		plt.savefig(save, bbox_inches="tight")

	if close:
		plt.close(fig)

	return ax