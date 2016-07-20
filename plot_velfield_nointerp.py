## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20151015 Routine to plot velocity (or similar) fields using
## methods from Michele Cappellari's plot_velfield routine and his 
## voronoi binning routine (both contained in the voronoi binning package),
## with some of my own inventions.
# warrenj 20160209 Added show_bin_num keyword.
# warrenj 20160413 Added flux_type keyword



import numpy as np # for reading files
from scipy.spatial import distance
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from sauron_colormap import sauron
import math
import matplotlib.pyplot as plt # used for plotting
import os



def plot_velfield_nointerp(x_pix, y_pix, bin_num, xBar_pix, yBar_pix, vel, 
    vmin=None, vmax=None, nodots=False, colorbar=False, label=None, flux=None, 
    flux_unbinned=None, galaxy = None, redshift = None, nticks=7, 
    ncolors=64, title=None, save=None, show_bin_num=False, flux_type='mag',
    ax = None, **kwargs):

    kwg = {}
    kwg.update(kwargs)

    
    if len(vel) != max(bin_num)+1:
        print "Not enough bins provided to vel keyword"
        return

    if ax is None:
        fig, ax = plt.subplots(nrows=1,ncols=1)
    
    if title is not None:
        plt.title(title, y=1.1)

    if vmin is None:
        vmin = np.min(vel)

    if vmax is None:
        vmax = np.max(vel)

    xBar_pix, yBar_pix, vel = map(np.ravel, [xBar_pix, yBar_pix, vel])
    levels = np.linspace(vmin, vmax, ncolors)

    res = 0.67 #arcsec per pixel
    xBar = xBar_pix*res
    yBar = yBar_pix*res
    x = x_pix*res
    y = y_pix*res

    x -= max(x)/2
    y -= max(y)/2
    axis_label = "Angular Size (arcsec)"


    im_xBar = np.copy(xBar)
    im_yBar = np.copy(yBar)
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

    cs = ax.imshow(np.rot90(img[:,:]), interpolation='none', 
        cmap=cmap,extent=[xmin - pixelSize/2, 
        xmax + pixelSize/2, ymin - pixelSize/2, ymax + pixelSize/2],
        clim = (vmin,vmax))
#    ax.clim(vmin,vmax)
    ax.invert_xaxis()

#    plt.gca().invert_yaxis()

    ax.set_ylabel(axis_label)
    ax.set_xlabel(axis_label)
    ax.axis('tight')  # Equal axes and no rescaling
    ax.minorticks_on()
    ax.tick_params(length=10, which='major')
    ax.tick_params(length=5, which='minor')

    ax2 = ax.twinx()
    ax3 = ax.twiny()

    xmin_sav, xmax_sav = ax.get_xlim()
    ymin_sav, ymax_sav = ax.get_ylim()
    xlim=np.array([xmin_sav, xmax_sav])
    ylim=np.array([ymin_sav, ymax_sav])



    if redshift is not None:
        c = 299792 #km/s
#        H = 67.8 #(km/s)/Mpc
#        H = 70.0 # value used by Bolonga group.
        xlim = np.radians(xlim/(60.0*60.0)) * redshift*c/100
        ylim = np.radians(ylim/(60.0*60.0)) * redshift*c/100
        xmax = xlim[1]
        ymax = ylim[1]
#        xlim -= xmax/2
#        ylim -= ymax/2
        axis_label = "Distance (Mpc/h)"
        if max(xlim) < 1.0:
            xlim *= 1000
            ylim *= 1000
            axis_label = "Distance (kpc/h)"
 

    ax2.minorticks_on()
    ax2.tick_params(length=10, which='major')
    ax2.tick_params(length=5, which='minor')
    ax2.set_ylim(ylim[0],ylim[1])
    ax2.set_ylabel(axis_label, rotation=270)

    ax3.minorticks_on()
    ax3.tick_params(length=10, which='major')
    ax3.tick_params(length=5, which='minor')
    ax3.set_xlim(xlim[0],xlim[1])
    ax3.set_xlabel(axis_label)





#    y1, y2 = ax.get_ylim()

        

    if galaxy is not None:
        plt.text(0.02,0.98, "Galaxy: " + galaxy, color='black',
            verticalalignment='top',transform=ax.transAxes)
        if redshift is not None:
            plt.text(0.02,0.93, "Redshift: " + str(round(redshift,3)), 
                color = 'black',verticalalignment='top',
                transform=ax.transAxes)



    if flux is not None:
        ax.tricontour(x[::-1], y[::-1],
                      -2.5*np.log10(flux.ravel()/np.max(flux)),
                      levels=np.arange(20), colors='k') # 1 mag contours

# NB: have assumed a square image!!!!
    if flux_unbinned is not None:
#        flux_unbinned[477]=0.001
        if flux_type == 'mag':
            contours = -2.5*np.log10(flux_unbinned.ravel()/
                                     np.max(flux_unbinned))
# 1 mag contours
            ax.tricontour(x[::-1], y[::-1], contours, levels=np.arange(20),
                      colors='k')
        else:
            ax.contour(np.reshape(y,np.shape(flux_unbinned)),
                       np.reshape(x,np.shape(flux_unbinned)), 
                       flux_unbinned, colors='k')

    if not nodots and not show_bin_num:
#**********************************################
#        xBar /= res # no idea why this needs removing...
#        yBar /= res
#**********************************################
        ax.plot(ax.get_xlim()[1]-xBar, ax.get_ylim()[1]-yBar, '.k',
                markersize=kwargs.get("markersize", 3))

    if show_bin_num:
        for i in range(len(xBar)):
            ax.text(ax.get_ylim()[0]+yBar[i],
                    ax.get_xlim()[0]+xBar[i], str(i), color='grey',
#                    fontsize='xx-small')
                    fontsize=5)

        
    if colorbar:
#        divider = make_axes_locatable(ax)
#        cax = divider.append_axes("right", size="5%", pad=0.7)
        fig.subplots_adjust(right=0.66, top=0.84)
        cax = fig.add_axes([0.75, 0.1, 0.02, 0.74])
## symmetric should make VD plots odd... ******************************
        ticks = MaxNLocator(nbins=nticks)#, symmetric=True)
        cbar = fig.colorbar(cs, cax=cax, ticks=ticks)
#        plt.clim(vmin,vmax)  # make color axis symmetrical
        if label:
            cbar.set_label(label, rotation=270)




    if save is not None:
        if not os.path.exists(os.path.dirname(save)):
            os.makedirs(os.path.dirname(save))
        plt.savefig(save, bbox_inches="tight")

    return ax



