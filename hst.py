## ==================================================================
## Unsharp mask HST images
## ==================================================================
## warrenj 20151020 Routine to apply an unsharp mask on various 
## scales to HST images.



import matplotlib.pyplot as plt # used for plotting
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import numpy as np # for array handling
from sauron_colormap import sauron
from matplotlib.ticker import MaxNLocator
from scipy import ndimage





#-----------------------------------------------------------------------------

def unsharpmask(galaxy, level, xshift, yshift, size, **kwargs):

    plot = True


    nticks=7

    img_file = glob.glob("/Data/vimosindi/hst/%s/*.fits" %(galaxy))

    img, hst_header = pyfits.getdata(img_file[0], 0, header=True)

    pa = hst_header['OORIENTA']
    img = ndimage.rotate(img, pa, reshape=False)

    s = np.shape(img)


    img_unsharp = np.zeros((size, size))
#    s=np.shape(img_unsharp)

    for x in range(np.shape(img_unsharp)[0]):
        for y in range(np.shape(img_unsharp)[1]):
            img_unsharp[x,y] = img[s[0]/2+xshift-size/2+x,
                                   s[1]/2+yshift-size/2+y] - \
                np.median(img[s[0]/2+xshift-size/2+x-level:
                              s[0]/2+xshift-size/2+x+level, 
                              s[1]/2+yshift-size/2+y-level:
                              s[1]/2+yshift-size/2+y+level])
    img_unsharp -= np.min(img_unsharp)
    img_unsharp /= np.max(img_unsharp)

    fig, ax = plt.subplots(nrows=1,ncols=1)
    imgplot=plt.imshow(np.rot90(np.log(img_unsharp)))#, interpolation='none', cmap=kwargs.get('cmap',sauron))

    if 'WFPC2,1,' in hst_header['PHOTMODE']:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
#        ax.set_xaxis(x*0.046) # res = 0.046"per px
        ax.set_yticks(np.arange(ymin*0.046, ymax*0.046,xmax*0.046/2))
        ax.set_xlabel("Angular Size (arcsec)")
        ax.set_ylabel("Angular Size (arcsec)")
    elif 'WFPC2,' in hst_header['PHOTMODE']:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim(xlim*0.1) # res = 0.1"per px
        ax.set_ylim(ylim*0.1)
        ax.set_xlabel("Angular Size (arcsec)")
        ax.set_ylabel("Angular Size (arcsec)")





    colorbar=False

    if colorbar:
#        divider = make_axes_locatable(ax)
#        cax = divider.append_axes("right", size="5%", pad=0.7)
        fig.subplots_adjust(right=0.66, top=0.84)
        cax = fig.add_axes([0.75, 0.1, 0.02, 0.74])
## symmetric should make VD plots odd... ******************************
        ticks = MaxNLocator(nbins=nticks)#, symmetric=True)
        cbar = fig.colorbar(imgplot, cax=cax, ticks=ticks)
#        plt.clim(vmin,vmax)  # make color axis symmetrical


#    imgplot.set_cmap('hot')
    plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), color='white')
#            verticalalignment='top',transform=ax.transAxes)

    if plot:
        plt.show()
    plt.savefig("/Data/vimosindi/hst/%s/%s_unsharp.png" %(galaxy,galaxy))

    plt.close('all')



    imgplot=plt.imshow(np.rot90(np.log(img)))
    plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), color='white')

    if plot:
        plt.show()
    plt.savefig("/Data/vimosindi/hst/%s/%s.png" %(galaxy,galaxy))    

    plt.close('all')






























##############################################################################

# Use of hst.py

if __name__ == '__main__':
    galaxies = ["ic1459", "ngc1316", "ic4296", "ngc3557", "ngc1399", "ngc0612"]
    levels = [6,6,6,6,6]
    xshifts = [10,-20,-10,20,-63]
    yshifts = [34,0,30,0,3]
    sizes = [50,50,50,60,50]
    
    for gal in range(1):


        unsharpmask(galaxies[gal], levels[gal], xshifts[gal], yshifts[gal], sizes[gal])
