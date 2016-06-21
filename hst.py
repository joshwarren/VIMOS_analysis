## ==================================================================
## Unsharp mask HST images
## ==================================================================
## warrenj 20151020 Routine to apply an unsharp mask on various 
## scales to HST images.



import matplotlib.pyplot as plt # used for plotting
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
import numpy as np # for array handling
from sauron_colormap import sauron
from matplotlib.ticker import MaxNLocator
from scipy import ndimage





#-----------------------------------------------------------------------------

def unsharpmask(galaxy, level, **kwargs):

    plot = False


    nticks=7
    snticks=5
    fov=20 #arcsec

    fig, ax = plt.subplots(nrows=1,ncols=1)
   
    img_file = glob.glob("/Data/hst/%s/*.fits" %(galaxy))

    img, hst_header = pyfits.getdata(img_file[0], 0, header=True)
  
    pa = hst_header['OORIENTA']
    img = ndimage.rotate(img, pa, reshape=False)
    img_sav = np.array(img)

    s = np.shape(img)
    
    dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
        "*_cube.fits" % (galaxy))
    galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)

    hst_x = np.arange(hst_header['NAXIS1'])[::-1]*hst_header['CD1_1']
    hst_y = np.arange(hst_header['NAXIS2'])*hst_header['CD2_2']

    hst_x += -hst_header['CD1_1']*hst_header['CRPIX1'] + header['CRVAL1'] - \
            hst_header['CRVAL1']
    hst_y -= hst_header['CD2_2']*hst_header['CRPIX2']+ header['CRVAL2'] - \
            hst_header['CRVAL2']


    hst_x *=60*60
    hst_y *=60*60

    f = np.where(np.abs(hst_x) < fov/2)[0]
    hst_x = hst_x[f]
    img = img[f]
    f2 = np.where(np.abs(hst_y) < fov/2)[0]
    hst_y = hst_y[f2]
    img = img[:,f2]
    
    imgplot=plt.imshow(np.rot90(np.log(img)))#, interpolation='none', cmap=kwargs.get('cmap',sauron))
    x = range(0, len(hst_x), len(hst_x)/snticks)
    plt.xticks(x, np.round(hst_x[x],1))
    y = range(0, len(hst_y), len(hst_y)/snticks)
    plt.yticks(y,np.round(-hst_y[y],1))
    ax.set_xlabel("Angular Size (arcsec)")
    ax.set_ylabel("Angular Size (arcsec)")
    
    plt.text(0.02,0.9, "Galaxy: " + galaxy.upper(), color='white')

    if plot:
        plt.show()
    plt.savefig("/Data/hst/%s/%s.png" %(galaxy,galaxy))    
    plt.close('all')




    img_unsharp = img_sav
    for x in range(level, s[0]-level):
        for y in range(level, s[1]-level):
            img_unsharp[x,y] -= \
                np.median(img_unsharp[x-level:x+level, y-level:y+level])
#    img_unsharp -= np.min(img_unsharp)
#    img_unsharp /= np.max(img_unsharp)

    img_unsharp = img_unsharp[f]
    img_unsharp = img_unsharp[:,f2]
    
    imgplot=plt.imshow(np.rot90(np.log(img_unsharp)))#, interpolation='none', cmap=kwargs.get('cmap',sauron))
    x = range(0, len(hst_x), len(hst_x)/snticks)
    plt.xticks(x, np.round(hst_x[x],1))
    y = range(0, len(hst_y), len(hst_y)/snticks)
    plt.yticks(y,np.round(-hst_y[y],1))

    ax.set_xlabel("Angular Size (arcsec)")
    ax.set_ylabel("Angular Size (arcsec)")


#    imgplot.set_cmap('hot')
    plt.text(0.02,0.9, "Galaxy: " + galaxy.upper(), color='white')
#            verticalalignment='top',transform=ax.transAxes)

    if plot:
        plt.show()
    plt.savefig("/Data/hst/%s/%s_unsharp.png" %(galaxy,galaxy))

    plt.close('all')




























##############################################################################

# Use of hst.py

if __name__ == '__main__':
# Removed ngc1316 (Fornax A) from the list and ngc0612
    galaxies = ["ic1459", "ic4296", "ngc3557", "ngc1399"] 
    levels = [6,6,6,6,6]
    xshifts = [10,-20,-10,20,-63]
    yshifts = [34,0,30,0,3]
    sizes = [50,50,50,60,50]
    
    for gal in range(6):

        print galaxies[gal]

        unsharpmask(galaxies[gal], levels[gal])
