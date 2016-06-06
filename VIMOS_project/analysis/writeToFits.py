import numpy as np
import matplotlib.pyplot as plt # used for plotting
from sauron_colormap import sauron
import glob # for searching for files
import pyfits # reads fits files (is from astropy)




def writeToFits(**kwargs):
    vLimit=2
    discard =2
    plot = "stellar_sigma"
    file_dir = "/Data/vimosindi/analysis_sav_2016-02-23/ic1459/results/4200-/gal_%s.dat" %(plot)


    v_binned, v_uncert_binned = np.loadtxt(file_dir, unpack=True)

# Read tessellation file
    tessellation_File = "/Data/vimosindi/analysis_sav_2016-02-23/ic1459/voronoi_2d_binning_output.txt"
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1) 
    n_spaxels = len(bin_num)
    number_of_bins = int(max(bin_num)+1)
    order = bin_num.argsort()



    
    img = np.zeros([36,36])
    dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
        "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % ("ic1459")) 


    galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)

    s = galaxy_data.shape
    rows_to_remove = range(discard)
    rows_to_remove.extend([s[1]-1-i for i in range(discard)])
    cols_to_remove = range(discard)
    cols_to_remove.extend([s[2]-1-i for i in range(discard)])

    galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
    galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)



    galaxy_data_unbinned = np.sum(galaxy_data, axis=0)
    

    # Asign v to every spaxel in bin
    v_unbinned = np.zeros(galaxy_data_unbinned.shape)
    v_uncert_unbinned = np.zeros(galaxy_data_unbinned.shape)
    for spaxel in range(n_spaxels):
        v_unbinned[x[spaxel],y[spaxel]] = v_binned[bin_num[spaxel]]
        v_uncert_unbinned[x[spaxel],y[spaxel]] = \
            v_uncert_binned[bin_num[spaxel]]

    if "vel" in plot:
        lwv = v_unbinned*galaxy_data_unbinned
        v_binned -= np.mean(lwv)*n_spaxels/np.sum(galaxy_data_unbinned)
                
    for i in range(len(x)):
        img[x[i],y[i]]=v_binned[bin_num[i]]

    vmax = max(v_binned)
    vmin = min(v_binned)
    v_sorted = sorted(np.unique(v_binned))
    if len(v_sorted) < 2*vLimit:
        v_sorted = sorted(v_binned)
    vmin = v_sorted[vLimit]
    vmax = v_sorted[-vLimit-1]
# Make velocity fields symmetric
    if "vel" in plot:
        if abs(vmin)>vmax:
            if vmin*vmax>0:
                vmin, vmax = vmax, vmin
            else:
                vmin=-abs(vmax)
        else:
            vmax=abs(vmin)

                
#    plt.imshow(np.rot90(img), interpolation='none', cmap=kwargs.get('cmap',sauron))
#    plt.clim(vmin,vmax)

#    plt.show()



    hdu = pyfits.PrimaryHDU(np.rot90(img))
    hdu.writeto("/Data/vimosindi/analysis_sav_2016-02-23/ic1459/results/4200-/gal_%s.fits" % (plot))

























##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

    writeToFits()




