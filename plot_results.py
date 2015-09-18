## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.

from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
#-----------------------------------------------------------------------------

wav_range=None
wav_range="4200-"
show_plot = True

#galaxy = "ngc3557"
galaxy="ic1459"
discard = 2 # rows of pixels to discard- must have been the same 
            #    for all routines 
vLimit = 2 # limit the velocity plotted: this is the number of 
           #     unique velocity values that will be ignored  

if wav_range:
    wav_range_dir = wav_range + "/"
else:
    wav_range_dir = ""

tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output.txt"
tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output2.txt"
output_v = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_vel.dat" % (wav_range_dir)
output_temp_weighting = "/Data/vimosindi/analysis/%s/" % (galaxy) +\
"results/%stemplate_weighting.dat" % (wav_range_dir)
output_sigma = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_sigma.dat" % (wav_range_dir)
output_h3 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_h3.dat" % (wav_range_dir)
output_h4 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_h4.dat" % (wav_range_dir)
output_h5 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_h5.dat" % (wav_range_dir)
output_h6 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_h6.dat" % (wav_range_dir)
output_Chi = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_Chi.dat" % (wav_range_dir)
output_OIII = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_OIII.dat" % (wav_range_dir)
output_NI = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_NI.dat" % (wav_range_dir)
output_Hb = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_Hb.dat" % (wav_range_dir)
output_Hd = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_Hd.dat" % (wav_range_dir)

plots = {"v" : output_v, "sigma" : output_sigma, "h3" : output_h3, "h4" : output_h4, "OIII" : output_OIII, "NI" : output_NI, "Hb" : output_Hb, "Hd" : output_Hd}


# Read tessellation file
x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
    skiprows = 1) 
n_spaxels = len(bin_num)
number_of_bins = int(max(bin_num)+1)
order = bin_num.argsort()




# ------------========== Total flux per bin ===========----------
# ----------========= Reading the spectrum  =============---------

# FILE_SEARCH returns an array even in cases where it only returns
# one result. This is NOT equivalent to a scalar. 
dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
    "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % (galaxy)) 

## Directory for plotting single quadrant
#dataCubeDirectory = glob.glob("/Data/vimosindi/%s-3/Q2/calibrated/cube/" \
#    "*_fluxcal_cube.fits" % (galaxy)) 

galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)

s = galaxy_data.shape
rows_to_remove = range(discard)
rows_to_remove.extend([s[1]-1-i for i in range(discard)])
cols_to_remove = range(discard)
cols_to_remove.extend([s[2]-1-i for i in range(discard)])

galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)







# ------------========== Spatially binning ===========----------
xBar, yBar = np.loadtxt(tessellation_File2, unpack=True, 
    skiprows = 1) 
flux_bar_binned = np.zeros((number_of_bins))
n_spaxels_in_bin = np.zeros((number_of_bins))

for spaxel in range(n_spaxels):
    flux_bar_binned[int(bin_num[spaxel])] += np.sum(
        galaxy_data[:,y[spaxel],x[spaxel]])
    n_spaxels_in_bin[int(bin_num[spaxel])] += 1

for bin in range(number_of_bins):
    flux_bar_binned[bin] = flux_bar_binned[bin]/n_spaxels_in_bin[bin]

flux_bar_binned = flux_bar_binned/np.median(flux_bar_binned)


# ------------============= Read fields =============----------
# Read results files - each entry in array corresponds to a bin (not
# a spaxel)
for plot in plots:
    print plot
    v_binned = np.loadtxt(plots[plot])

    if plot=="v":
        v_binned += -1.1*np.median(v_binned)
    if plot=="OIII" or plot=="NI" or plot=="Hb" or plot=="Hd":
        v_binned += -np.median(v_binned)

# ------------============ Setting v range =============----------
    vmax = max(v_binned)
    vmin = min(v_binned)
    v_sorted = sorted(np.unique(v_binned))
    if len(v_sorted) < 2*vLimit:
        v_sorted = sorted(v_binned)
    vmin = v_sorted[vLimit]
    vmax = v_sorted[-vLimit-1]


# ------------============= Plot velfield ==============----------
# automatically uses sauron colormap
    plt.clf()
    if plot=="v" or plot=="OIII" or plot=="NI" or plot=="Hb" or plot=="Hd":
        plt.title('Velocity Map')
        CBLabel = "LOSV (km s$^{-1}$)"
    elif plot=="sigma":
        plt.title('Velocity Dispersion Map')
        CBLabel = "LOSVD (km s$^{-1}$)"
    else:
        plt.title(plot + ' Map')
        CBLabel = ""

    plot_velfield(xBar, yBar, v_binned, vmin=vmin, vmax=vmax, 
        nodots=False, colorbar=True, label=CBLabel, flux=flux_bar_binned)

    plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%splots/%s_field_%s.png" % (wav_range_dir, plot, wav_range), \
        bbox_inches="tight")
    if show_plot:
        plt.show()


