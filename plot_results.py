## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to read the cube format and print it in
## table form into a text file.
## testing a change at NAM2015
## warrenj 20150727 Changing to a python script

from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
#-----------------------------------------------------------------------------

wav_range=None

plot = "v"
#plot ="sigma"
#plot ="h3"
#plot="h4"
wav_range="4200-"


galaxy = "ngc3557"
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

# Read tessellation file
x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
    skiprows = 1) 
n_spaxels = len(bin_num)
number_of_bins = int(max(bin_num)+1)
order = bin_num.argsort()

# Read results files - each entry in array corresponds to a bin (not
# a spaxel)
if plot=="v":
    v_binned = np.loadtxt(output_v)
    print "Velocity"
    v_binned += -1.1*np.mean(v_binned)
if plot=="sigma":
    v_binned = np.loadtxt(output_sigma)
    print "Velcity dispersion"
if plot=="h3":
    v_binned = np.loadtxt(output_h3)
    print "h3"
if plot=="h4":
    v_binned = np.loadtxt(output_h4)
    print "h4"

#v_binned += -np.median(v_binned)


# ------------========== Total flux per bin ===========----------
# ----------========= Reading the spectrum  =============---------

# FILE_SEARCH returns an array even in cases where it only returns
# one result. This is NOT equivalent to a scalar. 
dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
    "*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits" % (galaxy)) 

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






# ------------=========== unbinned version ============----------
v = []
flux_unbinned = []
for spaxel in range(n_spaxels):
     v.append(v_binned[bin_num[spaxel]])
     flux_unbinned.append(np.sum(galaxy_data[:,y[spaxel],x[spaxel]]))

flux_unbinned = flux_unbinned/np.median(flux_unbinned)




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


# ------------============ Setting v range =============----------
vmax = max(v_binned)
vmin = min(v_binned)
v_sorted = sorted(np.unique(v_binned))
#v_sorted = sorted(v_binned)
vmin = v_sorted[vLimit]
vmax = v_sorted[-vLimit-1]


# ------------============= Plot velfield ==============----------
# automatically uses sauron colormap
plt.clf()
if plot=="v":
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
#plot_velfield(x, y, v, vmin=vmin, vmax=vmax, 
#    nodots=False, flux=flux_unbinned)


plt.savefig("/home/warrenj/Desktop/" + plot + "_field_" + wav_range + ".png", bbox_inches="tight")
plt.show()


