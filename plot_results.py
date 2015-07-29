## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to read the cube format and print it in
## table form into a text file.
## testing a change at NAM2015
## warrenj 20150727 Changing to a python script

from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for reading files
#import numpy
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
#-----------------------------------------------------------------------------

galaxy = "ngc3557"
discard = 2

tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output.txt"
tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output2.txt"
output_v = "/Data/vimosindi/analysis/%s/results/gal_vel.dat" % (galaxy)
output_temp_weighting = "/Data/vimosindi/analysis/%s/" % (galaxy) +\
"results/template_weighting.dat" 
output_sigma = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"gal_sigma.dat"
output_h3 = "/Data/vimosindi/analysis/%s/results/gal_h3.dat" % (galaxy)
output_h4 = "/Data/vimosindi/analysis/%s/results/gal_h4.dat" % (galaxy)
output_h5 = "/Data/vimosindi/analysis/%s/results/gal_h5.dat" % (galaxy)
output_h6 = "/Data/vimosindi/analysis/%s/results/gal_h6.dat" % (galaxy)
output_Chi = "/Data/vimosindi/analysis/%s/results/gal_Chi.dat" % (galaxy)

# Read tessellation file
x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
    skiprows = 1) 
n_spaxels = len(bin_num)
number_of_bins = int(max(bin_num)+1)
order = bin_num.argsort()

# Read results files - each entry in array corresponds to a bin (not
# a spaxel)
v_binned = np.loadtxt(output_v)



# ------------========== Total flux per bin ===========----------
# ----------========= Reading the spectrum  =============---------

# FILE_SEARCH returns an array even in cases where it only returns
# one result. This is NOT equivalent to a scalar. 
dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
    "*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits" % (galaxy)) 

galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)

s = galaxy_data.shape
rows_to_remove = range(discard)
rows_to_remove.extend([s[1]-1-i for i in range(discard)])
cols_to_remove = range(discard)
cols_to_remove.extend([s[2]-1-i for i in range(discard)])

galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)



## +
#flux = np.zeros((number_of_bins))
## ----------========== Spatially Binning =============---------
#for spaxel in range(n_spaxels):
#    flux[bin_num[spaxel]] += np.sum(galaxy_data[:,y[spaxel],x[spaxel]])
#
#flux = flux/np.median(flux)
#
#
#
#b = []
## b contains the indices of bin_num which correspond the first time that a bin number occurs. 
#[b.append(i) for i in range(len(bin_num)) if bin_num[i] not in bin_num[b]]
#
#xNode = []
#yNode = []
#
#for i in b:
#    xNode.append(xBin[i])
#    yNode.append(yBin[i])
## -


#v_binned = v_binned - np.median(v_binned)




# ------------=========== unbinned version ============----------
v = []
flux_unbinned = []
for spaxel in range(n_spaxels):
     v.append(v_binned[bin_num[spaxel]])
     flux_unbinned.append(np.sum(galaxy_data[:,y[spaxel],x[spaxel]]))

flux_unbinned = flux_unbinned/np.median(flux_unbinned)




# ------------========== Different binning ===========----------
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





# ------------============= Plot velfield ==============----------
# automatically uses sauron colormap
plt.clf()
plt.title('Velocity')
#plot_velfield(xNode, yNode, v_binned, nodots=True, flux=flux)
plot_velfield(xBar, yBar, v_binned, nodots=False, flux=flux_bar_binned)
#plot_velfield(x, y, v, nodots=False, flux=flux_unbinned)
plt.show()


