## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to read the cube format and print it in
## table form into a text file.
## testing a change at NAM2015
## warrenj 20150727 Changing to a python script

import cap_plot_velfield# as plot_velfield
import numpy as np # for reading files
#import numpy
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
#-----------------------------------------------------------------------------

galaxy = "ngc3557"

tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output.txt"
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
number_of_bins = max(bin_num)
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

#n_spaxels = len(galaxy_data) * len(galaxy_data[0])




flux = np.zeros((number_of_bins+1))
# ----------========== Spatially Binning =============---------
for spaxel in range(n_spaxels):
    flux[bin_num[spaxel]] += np.sum(galaxy_data[:,x[spaxel],y[spaxel]])

flux = flux/np.median(flux)
#print flux



        

#v_binned = v_binned - MEDIAN(v_binned)



b = []
[b.append(i) for i in range(len(bin_num)) if bin_num[i] not in bin_num[b]]
# b contains the indices of bin_num which correspond the first time that a bin number occurs. 

xNode = []
yNode = []

for i in b:
    xNode.append(xBin[i])
    yNode.append(yBin[i])



#v = np.zeros((n_spaxels))
v = []
flux_unbinned = []
for spaxel in range(n_spaxels):
#    v[spaxel] = v_binned[bin_num[spaxel]]
     v.append(v_binned[bin_num[spaxel]])
     flux_unbinned.append(np.sum(galaxy_data[:,x[spaxel],y[spaxel]]))

flux_unbinned = flux_unbinned/np.median(flux_unbinned)


# automatically uses sauron colormap
plt.clf()
plt.title('Velocity')
#cap_plot_velfield.plot_velfield(xNode, yNode, v_binned, nodots=True, flux=flux)
cap_plot_velfield.plot_velfield(x, y, v, nodots=True, flux=flux_unbinned)
plt.show()




##example use of contour
# Create a simple dataset:
#data = RANDOMU(seed, 9, 9)

## Plot the unsmoothed data:
#unsmooth = CONTOUR(data, TITLE='Unsmoothed', $
#   LAYOUT=[2,1,1], RGB_TABLE=13, /FILL, N_LEVELS=10)
## Draw the outline of the 10 levels
#outline1 = CONTOUR(data, N_LEVELS=10, /OVERPLOT)
 
## Plot the smoothed data:
#smooth = CONTOUR(MIN_CURVE_SURF(data), TITLE='Smoothed', $
#   /CURRENT, LAYOUT=[2,1,2], RGB_TABLE=13, $
#   /FILL, N_LEVELS=10)
# Draw the outline of the 10 levels
#outline2 = CONTOUR(MIN_CURVE_SURF(data), $
#   N_LEVELS=10, /OVERPLOT)

