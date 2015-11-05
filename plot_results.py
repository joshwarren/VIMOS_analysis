## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.

from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for array handling
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 



#-----------------------------------------------------------------------------

def plot_results(galaxy, discard=0, wav_range="", vLimit=2, plots=False, 
    nointerp=False, **kwargs):

    data_file =  "/Data/vimosindi/analysis/galaxies.txt"
    # different data types need to be read separetly
    z_gals = np.loadtxt(data_file, skiprows=1, usecols=(1,))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    z = z_gals[i_gal]

# Normalisation methods:
#  lum: velocity of brightest pixel is set to 0km/s
#  lwv: luminosity weighted velocity mean is zero.
    norm = "lwv" 
#    norm = "lum"



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
    output_v_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_vel_uncert.dat" % (wav_range_dir)
    output_temp_weighting = "/Data/vimosindi/analysis/%s/" % (galaxy) +\
        "results/%stemplate_weighting.dat" % (wav_range_dir)
    output_sigma = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_sigma.dat" % (wav_range_dir)
    output_sigma_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_sigma_uncert.dat" % (wav_range_dir)
    output_h3 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h3.dat" % (wav_range_dir)
    output_h3_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h3_uncert.dat" % (wav_range_dir)
    output_h4 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h3.dat" % (wav_range_dir)
    output_h4_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h4_uncert.dat" % (wav_range_dir)
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
    

    outputs = {"v" : output_v, "v_uncert" : output_v_uncert, 
        "sigma" : output_sigma, "sigma_uncert" : output_sigma_uncert, 
        "h3" : output_h3, "h3_uncert" : output_h3_uncert, 
        "h4" : output_h4, "h4_uncert" : output_h4_uncert, 
        "OIII" : output_OIII, "NI" : output_NI, 
        "Hb" : output_Hb, "Hd" : output_Hd}
#    outputs = {"v":output_v}
    outputs = {"h3_uncert":output_h3_uncert}

# Read tessellation file
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1) 
    n_spaxels = len(bin_num)
    number_of_bins = int(max(bin_num)+1)
    order = bin_num.argsort()

# Read galaxies.txt file
    data_file =  "/Data/vimosindi/analysis/galaxies.txt"
    # different data types need to be read separetly
    x_gals, y_gals = np.loadtxt(data_file, 
        unpack=True, skiprows=1, usecols=(4,5))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]

    center_bin = bin_num[x_gals[i_gal]*(max(y)+1) + y_gals[i_gal]]

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



    galaxy_data_unbinned = np.sum(galaxy_data, axis=0)
#    galaxy_data_unbinned = galaxy_data_unbinned.flatten()

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
    for plot in outputs:
        print plot
        if plot=="v" or plot=="sigma" or plot=="h3" or plot=="h4" or \
            plot=="OIII" or plot=="NI" or plot=="Hb" or plot=="Hd":
            v_binned = np.loadtxt(outputs[plot])#, skiprows=1)
        else:
            v_binned = np.loadtxt(outputs[plot], usecols=(1,), unpack=True)#, skiprows=1)

#        if plot=="v":
#            v_binned += -1.1*np.median(v_binned)
#        if plot=="OIII" or plot=="NI" or plot=="Hb" or plot=="Hd":
#            v_binned += -np.median(v_binned)
# Asign v to every spaxel in bin
        v_unbinned = np.zeros(galaxy_data_unbinned.shape)
        for spaxel in range(n_spaxels):
            v_unbinned[x[spaxel],y[spaxel]] = v_binned[bin_num[spaxel]]
# ------------============ Setting v range =============----------
        if plot=="v" or plot=="v_uncert" or plot=="OIII" or plot=="NI" \
            or plot=="Hb" or plot=="Hd":
            if norm == "lum":
                v_binned -= v_binned[center_bin]
            if norm == "lwv":
#                galaxy_data_unbinned1=galaxy_data_unbinned/np.median(
#                    galaxy_data_unbinned)
                lwv = v_unbinned*galaxy_data_unbinned

                v_binned -= np.mean(lwv)*n_spaxels/np.sum(galaxy_data_unbinned)
#                v_binned -= np.mean(v_binned)


        vmax = max(v_binned)
        vmin = min(v_binned)
        v_sorted = sorted(np.unique(v_binned))
        if len(v_sorted) < 2*vLimit:
            v_sorted = sorted(v_binned)
        vmin = v_sorted[vLimit]
        vmax = v_sorted[-vLimit-1]
        if plot=="v" or plot=="v_uncert" or plot=="OIII" or plot=="NI" \
            or plot=="Hb" or plot=="Hd":
            if abs(vmin)<vmax:
                vmin=-vmax
            else:
                vmax=-vmin

# ------------============= Plot velfield ==============----------
# automatically uses sauron colormap
        if plot=="v":
            title = 'Stellar Velocity Map'
            CBLabel = "LOSV (km s$^{-1}$)"
        elif plot=="v_uncert":
            title = "Stellar Velocity Uncertainty Map"
            CBLabel = "LOSV (km s$^{-1}$)"
        elif plot=="OIII" or plot=="NI" or plot=="Hb" or plot=="Hd":
            title = plot + ' Velocity Map'
            CBLabel = "LOSV (km s$^{-1}$)"
        elif plot=="sigma":
            title = 'Velocity Dispersion Map'
            CBLabel = "LOSVD (km s$^{-1}$)"
        elif plot=="v_uncert":
            title = "Stellar Velocity Dispersion Uncertainty Map"
            CBLabel = "LOSV (km s$^{-1}$)"
        elif plot=="h3_uncert":
            title = "h3 Uncertainty Map"
        elif plot=="h4_uncert":
            title = "h4 Uncertainty Map"
        else:
            title = plot + ' Map'
            CBLabel = ""
  

# ------------===== Plot velfield - no interperlation ======----------
        if nointerp:
            saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
                plot, wav_range)
            plot_velfield_nointerp(x, y, bin_num, xBar, yBar, v_binned, 
                vmin=vmin, vmax=vmax, nodots=False, colorbar=True, 
                label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
                galaxy = galaxy.upper(), redshift = z, title=title, 
                save=saveTo)
#            plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
#                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
#                plot, wav_range), bbox_inches="tight")

# ------------===== Plot velfield - with interperlation ====----------
        else:
           saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/%s_field_%s.png" % (wav_range_dir, plot, wav_range)
           plot_velfield(xBar, yBar, v_binned, vmin=vmin, vmax=vmax, 
                nodots=False, colorbar=True, label=CBLabel, 
                flux_unbinned=galaxy_data_unbinned, galaxy = galaxy.upper(),
                redshift = z, title=title, save=saveTo)
#           plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
#                "%splots/%s_field_%s.png" % (wav_range_dir, plot, wav_range), 
#                bbox_inches="tight")


# ------------=========== Save and display plot =============----------
       
        if plots:
            plt.show()















##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    wav_range="4200-"
    galaxy = "ngc3557"
#    galaxy = "ngc7075"
#    galaxy = "ic1459"
    discard = 2 # rows of pixels to discard- must have been the same 
            #    for all routines 
    vLimit = 2 #

    plot_results(galaxy, discard=discard, vLimit=vLimit, 
        wav_range=wav_range, plots=True, nointerp = True)




