## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.
## warrenj 20151216 Added section to plot residuals.
## warrenj 20160111 Add section to plot histgrams of the fields.

## *************************** KEYWORDS ************************* ##
# galaxy 		Name of the galaxy being plotted: used to find 
#				correct files and to print onto the plot.
# discard	0	Interger giving the number of rows and columns 
#				to be removed from the plot to remove edge 
#				effects.
# wav_range 	null	Imposed wavelength range on top of the automated 
#				limits.	
# vLimit 	2 	Integer giving the number of lowest and highest 
#				results in the plot to be discarded. Defualt 
#				ignores 2 highest and 2 lowest bins.
# norm		"lwv"	Normalisation methods for velocity fields:
#				lwv: luminosity weighted mean of the whole 
#				field is set to 0.
#				lum: velocity of the brightest spaxel is set 
#				to 0.
# plots 	False   Boolean to show plots as routine runs.
# nointerp 	False 	Boolean to use interpolation between bins in 
#				plots or not.
# residual 	False	Method to measure the residuals:
#			mean: use the mean of the residuals in each 
#				bin.
#			median: use the median of the residuals in 
#				each bin.
#			max: use the maximum of the residuals in 
#				each bin.
#			False: do not calculate and produce plot of 
#				residuals.
## ************************************************************** ##


from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for array handling
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from plot_histogram import plot_histogram



#-----------------------------------------------------------------------------

def plot_results(galaxy, discard=0, wav_range="", vLimit=2, norm="lwv", 
    plots=False, nointerp=False, residual=False, **kwargs):

    data_file =  "/Data/vimosindi/analysis/galaxies.txt"
    # different data types need to be read separetly
    z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,4,5))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    z = z_gals[i_gal]



    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""

    tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output2.txt"

    output = "/Data/vimosindi/analysis/%s/results/%s" % (galaxy,wav_range_dir)

#    outputs = glob.glob(output+'gal_*.dat')
    outputs = glob.glob(output+'gal_stellar_vel*.dat')


# Read tessellation file
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1)
#    x = x[::-1]
#    xBin = xBin[::-1]
#    y = y[::-1]
#    yBin = yBin[::-1]
    n_spaxels = len(bin_num)
    number_of_bins = int(max(bin_num)+1)
    order = bin_num.argsort()


    center_bin = bin_num[x_gals[i_gal]*(max(y)+1) + y_gals[i_gal]]

# ------------========== Total flux per bin ===========----------
# ----------========= Reading the spectrum  =============---------

# FILE_SEARCH returns an array even in cases where it only returns
# one result. This is NOT equivalent to a scalar. 
    dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
        "*_cube.fits" % (galaxy)) 

## Directory for plotting single quadrant
#dataCubeDirectory = glob.glob("/Data/vimosindi/%s-3/Q2/calibrated/cube/" \
#    "*_fluxcal_cube.fits" % (galaxy)) 

    galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)

#    galaxy_data = np.rot90(galaxy_data,2)
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
#    xBar = yBar1[::-1]
#    yBar = xBar1[::-1]
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
        plot_title = plot.split('gal_')[-1].split('.')[0]
        print "       ", plot_title
        plt.close('all')
        v_binned, v_uncert_binned = np.loadtxt(plot, unpack=True)


# Asign v to every spaxel in bin
        v_unbinned = np.zeros(galaxy_data_unbinned.shape)
        v_uncert_unbinned = np.zeros(galaxy_data_unbinned.shape)
        for spaxel in range(n_spaxels):
            v_unbinned[x[spaxel],y[spaxel]] = v_binned[bin_num[spaxel]]
            v_uncert_unbinned[x[spaxel],y[spaxel]] = \
                v_uncert_binned[bin_num[spaxel]]
# ------------============ Setting v range =============----------
        if "vel" in plot_title:
            if norm == "lum":
                v_binned -= v_binned[center_bin]
            if norm == "lwv":
                lwv = v_unbinned*galaxy_data_unbinned
                v_binned -= np.nanmean(lwv)*n_spaxels/np.nansum(galaxy_data_unbinned)

# Limits on field
        vmax = max(v_binned)
        vmin = min(v_binned)
        v_sorted = sorted(np.unique(v_binned))
        if len(v_sorted) < 2*vLimit:
            v_sorted = sorted(v_binned)
        vmin = v_sorted[vLimit]
        vmax = v_sorted[-vLimit-1]
# Make velocity fields symmetric
        if "vel" in plot_title:
            if abs(vmin)>vmax:
                if vmin*vmax>0:
                    vmin, vmax = vmax, vmin
                else:
                    vmin=-abs(vmax)
            else:
                vmax=abs(vmin)


# Limits on uncertainty field
        v_uncert_max = max(v_uncert_binned)
        v_uncert_min = min(v_uncert_binned)
        v_uncert_sorted = sorted(np.unique(v_uncert_binned))
        if len(v_uncert_sorted) < 2*vLimit:
            v_uncert_sorted = sorted(v_uncert_binned)
        v_uncert_min = v_uncert_sorted[vLimit]
        v_uncert_max = v_uncert_sorted[-vLimit-1]
        mean_w_uncert = np.mean(v_uncert_binned)
        d_uncert = v_uncert_max-mean_w_uncert
        v_uncert_min = mean_w_uncert-d_uncert
        if v_uncert_min < 0:
            v_uncert_min=0

# ------------============= Plot velfield ==============----------
        CBLabel = None
        if "vel" in plot_title:
            title = 'Velocity'
            CBLabel = "LOSV (km s$^{-1}$)"
#	    vmin=-35
#	    vmax=-vmin
#        vmax=375
#        vmin=240
#        v_binned += 90
#        v_uncert_max = 100
        if "sigma" in plot_title:
            title = 'Velocity Dispersion'
            CBLabel = "LOSVD (km s$^{-1}$)"
        if "h3" in plot_title: title = 'h3'
        if "h4" in plot_title: title = 'h4'



        if "stellar" in plot_title:
            utitle = "Stellar Uncertainty " + title + " Map"
            htitle = "Stellar " + title + " Histogram"
            uhtitle = "Stellar Uncertainty " + title + " Histogram"
            title = "Stellar " + title + " Map"
        else:
            utitle = "Hot Gas Uncertainty " + title + " Map"
            htitle = "Hot Gas " + title + " Histogram"
            uhtitle = "Hot Gas Uncertainty " + title + " Histogram"
            title = "Hot Gas " + title + " Map"

  
# ------------================= Plot Histogram ===============----------
# Field histogram

        saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
            "%splots/%s_hist_%s.png" % (wav_range_dir, plot_title, wav_range)
        plot_histogram(v_binned, galaxy=galaxy.upper(), redshift=z, vmin=vmin,vmax=vmax, weights=n_spaxels_in_bin, title=htitle, xaxis=CBLabel, save=saveTo)
# Uncertainty histogram
        saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
            "%splots/%s_hist_%s.png" % (wav_range_dir, plot_title+'_uncert', 
            wav_range)
        plot_histogram(v_uncert_binned, galaxy=galaxy.upper(), redshift=z, vmin=v_uncert_min,vmax=v_uncert_max, weights=n_spaxels_in_bin, title=uhtitle, xaxis=CBLabel, save=saveTo)

        if plots:
            plt.show()

# ------------===== Plot velfield - no interperlation ======----------
        if nointerp:
# Field plot
            saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
                plot_title, wav_range)
            plot_velfield_nointerp(x, y, bin_num, xBar, yBar, v_binned, 
                vmin=vmin, vmax=vmax, 
                nodots=True, show_bin_num=True, colorbar=True, 
                label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
                galaxy = galaxy.upper(), redshift = z, title=title, 
                save=saveTo)
# Uncertainty plot
            saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
                plot_title+'_uncert', wav_range)
            plot_velfield_nointerp(x, y, bin_num, xBar, yBar, v_uncert_binned, 
                vmin=v_uncert_min, vmax=v_uncert_max, 
                nodots=True, show_bin_num=True, colorbar=True, 
                label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
                galaxy = galaxy.upper(), redshift = z, title=utitle, 
                save=saveTo)

# ------------===== Plot velfield - with interperlation ====----------
        else:
# Field plot
           saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/%s_field_%s.png" % (wav_range_dir, plot_title, 
                wav_range)
           plot_velfield(xBar, yBar, v_binned, vmin=vmin, vmax=vmax, 
                nodots=False, colorbar=True, label=CBLabel, 
                flux_unbinned=galaxy_data_unbinned, galaxy = galaxy.upper(),
                redshift = z, title=title, save=saveTo)
# Uncertainty plot
           saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
                "%splots/%s_field_%s.png" % (wav_range_dir, 
                plot_title+'_uncert', wav_range)
           plot_velfield(xBar, yBar, v_uncert_binned, vmin=v_uncert_min, 
                vmax=v_uncert_max, nodots=False, colorbar=True, label=CBLabel, 
                flux_unbinned=galaxy_data_unbinned, galaxy = galaxy.upper(),
                redshift = z, title=utitle, save=saveTo)

# ------------=========== Save and display plot =============----------
       
        if plots:
            plt.show()


# ------------================= Plot residuals ===============----------

    if residual:
        print "        " + residual + " residuals"
        bestfit_dir = "/Data/vimosindi/analysis/%s/gas_MC/" % (galaxy) +\
            "bestfit/"
        data_dir = "/Data/vimosindi/analysis/%s/gas_MC/" % (galaxy) +\
            "input/"

        average_residuals = np.zeros(number_of_bins)

        for i in range(number_of_bins):
            bestfit = np.loadtxt(bestfit_dir +'%d.dat' % (i))
            spectrum = np.loadtxt(data_dir +'%d.dat' % (i))
            residuals = np.abs(spectrum - bestfit)
# remove edge pixels
            residuals = np.delete(residuals, [np.arange(5), 
	        len(residuals)+np.arange(-5,0)], axis=0)

            if residual=="mean":
                average_residuals[i] = np.mean(residuals)
            elif residual=="median":
                average_residuals[i] = np.median(residuals)
            elif residual=="max":
                average_residuals[i] = np.max(np.abs(residuals))

        res_sorted = sorted(np.unique(average_residuals))
        maxres = res_sorted[-vLimit-1]
#        minres = res_sorted[vLimit]

        mean_w = np.mean(average_residuals)
        d = maxres-mean_w
        minres = mean_w-d
        if minres < 0:
            minres=0

        CBLabel = "Residuals"
        title = str.capitalize(residual) + \
	    " Residuals of Bestfit to Normalised Spectrum"
        saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
            "%splots/notinterpolated/%s_residual_%s.png" % (wav_range_dir, 
            residual, wav_range)

        plot_velfield_nointerp(x, y, bin_num, xBar, yBar, average_residuals, 
            vmin=minres, vmax=maxres, 
            nodots=True, show_bin_num=True, colorbar=True, 
            label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
            galaxy = galaxy.upper(), redshift = z, title=title, 
            save=saveTo)
        if plots:
            plt.show()

# ------------================= Plot Chi2/DOF ===============----------
    print "        chi2"
    chi2_dir = "/Data/vimosindi/analysis/%s/gas_MC/chi2/" % (galaxy)
    chi2 = np.zeros(number_of_bins)
    for i in range(number_of_bins):
        chi2[i] = np.loadtxt("%s%d.dat" % (chi2_dir, i))

    chi2_sorted = sorted(np.unique(chi2))
    maxchi2 = chi2_sorted[-vLimit-1]
    minchi2 = 0 # chi2_sorted[vLimit]
    
    CBLabel = "Chi2/DOF"
    title = "Chi2/DOF of the bestfit"
    saveTo = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%splots/notinterpolated/chi2_%s.png" % (wav_range_dir, wav_range)

    plot_velfield_nointerp(x, y, bin_num, xBar, yBar, chi2, 
        vmin=minchi2, vmax=maxchi2, 
        nodots=True, show_bin_num=False, colorbar=True, 
        label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
        galaxy = galaxy.upper(), redshift = z, title=title, 
        save=saveTo)
    if plots:
        plt.show()

##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[5]

    wav_range="4200-"
    discard = 2 # rows of pixels to discard- must have been the same 
            #    for all routines 
    vLimit = 2 #

    plot_results(galaxy, discard=discard, vLimit=vLimit, 
        wav_range=wav_range, plots=False, nointerp = True, residual = "median")




