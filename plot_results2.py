## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.
## warrenj 20151216 Added section to plot residuals.
## warrenj 20160111 Add section to plot histgrams of the fields.
## warrenj 20160405 Added keyword CO to overlay CO maps from ALMA if avaible.
## This supersedes plot_results_CO.py

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
#				sig: Noralised to the mean velocity of 5 bins with the
#				highest LOSVD.
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
# CO       False	Boolean to show ALMA CO plots overlaied (if they exist)
## ************************************************************** ##

#import matplotlib # 20160202 JP to stop lack-of X-windows error
#matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from plot_histogram import plot_histogram
import ppxf_util as util
from numpy.polynomial import legendre
import os
import colormaps as cm
from sauron_colormap import sauron

#-----------------------------------------------------------------------------
def set_lims(galaxy, vmin, vmax, plot_species, plot_type):
    if 'stellar' in plot_species:
        plot_species = 'stellar'
    if 'Hbeta' in plot_species:
        plot_species = 'Hbeta'
    if 'Hgamma' in plot_species:
        plot_species = 'Hgamma'
    if 'OIII' in plot_species:
        plot_species = 'OIII'

    if 'vel' in plot_type or 'velocity' in plot_type:
        plot_type = 'velocity'
    if 'sigma' in plot_type:
        plot_type = 'sigma'
    if 'h3' in plot_type:
        plot_type = 'h3'
    if 'h4' in plot_type:
        plot_type = 'h4'
    if 'image' in plot_type or 'flux' in plot_type:
        plot_type = 'image'
    if 'equivalent width' in plot_type:
        plot_type = 'equivalent_width'
    if 'line ratio' in plot_type and 'OIII' in plot_type:
        plot_type = 'lrOIII'
    elif 'line ratio' in plot_type and 'Hbeta' in plot_type:
        plot_type = 'lrHbeta'

    f = '/Data/vimos/analysis/' + galaxy + '/limits.dat'
    species, types = np.loadtxt(f, unpack=True, dtype=str, usecols=(0,1),
                                skiprows=1)
    mins, maxs = np.loadtxt(f, unpack=True, usecols=(2,3), skiprows=1)


    s = np.where(species == plot_species) ### needs more than this to search plot_type too.
    n = np.where(types[s] == plot_type)

    mins = mins[s]
    maxs = maxs[s]

    if np.isnan(mins[n]): mins[n] = vmin
    if np.isnan(maxs[n]): maxs[n] = vmax

    
    if len(n) == 0 or len(mins[n]) == 0:
        return vmin, vmax
    else:
        if mins[n] < maxs[n]:
            return mins[n], maxs[n]
        else:
            print "min must be less than max limit."


    

#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
def use_templates(galaxy, glamdring=False):
    if glamdring:
        template_weighting = '/users/warrenj/analysis/' + galaxy + \
	    '/templates.txt' 
    else:
        template_weighting = '/Data/vimos/analysis/' + galaxy + \
	    '/templates.txt' 

    templatesToUse = np.loadtxt(template_weighting, usecols=(0,), dtype='i')
    return templatesToUse
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
def set_ax_y(plt_title):
    if "gas" in plt_title:
        ax_y=2
    elif "SF" in plt_title:
        ax_y=2
    elif "Shocks" in plt_title:
        ax_y=4
    elif 'Hbeta' in plt_title:
        ax_y=4
    elif 'Hgamma' in plt_title:
        ax_y=6
    elif 'OIII' in plt_title:
        ax_y=2
    elif 'stellar' in plt_title:
       ax_y=0
       
    return ax_y
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------

def plot_results(galaxy, discard=0, wav_range="", vLimit=2, norm="lwv", 
    plots=False, nointerp=False, residual=False, CO=False, show_bin_num=False,
    **kwargs):

    data_file =  "/Data/vimos/analysis/galaxies.txt"
    # different data types need to be read separetly
    z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,4,5))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    z = z_gals[i_gal]



    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""

    tessellation_File = "/Data/vimos/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = "/Data/vimos/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output2.txt"

    output = "/Data/vimos/analysis/%s/results/%s" % (galaxy,wav_range_dir)

    outputs = glob.glob(output+'gal_*.dat')

# Create figure and array of axes
    if any('OIII' in o for o in outputs):
        n_rows = len(outputs)/2 +1
    else:
        n_rows = len(outputs)/2
    f, ax_array = plt.subplots(n_rows, 3, sharex='col', sharey='row')
    f.set_size_inches(8,n_rows*1.8)
    f.tight_layout(h_pad=2.0)#pad=0.4, w_pad=0.5, h_pad=1.0)
    ax_array[1,0].axis('off')
    ax_array[0,0].invert_xaxis()
    f.subplots_adjust(top=0.94)
    f.suptitle(galaxy.upper())


# Read tessellation file
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1)
    n_spaxels = len(bin_num)
    number_of_bins = int(max(bin_num)+1)
    order = bin_num.argsort()


    center_bin = bin_num[x_gals[i_gal]*(max(y)+1) + y_gals[i_gal]]

# ------------========== Total flux per bin ===========----------
# ----------========= Reading the spectrum  =============---------

# FILE_SEARH returns an array even in cases where it only returns
# one result. This is NOT equivalent to a scalar. 
    dataCubeDirectory = glob.glob("/Data/vimos/cubes/%s.cube.combined.fits" \
        % (galaxy)) 

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
    flux_bar_binned = np.zeros((number_of_bins))
    n_spaxels_in_bin = np.zeros((number_of_bins))


    for spaxel in range(n_spaxels):
        flux_bar_binned[int(bin_num[spaxel])] += np.sum(
            galaxy_data[:,y[spaxel],x[spaxel]])
        n_spaxels_in_bin[int(bin_num[spaxel])] += 1

        


    for bin in range(number_of_bins):
        flux_bar_binned[bin] = flux_bar_binned[bin]/n_spaxels_in_bin[bin]

#    flux_bar_binned = flux_bar_binned/np.median(flux_bar_binned)


# ------------============= Read fields =============----------
# Read results files - each entry in array corresponds to a bin (not
# a spaxel)
    for plot in outputs:
        plot_title = plot.split('gal_')[-1].split('.')[0]
        print "       ", plot_title
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
            norm='lwv'
            if norm == "lum":
                v_binned -= v_binned[center_bin]
            if norm == "lwv":
                lwv = v_unbinned*galaxy_data_unbinned
                v_binned -= np.nanmean(lwv)*n_spaxels/np.nansum(
                    galaxy_data_unbinned)
            if norm == "sig":
                sig_file = glob.glob(output+'gal_stellar_sigma*.dat')
                s_binned, s_uncert_binned = np.loadtxt(sig_file[0], unpack=True)
                s_sort = sorted(np.unique(s_binned))
                c = np.where(s_binned > s_sort[-6])
                v_binned -= np.mean(v_binned[c[0]])




                
#            v_binned -=5

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
        im_type = plot_title.split('_')[0]
        if im_type == "gas":
            im_type=""
            ax_y=set_ax_y(plot_title)
        elif im_type == "SF":
            im_type=" (Star Forming)"
            ax_y=set_ax_y(plot_title)
        elif im_type == "Shocks":
            im_type=" (Shocking)"
            ax_y=set_ax_y(plot_title)
        elif 'Hbeta' in im_type:
            im_type=" ("+r'H$_\beta$'+")"
            ax_y=set_ax_y(plot_title)
        elif 'Hgamma' in im_type:
            im_type=" ("+r'H$_\gamma$'+")"
            ax_y=set_ax_y(plot_title)
        elif 'OIII' in im_type:
            im_type=" (OIII)"
            ax_y=set_ax_y(plot_title)
        else:
            im_type=" (" + im_type + ")"
            ax_y=set_ax_y(plot_title)

            
        CBLabel = None
        if "vel" in plot_title:
            ax_x=1
            title = 'Velocity'
            CBLabel = "V (km s$^{-1}$)"
            cmap = sauron

        else:
            cmap = sauron#cm.blue
        if "sigma" in plot_title:
            ax_x=2
            title = 'Velocity Dispersion'
            CBLabel = r'$\mathrm{\sigma}$ (km s$^{-1}$)'

        if "h3" in plot_title:
            ax_x=1
            ax_y+=1
            title = 'h3'

        if "h4" in plot_title:
            ax_x=2
            ax_y+=1
            title = 'h4'        


 



        if "stellar" in plot_title:
            utitle = "Stellar Uncertainty " + title + " Map"
            htitle = "Stellar " + title + " Histogram"
            uhtitle = "Stellar Uncertainty " + title + " Histogram"
            title = "Stellar " + title + " Map"
        else:
            utitle = "Ionised" + im_type + " Gas Uncertainty " + title + " Map"
            htitle = "Ionised" + im_type + " Gas " + title + " Histogram"
            uhtitle = "Ionised" + im_type + " Gas Uncertainty " + title + \
                " Histogram"
            title = "Ionised" + im_type + " Gas\n" + title + " Map"

        if CO:
            galaxy_data_unbinned_sav = galaxy_data_unbinned
            galaxy_data_unbinned = None

        vmin, vmax = set_lims(galaxy, vmin, vmax, plot_title, plot_title)

            
  
# ------------===== Plot velfield - no interperlation ======----------
        if nointerp:
# Field plot
            saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
                plot_title, wav_range)
#            ax_array[ax_y,ax_x] = plot_velfield_nointerp(x, y, bin_num, xBar,
#                yBar, v_binned, vmin=vmin, vmax=vmax, #flux_type='notmag',
#                nodots=True, show_bin_num=show_bin_num, colorbar=True, 
#                label=CBLabel, #flux_unbinned=galaxy_data_unbinned, 
#                title=title, cmap=cmap, save=saveTo,
#                CO=CO, galaxy=galaxy.upper(), redshift=z, header=header)
            plot_velfield_nointerp(x, y, bin_num, xBar,
                yBar, v_binned, vmin=vmin, vmax=vmax, #flux_type='notmag',
                nodots=True, show_bin_num=show_bin_num, colorbar=True, 
                label=CBLabel, #flux_unbinned=galaxy_data_unbinned, 
                title=title, cmap=cmap, save=saveTo,
                CO=CO, galaxy=galaxy.upper(), redshift=z, header=header).set_figure(f)


# Uncertainty plot
            saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/%s_field_%s.png" % (wav_range_dir, 
                plot_title+'_uncert', wav_range)
            ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
                v_uncert_binned, vmin=v_uncert_min, vmax=v_uncert_max,
                flux_type='notmag', nodots=True, show_bin_num=show_bin_num,
                colorbar=True, label=CBLabel, galaxy = galaxy.upper(),
                redshift = z, title=utitle, save=saveTo, CO=CO, header=header)
              #, cmap=cm.blue)

                
                
# ------------=========== Save and display plot =============----------
        if plots:
            plt.show()
        if CO:
            galaxy_data_unbinned = galaxy_data_unbinned_sav
    if CO:
        galaxy_data_unbinned = None
# ------------================== Plot image ==================----------

    print "        Image"
    fmax = max(flux_bar_binned)
    fmin = min(flux_bar_binned)
    f_sorted = sorted(np.unique(flux_bar_binned))
    if len(f_sorted) < 2*vLimit:
        v_sorted = sorted(f_bar_binned)
    fmin = f_sorted[vLimit]
    fmax = f_sorted[-vLimit-1]
    
    title = "Total Flux"
    CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
    saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
        "%splots/notinterpolated/total_image_%s.png" % (wav_range_dir,
        wav_range)

    fmin, fmax = set_lims(galaxy, fmin, fmax, 'stellar', title)
        
    ax_array[0,0] = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
        flux_bar_binned, vmin=fmin, vmax=fmax, nodots=True,
        show_bin_num=show_bin_num, colorbar=True, label=CBLabel,
        title=title, cmap="gist_yarg", save=saveTo, CO=CO, redshift=z,
        galaxy=galaxy.upper(), header=header)
    if plots:
        plt.show()


# ------------================= Plot residuals ===============----------

    if residual:
        print "        " + residual + " residuals"
        bestfit_dir = "/Data/vimos/analysis/%s/gas_MC/" % (galaxy) +\
            "bestfit/"
        data_dir = "/Data/vimos/analysis/%s/gas_MC/" % (galaxy) +\
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
        saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
            "%splots/notinterpolated/%s_residual_%s.png" % (wav_range_dir, 
            residual, wav_range)

        ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            average_residuals, vmin=minres, vmax=maxres, flux_type='notmag',
            nodots=True, show_bin_num=show_bin_num, colorbar=True, 
            label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
            galaxy = galaxy.upper(), redshift = z, title=title, 
            save=saveTo, CO=CO, header=header)#, cmap = cm.blue)
        if plots:
            plt.show()
  
# ------------================= Plot Chi2/DOF ===============----------
    print "        chi2"
    chi2_dir = "/Data/vimos/analysis/%s/gas_MC/chi2/" % (galaxy)
    chi2 = np.zeros(number_of_bins)
    for i in range(number_of_bins):
        chi2[i] = np.loadtxt("%s%d.dat" % (chi2_dir, i))

    chi2_sorted = sorted(np.unique(chi2))
    maxchi2 = chi2_sorted[-vLimit-1]
    minchi2 = 0 # chi2_sorted[vLimit]
    
    CBLabel = "Chi2/DOF"
    title = "Chi2/DOF of the bestfit"
    saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
        "%splots/notinterpolated/chi2_%s.png" % (wav_range_dir, wav_range)

    ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar, chi2, 
        vmin=minchi2, vmax=maxchi2, flux_type='notmag',
        nodots=True, show_bin_num=show_bin_num, colorbar=True, 
        label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
        galaxy = galaxy.upper(), redshift = z, title=title, 
        save=saveTo, CO=CO, header=header)#, cmap=cm.blue)
    if plots:
        plt.show()
        
# ------------============ Plot intensity (& EW) ===============----------
    print "        gas map(s) and equivalent widths"


## Getting the gas templates used.
    FWHM_gal = 4*0.71

    degree = 4 # in all our analysis.
    


    temp_weights = []
    temp_name = []
    bestfit = []
    lam = []
    emission_lines = []
    for i in range(number_of_bins):
        lam_dir = "/Data/vimos/analysis/%s/gas_MC/lambda/%s.dat" % (
            galaxy,str(i))
        lam_bin = np.loadtxt(lam_dir, unpack=True)
        lam.append(lam_bin)
## Getting the contiuum model used in ppxf
        l = np.linspace(-1, 1, len(lam))


        loglam_bin = np.log(lam_bin)
        emission_lines_bin, line_name, line_wav = util.emission_lines(
            loglam_bin, [lam_bin[0],lam_bin[-1]], FWHM_gal, quiet=True)
        emission_lines.append(emission_lines_bin)
        
        weights_dir = "/Data/vimos/analysis/%s/gas_MC/temp_weights/%s.dat"\
            % (galaxy,str(i))
        temp_weights_temp = np.loadtxt(weights_dir, unpack=True, usecols=(1,))
        temp_name_temp=np.loadtxt(weights_dir, unpack=True, usecols=(0,), dtype=str)
        if '[OIII]5007d' not in temp_name_temp: np.concatenate((temp_name_temp, [0]))
        temp_weights.append(temp_weights_temp)
        temp_name.append(temp_name_temp)





        

        bestfit_dir = "/Data/vimos/analysis/%s/gas_MC/bestfit/%s.dat" % (
            galaxy,str(i))
        bestfit_bin = np.loadtxt(bestfit_dir, unpack=True)
        bestfit.append(bestfit_bin)

        
    temp_weights = np.array(temp_weights)
    bestfit = np.array(bestfit)
    lam = np.array(lam)
    emission_lines = np.array(emission_lines)

    

    
    components = []
    temp_name_temp = [item for sublist in temp_name for item in sublist]
    temp_name_lines = [i for i in list(set(temp_name_temp)) if not i.isdigit()]
    check_bins = {i:temp_name_temp.count(i) for i in temp_name_lines}
    
    components = [i for i in temp_name_lines if check_bins[i]==number_of_bins]
    
#    weights_dir = "/Data/vimos/analysis/%s/gas_MC/temp_weights/%s.dat" % (
#        galaxy,str(0))
#    temp_name = np.loadtxt(weights_dir, unpack=True, usecols=(0,),dtype=str)
    
    for c in components:
        i = np.where(line_name == c)[0][0]
        temp_flux = np.trapz(emission_lines[0][:,i], x=lam[0])
        wav = line_wav[i]
        flux = []
        for i in range(len(temp_weights)):
            j = np.where(temp_name[i] == c)[0][0]
            flux.append(temp_weights[i][j]*temp_flux)
        
        f_max = max(flux)
        f_min = min(flux)
        f_sorted = sorted(np.unique(flux))
        f_min = f_sorted[vLimit]
        f_max = f_sorted[-vLimit-1]


#        if 'OIII' in c:
#            f_min = 1
#            f_max = 10
#
#            eq_min = 0.6
#            eq_max = 3.7
#            
#        if 'Hbeta' in c:
#            f_min = 0.5
#            f_max = 3.7
#
#            eq_min = 0.3
#            eq_max = 2.1
#
#        if 'Hgamma' in c:
#            f_min = 0.1
#            f_max = 1.8
#
#            eq_min = 0.1
#            eq_max = 1.9


 

        if 'OIII' in c:
            c_title = '[OIII]'
        elif 'Hbeta' in c:
            c_title = r'H$_\beta$'
        elif 'Hgamma' in c:
            c_title = r'H$_\gamma$'
        else:
            c_title = c

        f_title = "%s Flux" % (c_title)
## from header
        fCBtitle = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
        f_min, f_max = set_lims(galaxy, f_min, f_max, c, f_title)
    
        ax_y = set_ax_y(c)
        saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
        "%splots/notinterpolated/%s_img_%s.png" % (wav_range_dir, c, wav_range)
        
        ax_array[ax_y,0] = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            flux, vmin=f_min, vmax=f_max, colorbar=True, nodots=True,
            label=fCBtitle, title=f_title, cmap = 'gist_yarg', save=saveTo,
            CO=CO, galaxy=galaxy.upper(), redshift=z, header=header)

        if plots: plt.show()

#        equiv_width = flux/continuum[:,np.argmin(np.abs(lam-wav))]
        continuum = np.zeros(number_of_bins)
        for k in range(number_of_bins):
            continuum[k] = bestfit[k][np.argmin(np.abs(lam[k]-wav))] - \
                np.max(emission_lines[k][i])*temp_weights[k][j]
        equiv_width = flux/continuum#(bestfit[0][np.argmin(np.abs(lam[0]-wav))])#-np.max(emission_lines[:,i])*temp_weights[:,j])
        
        eq_max = max(equiv_width)
        eq_min = min(equiv_width)
        eq_sorted = sorted(np.unique(equiv_width))
        eq_min = eq_sorted[vLimit]
        eq_max = eq_sorted[-vLimit-1]
        
        eq_title = "%s Equivalent Width" % (c_title)
        eqCBtitle = r"Equivalent Width ($\AA$)"
        eq_min, eq_max = set_lims(galaxy, eq_min, eq_max, c, eq_title)
    
        saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
            "%splots/notinterpolated/%s_equiv_width_%s.png" % (wav_range_dir,
            c, wav_range)

        ax_array[ax_y+1,0] = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            equiv_width, vmin=eq_min, vmax=eq_max, colorbar=True, nodots=True,
            label=eqCBtitle, title=eq_title, save=saveTo, CO=CO, redshift=z,
            galaxy=galaxy.upper(), header=header)

  
# ------------============== Line ratio maps ==================----------
    if any('OIII' in o for o in outputs):
        print "        line ratios"

        t_num = (len(components)-1)*len(components)/2
        for n in range(t_num):
            i = 0
            m = t_num
            while m > n:
                i += 1
                m -= i

            plotA = len(components)-i-1
            plotB = len(components)-i+n-m

            cA = line_name[plotA]
            cB = line_name[plotB]
        
            iA = np.where(line_name == cA)[0][0]
            temp_fluxA = np.trapz(emission_lines[0][:,iA], x=lam[0])
            iB = np.where(line_name == cB)[0][0]
            temp_fluxB = np.trapz(emission_lines[0][:,iB], x=lam[0])

            fluxA = []
            fluxB = []
            for i in range(len(temp_weights)):
                iA = np.where(temp_name[i] == cA)[0][0]
                fluxA.append(temp_weights[i][iA]*temp_fluxA)
                iB = np.where(temp_name[i] == cB)[0][0]
                fluxB.append(temp_weights[i][iB]*temp_fluxB)

            fluxA=np.array(fluxA)+0.00001
            fluxB=np.array(fluxB)+0.00001


## Plotting the inverse of the above
#            vLimit = 3 
            line_ratio = np.log10(fluxB/fluxA)

            if 'OIII' in cA:
                cA_title = '[OIII]'
            elif 'Hbeta' in cA:
                cA_title = r'H$_\beta$'
            elif 'Hgamma' in cA:
                cA_title = r'H$_\gamma$'
            else:
                cA_title = cA

            if 'OIII' in cB:
                cB_title = '[OIII]'
            elif 'Hbeta' in cB:
                cB_title = r'H$_\beta$'
            elif 'Hgamma' in cB:
                cB_title = r'H$_\gamma$'
            else:
                cB_title = cB
                
            lr_title = "%s/%s Line Ratio" % (cB_title, cA_title)
            lrCBtitle = r"log$_{10}$ (%s/%s)" %(cB_title,cA_title)

            lr_max = max(line_ratio)
            lr_min = min(line_ratio)
            lr_sorted = sorted(np.unique(line_ratio))
            lr_min = lr_sorted[vLimit]
            lr_max = lr_sorted[-vLimit-1]

            lr_min, lr_max = set_lims(galaxy, lr_min, lr_max, cA, lr_title)

#            if 'OIII' in cB and 'Hbeta' in cA:
#                lr_min = -0.5
#                lr_max = 0.75
#            if 'OIII' in cB and 'Hgamma' in cA:
#                lr_min = -0.5
#                lr_max = 2.1
#            if 'Hbeta' in cB and 'Hgamma' in cA:
#                lr_min = -0.5
#                lr_max = 1
            saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
                "%splots/notinterpolated/lineratio/" % (wav_range_dir) + \
                "%s_%s_line_ratio_%s.png" % (cB, cA, wav_range)

            ax_array[-1,n] = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
                line_ratio, vmin=lr_min, vmax=lr_max, colorbar=True,
                nodots=True, title=lr_title, label=lrCBtitle, save=saveTo,
                CO=CO, galaxy=galaxy.upper(), redshift=z, header=header)

    saveTo = "/Data/vimos/analysis/%s/results/" % (galaxy) + \
            "%splots/grid_%s.pdf" % (wav_range_dir, wav_range)

    f.savefig(saveTo, bbox_inches="tight",format='pdf')
    

    
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




