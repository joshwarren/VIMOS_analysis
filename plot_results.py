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
## warrenj 20160726 Now plots in a more object orientated way and creates a
## grid of plots too. This supersedes plot_results2.py

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
#from cap_plot_velfield import plot_velfield #as plot_velfield
import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from plot_histogram import plot_histogram
import ppxf_util as util
from numpy.polynomial import legendre
import os
import colormaps as cm
from sauron_colormap import sauron
from checkcomp import checkcomp
cc = checkcomp()

# Give axes a saveTo property
plt.axes.saveTo = property(lambda self:str())
# Give axes an x and y on figure grid property
plt.axes.figx = property(lambda self:int())
plt.axes.figy = property(lambda self:int())
# give axes a property to hold a colorbar axes
plt.axes.cax = property(lambda self:plt.axes())
# give axes a property to hold 2 additional axes for showing other axis
plt.axes.ax2 = property(lambda self:plt.axes())
plt.axes.ax3 = property(lambda self:plt.axes())

vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
vin_dir_cube = '%s/Data/vimos/cubes' % (cc.base_dir)
ain_dir = '%s/Data/alma' % (cc.base_dir)
out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)


#-----------------------------------------------------------------------------
def set_lims(galaxy, v, vLimit, plot_species, plot_type, mean_centered=False,
	positive=False, symmetric=False):
    if 'uncert' in plot_type:
        uncert = True
        positive = True
    else:
        uncert = False

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
        mean_centered = False
        symmetric = True
    if 'sigma' in plot_type:
        plot_type = 'sigma'
    if 'h3' in plot_type:
        plot_type = 'h3'
        symmetric = True
    if 'h4' in plot_type:
        plot_type = 'h4'
    if 'image' in plot_type or 'Flux' in plot_type:
        plot_type = 'image'
        postive = True
    if 'Equivalent Width' in plot_type:
        plot_type = 'equivalent_width'
        positive = True
    if 'line ratio' in plot_type and 'OIII' in plot_type:
        plot_type = 'lrOIII'
        positive = True
    elif 'line ratio' in plot_type and 'Hbeta' in plot_type:
        plot_type = 'lrHbeta'
        positive = True

    # Read Limits file
    f = '%s/%s/limits.dat' % (vin_dir, galaxy)
    species, types = np.loadtxt(f, unpack=True, dtype=str, usecols=(0,1),
                                skiprows=1)
    if uncert:
        c = (4,5)
    else:
        c = (2,3)
    mins, maxs = np.genfromtxt(f, unpack=True, usecols=c, skip_header=1, missing_values='nan', filling_values=np.nan)
    
    s = np.where(species == plot_species)[0]
    n = s[np.where(types[s] == plot_type)[0]]

    v_sorted = np.array(sorted(np.unique(v)))
    if len(v_sorted) < 2*vLimit:
        v_sorted = np.array(sorted(v))
    # remove nan from v_sorted
    v_sorted = v_sorted[~np.isnan(v_sorted)]

    # If limits nan in limits file.
    if np.isnan(mins[n]): mins[n] = v_sorted[vLimit]
    if np.isnan(maxs[n]): maxs[n] = v_sorted[-vLimit-1]
    
    # If plot limits not in limits file
    if len(n) == 0: # or len(mins[n]) == 0:
        vmin, vmax = v_sorted[vLimit], v_sorted[-vLimit-1]
    else:
        # min must be less than max!
        if mins[n] < maxs[n]:
            vmin, vmax = mins[n], maxs[n]
        else:
        	raise ValueError("min must be less than max limit.")
        	vmin, vmax = maxs[n], mins[n]

    # center the range on the mean.
    if mean_centered:
        mean = np.mean(v)
        d = 2*mean - vmax
        vmin = mean-d

    # Uncertinty should always positive.
    if positive and vmin < 0:
        vmin = 0

	# Find next entry above (for min) and below (for max)
         # may throw exception if nan in v_sorted, but handles fine.
    vmin = v_sorted[min(np.where(v_sorted >= vmin)[0])]
    vmax = v_sorted[max(np.where(v_sorted <= vmax)[0])]

    # Velocity plots should be symmetric about 0
    if symmetric:
         if abs(vmin)>vmax:
             # Check if both are < or > 0.
             vmin=-abs(vmax)
         else:
             vmax=abs(vmin)

    return vmin, vmax
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#def use_templates(galaxy, glamdring=False):
#    if cc.device == 'glamdring':
#        template_weighting = '/users/warrenj/analysis/' + galaxy + \
#	    '/templates.txt' 
#    else:
#        template_eighting = '%s/%s/templates.txt' % (vin_dir, galaxy)
#
#    templatesToUse = np.loadtxt(template_weighting, usecols=(0,), dtype='i')
#    return templatesToUse
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
def add_CO(ax, galaxy, header, close=False):
    CO_image_dir="%s/%s-mom0.fits" % (ain_dir, galaxy)
    if os.path.exists(CO_image_dir):
        CO_image, CO_header = pyfits.getdata(CO_image_dir, 0, header=True)

        #remove random extra dimenisons.
        CO_image = np.sum(np.sum(CO_image,axis=0), axis=0)

        CO_x = np.arange(CO_header['NAXIS1'])*CO_header['CDELT1']*60*60
        CO_y = np.arange(CO_header['NAXIS2'])*CO_header['CDELT2']*60*60

        #x += max(ax.get_xlim())
        #y -= max(ax.get_ylim())

        # Coordinates of VIMOS pointing
        vc = SkyCoord(header['HIERARCH CCD1 ESO INS IFU RA'], 
            header['HIERARCH CCD1 ESO INS IFU DEC'], 
            unit=(u.deg, u.deg))

        # Coordinates of ALMA pointing
        ac = SkyCoord(CO_header['CRVAL1'], CO_header['CRVAL2'],
            unit=(u.deg, u.deg))

        # Offset between the two pointings
        CO_x -= ((vc.ra.degree - header['CRPIX1']*header['CDELT1']/(60*60)) -
            (ac.ra.degree +
            CO_header['CRPIX1']*CO_header['CDELT1']/(60*60)))*60*60
                
        CO_y += ((vc.dec.degree - header['CRPIX2']*header['CDELT2']/(60*60)) -
            (ac.dec.degree +
            CO_header['CRPIX2']*CO_header['CDELT2']/(60*60)))*60*60
            
        cs = ax.contour(CO_x,CO_y,CO_image, colors='k')

        saveTo = os.path.dirname(ax.saveTo)+"/withCO/" + \
            os.path.basename(ax.saveTo)
        if not os.path.exists(os.path.dirname(saveTo)):
            os.makedirs(os.path.dirname(saveTo))
        plt.savefig(saveTo, bbox_inches="tight")

        if close:
            plt.close()
        else:
            # Make lines thinner for pdf by finding the line objects
            for o in ax.get_children():
                if type(o) is LineCollection:
                    o.set_linewidth(0.3)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def plot_results(galaxy, discard=0, wav_range="", vLimit=2, norm="lwv", 
    plots=False, nointerp=False, residual=False, CO=False, show_bin_num=False,
    **kwargs):    

    data_file =  "%s/galaxies.txt" % (vin_dir)
    # different data types need to be read separetly
    z_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
        usecols=(1,4,5))
    x_gals, y_gals = x_gals.astype(int), y_gals.astype(int)
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    z = z_gals[i_gal]


    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""

    tessellation_File = "%s/%s/voronoi_2d_binning_output.txt" % (vin_dir, galaxy)
    tessellation_File2 = "%s/%s/voronoi_2d_binning_output2.txt" %(vin_dir, galaxy)
    dataCubeDirectory = "%s/%s.cube.combined.fits" % (vin_dir_cube, galaxy)
    output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
    out_plots = "%splots" % (output)
    out_nointerp = "%s/notinterpolated" % (out_plots)
    vin_dir_gasMC = "%s/%s/gas_MC" % (vin_dir, galaxy)

    # lists the files produced by man_errors[2].py
    outputs = glob.glob(output+'gal_*.dat')
    #outputs = []

    # Create figure and array for axes
    if any('OIII' in o for o in outputs):
        n_rows = len(outputs)/2 +1
    else:
        n_rows = len(outputs)/2
    f = plt.figure(frameon=False)
    ax_array = []

    


    # Read tessellation file
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1 )
    x, y, bin_num = x.astype(int), y.astype(int), bin_num.astype(int)
    n_spaxels = len(bin_num)
    number_of_bins = int(max(bin_num)+1)
    order = bin_num.argsort()


    center_bin = bin_num[x_gals[i_gal]*(max(y)+1) + y_gals[i_gal]]
# ----------========= Reading the spectrum  =============--------- 
    galaxy_data, header = pyfits.getdata(dataCubeDirectory, 0, header=True)

    #galaxy_data = np.rot90(galaxy_data,2)
    s = galaxy_data.shape
    rows_to_remove = range(discard)
    rows_to_remove.extend([s[1]-1-i for i in range(discard)])
    cols_to_remove = range(discard)
    cols_to_remove.extend([s[2]-1-i for i in range(discard)])

    galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
    galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)



    galaxy_data_unbinned = np.sum(galaxy_data, axis=0)
    #galaxy_data_unbinned = galaxy_data_unbinned.flatten()

# ----------============ Spatially binning ==============----------
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

    #flux_bar_binned = flux_bar_binned/np.median(flux_bar_binned)

# ------------=============== Plot image ================----------
    if CO:
        galaxy_data_unbinned = None

    print "        Image"
    
    title = "Total Flux"
    CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"

    ax = f.add_subplot(111, aspect='equal')
    saveTo = "%s/total_image_%s.png" % (out_nointerp, wav_range)
    ax.saveTo = saveTo
    ax.figx, ax.figy = 0, 0

    fmin, fmax = set_lims(galaxy, flux_bar_binned, vLimit, 'stellar', title)
        
    ax = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
        flux_bar_binned, vmin=fmin, vmax=fmax, nodots=True,
        show_bin_num=show_bin_num, colorbar=True, label=CBLabel,
        title=title, cmap="gist_yarg", ax=ax)
    ax_array.append(ax)
    f.delaxes(ax)
    f.delaxes(ax.cax)
    if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
    if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
    
    if plots:
        plt.show()

# ------------========= Plot intensity (& EW) ===========----------
    print "        gas map(s) and equivalent widths"


    # Getting the gas templates used.
    FWHM_gal = 4*0.71
    degree = 4 # in all our analysis.

    temp_weights = []
    temp_name = []
    bestfit = []
    lam = []
    emission_lines = []
    for i in range(number_of_bins):
        lam_bin = np.loadtxt("%s/lambda/%d.dat" % (vin_dir_gasMC, i),
        	unpack=True)
        lam.append(lam_bin)
        # Getting the contiuum model used in ppxf
        l = np.linspace(-1, 1, len(lam))

        loglam_bin = np.log(lam_bin)
        emission_lines_bin, line_name, line_wav = util.emission_lines(
            loglam_bin, [lam_bin[0],lam_bin[-1]], FWHM_gal, quiet=True)
        emission_lines.append(emission_lines_bin)
        
        weights_dir = "%s/temp_weights/%d.dat" % (vin_dir_gasMC, i)
        temp_weights_temp = np.loadtxt(weights_dir, unpack=True, usecols=(1,))
        temp_name_temp=np.loadtxt(weights_dir, unpack=True, usecols=(0,), dtype=str)
        if '[OIII]5007d' not in temp_name_temp: np.concatenate((temp_name_temp, [0]))
        temp_weights.append(temp_weights_temp)
        temp_name.append(temp_name_temp)


        bestfit_bin = np.loadtxt("%s/bestfit/%d.dat" %(vin_dir_gasMC,i), unpack=True)
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
    for c in components:
        i = np.where(line_name == c)[0][0]
        temp_flux = np.trapz(emission_lines[0][:,i], x=lam[0])
        wav = line_wav[i]
        flux = []
        for i in range(len(temp_weights)):
            j = np.where(temp_name[i] == c)[0][0]
            flux.append(temp_weights[i][j]*temp_flux)

        if 'OIII' in c:
            c_title = '[OIII]'
        elif 'Hbeta' in c:
            c_title = r'H$_\beta$'
        elif 'Hgamma' in c:
            c_title = r'H$_\gamma$'
        else:
            c_title = c

        f_title = "%s Flux" % (c_title)
        fh_title = "%s Flux Histogram" % (c_title)
        # from header
        fCBtitle = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
        f_min, f_max = set_lims(galaxy, flux, vLimit, c, f_title)

        saveTo = "%s/%s_flux_hist_%s.png" % (out_plots, c, wav_range)
        plot_histogram(flux, galaxy=galaxy.upper(), redshift=z,
            vmin=f_min,vmax=f_max, weights=n_spaxels_in_bin, title=fh_title,
            xaxis=fCBtitle, save=saveTo)
        
        ax_y = set_ax_y(c)

        ax = f.add_subplot(111, aspect='equal')
        saveTo = "%s/%s_img_%s.png" % (out_nointerp, c, wav_range)
        ax.saveTo = saveTo
        ax.figx, ax.figy = 0, ax_y
        
        
        ax = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            flux, vmin=f_min, vmax=f_max, colorbar=True, nodots=True,
            label=fCBtitle, title=f_title, cmap = 'gist_yarg', ax=ax)
        ax_array.append(ax)
        f.delaxes(ax)
        f.delaxes(ax.cax)
        if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
        if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
            
        if plots: plt.show()


        
        #equiv_width = flux/continuum[:,np.argmin(np.abs(lam-wav))]
        continuum = np.zeros(number_of_bins)
        for k in range(number_of_bins):
######************* Check continuum calc here ***********############
            continuum[k] = bestfit[k][np.argmin(np.abs(lam[k]-wav))] - \
                np.max(emission_lines[k][i])*temp_weights[k][j]
######################################################################
        equiv_width = flux/continuum#(bestfit[0][np.argmin(np.abs(lam[0]-wav))])#-np.max(emission_lines[:,i])*temp_weights[:,j])
        
        eq_title = "%s Equivalent Width" % (c_title)
        eqh_title = "%s Equivalent Width Histogram" % (c_title)
        eqCBtitle = r"Equivalent Width ($\AA$)"

        eq_min, eq_max = set_lims(galaxy, equiv_width, vLimit, c, eq_title)

        saveTo = "%s/%s_eqWidth_hist_%s.png" % (out_plots, c, wav_range)
        plot_histogram(equiv_width, galaxy=galaxy.upper(), redshift=z,
            vmin=eq_min,vmax=eq_max, weights=n_spaxels_in_bin, title=eqh_title,
            xaxis=eqCBtitle, save=saveTo)
        
        ax = f.add_subplot(111, aspect='equal')
        saveTo = "%s/%s_equiv_width_%s.png" % (out_nointerp, c, wav_range)
        ax.saveTo = saveTo
        ax.figx, ax.figy = 0, ax_y+1


        ax = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            equiv_width, vmin=eq_min, vmax=eq_max, colorbar=True, nodots=True,
            label=eqCBtitle, title=eq_title, ax=ax)
        ax_array.append(ax)
        f.delaxes(ax)
        f.delaxes(ax.cax)
        if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
        if hasattr(ax,'ax3'): f.delaxes(ax.ax3)

# ------------============ Find gas masks! ==============----------
		amp = np.zeros(len(components), number_of_bins)
		for k in range(number_of_bins):
			i = np.where(line_name == c)[0][0]
####************** Not sure about structure of emission_lines ****############
        	amp[k] = max(emission_lines[0][:,i], x=lam[0])














			
  
# ------------============== Read fields ================----------
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

# ------------=========== Setting titles etc ============----------
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

# ------------============ Setting v range ==============----------
        if "vel" in plot_title:
            if norm == "lum":
                v_binned -= v_binned[center_bin]
            if norm == "lwv":
                lwv = v_unbinned*galaxy_data_unbinned
                v_binned -= np.nanmean(lwv)*n_spaxels/np.nansum(
                    galaxy_data_unbinned)
            if norm == "sig":
                sig_file = output+'gal_stellar_sigma.dat'
                s_binned, s_uncert_binned = np.loadtxt(sig_file, unpack=True)
                s_sort = sorted(np.unique(s_binned))
                c = np.where(s_binned > s_sort[-6])
                v_binned -= np.mean(v_binned[c[0]])

        vmin, vmax = set_lims(galaxy, v_binned, vLimit, plot_title, plot_title)

        v_uncert_min, v_uncert_max = set_lims(galaxy, v_uncert_binned, vLimit, 
            plot_title, utitle)

# ------------============== Plot Histogram =============----------
        # Field histogram

        saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title, wav_range)
        plot_histogram(v_binned, galaxy=galaxy.upper(), redshift=z,
            vmin=vmin,vmax=vmax, weights=n_spaxels_in_bin, title=htitle,
            xaxis=CBLabel, save=saveTo)
        # Uncertainty histogram
        saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title+'_uncert', wav_range)
        plot_histogram(v_uncert_binned, galaxy=galaxy.upper(), redshift=z,
            vmin=v_uncert_min,vmax=v_uncert_max, weights=n_spaxels_in_bin,
            title=uhtitle, xaxis=CBLabel, save=saveTo)

        if plots:
            plt.show()            
  
# ------------==== Plot velfield - no interperlation ====----------
        if CO:
            galaxy_data_unbinned_sav = galaxy_data_unbinned
            galaxy_data_unbinned = None
        if nointerp:
            # Field plot
            ax = f.add_subplot(111, aspect='equal')
            saveTo = ("%s/%s_field_%s.png" % (out_nointerp, plot_title, wav_range))
            ax.saveTo = saveTo
            ax.figx, ax.figy = ax_x, ax_y
           
            ax = plot_velfield_nointerp(x, y, bin_num, xBar,
                yBar, v_binned, vmin=vmin, vmax=vmax, #flux_type='notmag',
                nodots=True, show_bin_num=show_bin_num, colorbar=True, 
                label=CBLabel, #flux_unbinned=galaxy_data_unbinned, 
                title=title, cmap=cmap, ax=ax)
            ax_array.append(ax)
            f.delaxes(ax)
            f.delaxes(ax.cax)
            if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
            if hasattr(ax,'ax3'): f.delaxes(ax.ax3)

            # Uncertainty plot
            saveTo = "%s/%s_field_%s.png" % (out_nointerp,plot_title+'_uncert', 
                wav_range)
            ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
                v_uncert_binned, vmin=v_uncert_min, vmax=v_uncert_max,
                flux_type='notmag', nodots=True, show_bin_num=show_bin_num,
                colorbar=True, label=CBLabel, galaxy = galaxy.upper(),
                redshift = z, title=utitle, save=saveTo, close=not CO)#, cmap=cm.blue)
            if CO:
                ax1.saveTo = saveTo
                add_CO(ax1, galaxy, header, close=True)
            

        if plots:
            plt.show()
        if CO:
            galaxy_data_unbinned = galaxy_data_unbinned_sav

# ------------============= Plot residuals ==============----------
    if residual:
        print "        " + residual + " residuals"

        average_residuals = np.zeros(number_of_bins)
        for i in range(number_of_bins):
            bestfit = np.loadtxt('%s/bestfit/%d.dat' % (vin_dir_gasMC, i))
            spectrum = np.loadtxt('%s/input/%d.dat' % (vin_dir_gasMC, i))
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
                
        minres, maxres = set_lims(galaxy, average_residuals, vLimit, 'stellar', 
            'residuals', positive=True) #mean_centered=True,
        
        CBLabel = "Residuals"
        title = str.capitalize(residual) + \
	    " Residuals of Bestfit to Normalised Spectrum"
        saveTo = "%s/%s_residual_%s.png" % (out_nointerp, residual, wav_range)

        ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
            average_residuals, vmin=minres, vmax=maxres, flux_type='notmag',
            nodots=True, show_bin_num=show_bin_num, colorbar=True, 
            label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
            galaxy = galaxy.upper(), redshift = z, title=title, 
            save=saveTo, close=not CO)#, cmap = cm.blue)
        if plots:
            plt.show()
        if CO:
            ax1.saveTo = saveTo
            add_CO(ax1, galaxy, header, close=True)

# ------------=============== Plot Chi2/DOF =============----------
    print "        chi2"

    chi2 = np.zeros(number_of_bins)
    for i in range(number_of_bins):
        chi2[i] = np.loadtxt("%s/chi2/%d.dat" % (vin_dir_gasMC, i))

    minchi2, maxchi2 = set_lims(galaxy, chi2, vLimit, 'stellar', 'chi2',
    	positive = True)
    
    CBLabel = "Chi2/DOF"
    title = "Chi2/DOF of the bestfit"
    saveTo = "%s/chi2_%s.png" % (out_nointerp, wav_range)

    ax1 = plot_velfield_nointerp(x, y, bin_num, xBar, yBar, chi2, 
        vmin=minchi2, vmax=maxchi2, flux_type='notmag',
        nodots=True, show_bin_num=show_bin_num, colorbar=True, 
        label=CBLabel, flux_unbinned=galaxy_data_unbinned, 
        galaxy = galaxy.upper(), redshift = z, title=title, 
        save=saveTo, close=not CO)#, cmap=cm.blue)
    if plots:
        plt.show()
    if CO:
        ax1.saveTo = saveTo
        add_CO(ax1, galaxy, header, close=True)

# ------------============ Line ratio maps ==============----------
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


            # Plotting the inverse of the above
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

            lr_min, lr_max = set_lims(galaxy, line_ratio, vLimit, cA, lr_title)


            ax = f.add_subplot(111, aspect='equal')
            saveTo = "%s/lineratio/%s_%s_line_ratio_%s.png" % (out_nointerp, cB, cA, 
                wav_range)
            ax.saveTo = saveTo
            ax.figx, ax.figy = n, n_rows-1

            ax = plot_velfield_nointerp(x, y, bin_num, xBar, yBar,
                line_ratio, vmin=lr_min, vmax=lr_max, colorbar=True,
                nodots=True, title=lr_title, label=lrCBtitle, ax=ax)

            ax_array.append(ax)
            f.delaxes(ax)
            f.delaxes(ax.cax)
            if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
            if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
    
# ------------============= Plot and save ===============----------

    print "        Plotting and saving"

    for i, a in enumerate(ax_array):
        f.add_axes(a)
        #a.axis('tight')
        f.add_axes(a.cax)
        if hasattr(a,'ax2'): f.add_axes(a.ax2)
        if hasattr(a,'ax2'): f.add_axes(a.ax3)
        if not os.path.exists(os.path.dirname(a.saveTo)):
            os.makedirs(os.path.dirname(a.saveTo))  
        f.savefig(a.saveTo, bbox_inches="tight")

        if CO:
            add_CO(a, galaxy, header)

        f.delaxes(a)
        f.delaxes(a.cax)
        a.change_geometry(n_rows, 3, a.figy*3+a.figx+1)
        if hasattr(a,'ax2'): f.delaxes(a.ax2)        
        if hasattr(a,'ax2'): f.delaxes(a.ax3)



    for a in ax_array:
        f.add_axes(a)
        f.add_axes(a.cax)
        a.xaxis.set_visible(False)
        a.yaxis.set_visible(False)
        a.axis('off')
        a.autoscale(False)

    f.set_size_inches(8.5,n_rows*1.8)
    f.tight_layout(h_pad=0.5)#pad=0.4, w_pad=0.5, h_pad=1.0)
    f.subplots_adjust(top=0.94)
    f.suptitle(galaxy.upper())

    saveTo = "%s/grid_%s.pdf" % (out_plots, wav_range)
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



