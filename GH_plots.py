## ==================================================================
## Produce Gauss-Hermite plots
## ==================================================================
## warrenj 20150825 Routine to plot h_3 and h_4 vs v/sigma for all
## bins
## warrenj 20150917 Loop added to plot both without editing plot


import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting
from scipy.stats import gaussian_kde # for calc plot density
#-----------------------------------------------------------------------------


def GH_plots(galaxy, wav_range="", plots=False):

    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""
        wav_range = ""
    limit=1 # how many of the lowest/highest errors to clip
    scale=1 #scale the errors


    tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output2.txt"
    output_v = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_vel.dat" % (wav_range_dir)
    output_v_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_vel_uncert.dat" % (wav_range_dir)
    output_sigma = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_sigma.dat" % (wav_range_dir)
    output_sigma_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_sigma_uncert.dat" % (wav_range_dir)
    output_h3 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h3.dat" % (wav_range_dir)
    output_h3_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h3_uncert.dat" % (wav_range_dir)
    output_h4 = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h4.dat" % (wav_range_dir)
    output_h4_uncert = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_h4_uncert.dat" % (wav_range_dir)





    vel = np.loadtxt(output_v)
    sigma = np.loadtxt(output_sigma)
    h3 = np.loadtxt(output_h3)
    h4 = np.loadtxt(output_h4)

    vel_uncert = np.loadtxt(output_v_uncert, usecols=(1,))
    sigma_uncert = np.loadtxt(output_sigma_uncert, usecols=(1,))
    h3_uncert = np.loadtxt(output_h3_uncert, usecols=(1,))
    h4_uncert = np.loadtxt(output_h4_uncert, usecols=(1,))



    vel += -np.mean(vel)
    x = vel/sigma

    x_uncert = np.sqrt((vel_uncert/vel)**2 + (sigma_uncert/sigma)**2)*x/scale
    x_uncert_sorted = sorted(x_uncert)
    x_uncert_min = x_uncert_sorted[limit]
    x_uncert_max = x_uncert_sorted[-limit-1]

    plot = "h3"

    for i in [0,1]:
 
    
        if plot == "h3":
            y = h3
            y_uncert = h3_uncert/scale
            ytitle = r"$h_3$"
        else:
            y = h4
            y_uncert = h4_uncert/scale
            ytitle = r"$h_4$"




    # Calculate the point density
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    
    
    
    #plt.plot(x, y, 'bx')
        plt.title("Local " + ytitle + " - (v/sigma) relation")
        plt.xlabel(r"$v/\sigma$")
        plt.ylabel(ytitle)
        plt.errorbar(x, y, xerr=np.clip(x_uncert,x_uncert_min,x_uncert_max), 
            yerr=y_uncert, fmt='+',zorder=1)
        plt.scatter(x, y, c=z, s=50, edgecolor='',zorder=2)
        ax=plt.gca()
        plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), 
            verticalalignment='top',
            transform=ax.transAxes)
    #plt.savefig("/home/warrenj/Desktop/" + plot + "-(v-sigma)_" + wav_range + \
    #".png", bbox_inches="tight")
        plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
            "%s/plots/%s-(v-sigma)_%s.png" % (wav_range_dir, plot, wav_range), \
            bbox_inches="tight")
        if plots:
            plt.show()
    
    
        plot = "h4"






##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    galaxy = "ngc3557"
    wav_range="4200-"

    GH_plots(galaxy, wav_range, plots=True)

