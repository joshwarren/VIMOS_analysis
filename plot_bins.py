## ==================================================================
## 		Plot the histograms of indervidual bins
## ==================================================================
## warrenj 20150112 Routine to plot histograms of MC output of 
## indervidual bins.

from plot_histogram import plot_histogram
import numpy as np # for reading files
from checkcomp import checkcomp
cc = checkcomp()

def plot_bins():
    
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[8]
    wav_range = "4200-/"

    glamdring_dir = "%s/Data/vimosindi/analysis/%s/errors_results/" % (
        cc.base_dir, galaxy)
    tessellation_File = "%s/Data/vimosindi/analysis/%s/" %(cc.base_dir, 
        galaxy) + "voronoi_2d_binning_output.txt"
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1) 
    n_bins = int(max(bin_num)+1)


    for bin in range(n_bins):
        print(bin)
        glamdring_file = glamdring_dir + str(bin) + ".dat"
        glamdring_file2 = glamdring_dir + "/errors/" + str(bin) + ".dat"

        vel, sig, h3s, h4s  =np.loadtxt(glamdring_file, unpack=True)
        vel_e, sig_e, h3s_e, h4s_e  =np.loadtxt(glamdring_file2, unpack=True)

        saveTo = "%s/Data/vimosindi/analysis/%s/" % (base_dir, galaxy) + \
            "errors_results/errors/histograms/vel_%s.png" % (str(bin))

        h = vel_e
        limit = 1
        v_sorted = np.argsort(h)
	plot_histogram(h, galaxy=galaxy.upper(), title='Error' + str(bin), 
        save=saveTo, vmin=h[v_sorted[limit]], vmax=h[v_sorted[-limit-1]])





##############################################################################

# Use of plot_bins.py

if __name__ == '__main__':

    plot_bins()



