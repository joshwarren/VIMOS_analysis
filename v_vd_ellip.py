## ==================================================================
## Plotting v/sig vs ellip
## ==================================================================
## warrenj 20160224 Routine to plot v/sig vs ellipticity.

import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import matplotlib.pyplot as plt # used for plotting


def v_vd_ellip(wav_range=""):
    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""
        wav_range = ""
    discard =2
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    data_file2 =  "/Data/vimosindi/analysis/galaxies2.txt"

    ellip = np.loadtxt(data_file2, unpack=True, skiprows=1, usecols=(1,))


    

    v_vd = np.zeros(len(galaxies))
    for i in range(len(galaxies)):
        galaxy = galaxies[i]
        tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
            "voronoi_2d_binning_output_kin.txt"
        dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
            "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % (galaxy))
        v_file = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
            "%sgal_stellar_vel.dat" % (wav_range_dir)
        s_file = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
            "%sgal_stellar_sigma.dat" % (wav_range_dir)


# ------------========== Reading the data cube ===========----------
        galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)
        s = galaxy_data.shape
        rows_to_remove = range(discard)
        rows_to_remove.extend([s[1]-1-i for i in range(discard)])
        cols_to_remove = range(discard)
        cols_to_remove.extend([s[2]-1-i for i in range(discard)])

        galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
        galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)


        galaxy_data_error = pyfits.getdata(dataCubeDirectory[0], 1)
        galaxy_data_error = np.delete(galaxy_data_error, rows_to_remove, axis=1)
        galaxy_data_error = np.delete(galaxy_data_error, cols_to_remove, axis=2)

# Collapse spectrum
        galaxy_data = np.sum(galaxy_data, axis=0)
        galaxy_data_error = np.sum(galaxy_data_error, axis=0)
        ##galaxy_data_error += galaxy_data


# ******** handling of error ends here ****** (so far) *******
# binning flux
        x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
            skiprows = 1)
        n_spaxels = len(bin_num)
        number_of_bins = int(max(bin_num)+1)
        
        flux_binned = np.zeros((number_of_bins))
        for spaxel in range(n_spaxels):
            flux_binned[int(bin_num[spaxel])] += galaxy_data[y[spaxel],x[spaxel]]
#            n_spaxels_in_bin[int(bin_num[spaxel])] += 1

#        for bin in range(number_of_bins):
#            flux_bar_binned[bin] = flux_bar_binned[bin]/n_spaxels_in_bin[bin]

        flux_binned = flux_binned/np.median(flux_binned)
        
        
# ------------======== Reading the velocity field ========----------
# Read tessellation file

        vel, vel_uncert = np.loadtxt(v_file, unpack=True)
        vel -= np.median(vel)

# ------------======== Reading the velocity field ========----------
        sig, sig_uncert = np.loadtxt(s_file, unpack=True)


        v_vd[i] = np.sum(flux_binned*np.power(vel,2))/np.sum(flux_binned*np.power(sig,2))

    v_vd = np.sqrt(v_vd)

    plt.scatter(ellip,v_vd)
    plt.show()

         












##############################################################################

# Use of v_vd-ellip.py

if __name__ == '__main__':
    v_vd_ellip(wav_range="4200-")
