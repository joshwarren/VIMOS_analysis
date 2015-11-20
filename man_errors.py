## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 

import numpy as np # for reading files
import array
import sys # for early exit of python

def man_errors():
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[2]
    wav_range = "4200-/"

    glamdring_dir = "/Data/vimosindi/analysis/%s/montecarlo/" % (
        galaxy)
    tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1) 
    n_bins = int(max(bin_num)+1)

    v = np.zeros(n_bins)
    v_s = np.zeros(n_bins)
    s = np.zeros(n_bins)
    s_s = np.zeros(n_bins)
    h3 = np.zeros(n_bins)
    h3_s = np.zeros(n_bins)
    h4 = np.zeros(n_bins)
    h4_s = np.zeros(n_bins)

# ------------============= Reading files ==============----------
#********************************************************************#
#******************CURRENTLY MISSING FINAL BIN **********************#
#****NEED TO REMOVE -1 FROM RANGE(N_BINS) COMMAND WHEN CORRECTED ****#
#********************************************************************#
    for bin in range(n_bins-1):
        glamdring_file = glamdring_dir + str(bin) + ".dat"
        vel, sig, h3s, h4s  =np.loadtxt(glamdring_file, unpack=True)

# Error thrown if all 500 reps have not been completed for this bin        
        if len(vel) != 5000:
            print("Not all reps completed")
            sys.exit()




        v[bin] = np.mean(vel)
        v_s[bin] = np.std(vel)
        s[bin] = np.mean(sig)
        s_s[bin] = np.std(sig)
        h3[bin] = np.mean(h3s)
        h3_s[bin] = np.std(h3s)
        h4[bin] = np.mean(h4s)
        h4_s[bin] = np.std(h4s)


    

    v_file = "/Data/vimosindi/analysis/%s/results/%sgal_vel_uncert.dat" % (
        galaxy, wav_range)
    s_file = "/Data/vimosindi/analysis/%s/results/%sgal_sigma_uncert.dat" % (
        galaxy, wav_range)
    h3_file = "/Data/vimosindi/analysis/%s/results/%sgal_h3_uncert.dat" % (
        galaxy, wav_range)
    h4_file = "/Data/vimosindi/analysis/%s/results/%sgal_h4_uncert.dat" % (
        galaxy, wav_range)


    f_v = open(v_file, 'w')
    f_s = open(s_file, 'w')
    f_h3 = open(h3_file, 'w')
    f_h4 = open(h4_file, 'w')
    for i in range(len(v)):
        f_v.write(str(v[i]) + '    ' + str(v_s[i]) + '\n')
        f_s.write(str(s[i]) + '    ' + str(s_s[i]) + '\n')
        f_h3.write(str(h3[i]) + '    ' + str(h3_s[i]) + '\n')
        f_h4.write(str(h4[i]) + '    ' + str(h4_s[i]) + '\n')




##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    man_errors()
