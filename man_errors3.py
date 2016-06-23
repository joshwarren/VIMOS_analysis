## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20151015 Routine to manipulate the output from the 
## galamdring monte carlo simulations into a useful form.
## warrenj 20160114 Combined man_errors and man_errors2 with use of 
## keywords.
## warrenj 20160215 Removed options as they were incorrect methods
## and not in use any longer.
## warrenj 20160322 Use man_errors2 for when there are no rep
## performed by errors.py

## *************************** KEYWORDS ************************* ##
# Galaxy
# wav_range
## ************************************************************** ##

import numpy as np # for reading files
import array
import sys # for early exit of python
import os

def man_errors(galaxy, wav_range=""):
    dir = "/Data/vimosindi/analysis/%s/gas_MC/" % (galaxy)
    output_dir = "/Data/vimosindi/analysis/%s/results/%s/" % (galaxy, wav_range)
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

    componants = ["stellar1", "gas", "stellar2"] #,"SF_gas", "shock_g"] 
    dynamics = [v, s, h3, h4]
    dynamics_uncert = [v_s, s_s, h3_s, h4_s]

# ------------====== Reading Files and calc std_devs ==========----------
    for i in range(len(componants)):
        for bin in range(n_bins):
# Bestfit values
            glamdring_file = dir + str(bin) + ".dat"
            vel, sig, h3s, h4s = np.loadtxt(glamdring_file, unpack=True)
            v[bin] = vel[i]
            s[bin] = sig[i]
            h3[bin] = h3s[i]
            h4[bin] = h4s[i]


        v_s =v
        s_s=s
        h3_s=h3
        h4_s=h4



        if not os.path.exists(output_dir):
            os.makedirs(output_dir) 
            
        v_file = output_dir + "gal_" + componants[i] + "_vel.dat"
        s_file = output_dir + "gal_" + componants[i] + "_sigma.dat"
        h3_file = output_dir + "gal_" + componants[i] + "_h3.dat"
        h4_file = output_dir + "gal_" + componants[i] + "_h4.dat"
            



        f_v = open(v_file, 'w')
        f_s = open(s_file, 'w')
        f_h3 = open(h3_file, 'w')
        f_h4 = open(h4_file, 'w')
        for i in range(len(v)):
            f_v.write(str(v[i]) + '    ' + str(v_s[i]) + '\n')
            f_s.write(str(s[i]) + '    ' + str(s_s[i]) + '\n')
            f_h3.write(str(h3[i]) + '    ' + str(h3_s[i]) + '\n')
            f_h4.write(str(h4[i]) + '    ' + str(h4_s[i]) + '\n')

    print "Done"


# ------------============= Collecting Chi2/DOF =============----------
    chi2_file = output_dir + "chi2.dat"
    f_c2 = open(chi2_file, 'w')
    for bin in range(n_bins):
        chi2_bin_file = dir + "chi2/%s.dat" % (str(bin))
        f_c2.write(str(np.loadtxt(chi2_bin_file)) + '\n')

    

    

##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[0]
    wav_range = "4200-"

    man_errors(galaxy, wav_range=wav_range)
