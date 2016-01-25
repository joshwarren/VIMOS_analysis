## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20151015 Routine to manipulate the output from the 
## galamdring monte carlo simulations into a useful form.
## warrenj 20160114 Combined man_errors and man_errors2 with use of 
## keywords.

## ************************************************************** ##
## NB: only vel and vel uncertainty are affected by different methods.
## ************************************************************** ##


## *************************** KEYWORDS ************************* ##
# vel_method	median	Method of calculating velocity for bin:
#			median: use median of MC
#			mean: use mean (min and max value discarded)
#			None: bestfit results (no MC)
# err_method	median	Method of calculating uncertainty for bin:
#			median: use median of MC
#			mean: use mean (min and max value discarded)
#			std_vel: standard deviation of velocity in MC
#			None: bestfit results (no MC)
# bias		True	Boolean for if the bias keyword was set. NB: if
#			False then err_method must be std_vel (default) 
#			or None.
## ************************************************************** ##

import numpy as np # for reading files
import array
import sys # for early exit of python

def man_errors():
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[5]
    wav_range = "4200-/"

    dir = "/Data/vimosindi/analysis/%s/gas_MC/" % (galaxy)
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

    componants = [stellar, gas] #,SF_gas, shock_g] 
    dynamics = [v, s, h3, h4]
    dynamics_uncert = [v_s, s_s, h3_s, h4_s]

# ------------====== Reading Files and calc std_devs ==========----------
    for i in range(len(componants)):
        for bin in range(n_bins):

# Bestfit values
            glamdring_file = dir + str(bin) + ".dat"
            vel, sig, h3s, h4s = np.loadtxt(glamdring_file, unpack=True)
            v[bin] = vel[i]
            sig[bin] = sig[i]
            h3[bin] = h3s[i]
            h4[bin] = h4s[i]

# Calculating uncertainties
            glamdring_file = dir + componants[i] + "/" + str(bin) + ".dat"
            vel, sig, h3s, h4s =np.loadtxt(glamdring_file, unpack=True)

# Error thrown if all 500 reps have not been completed for this bin        
            if len(vel) != 5000:
                print("Not all reps completed")
                sys.exit()

            v_s[bin] = np.std(vel)
            s_s[bin] = np.std(sig)
            h3_s[bin] = np.std(h3s)
            h4_s[bin] = np.std(h4s)

        output_dir = "/Data/vimosindi/analysis/%s/results/%s" % (galaxy, 
            wav_range)
#        if componant == "stellar":
#            v_file = output_dir + "gal_vel.dat"
#            s_file = output_dir + "gal_sigma.dat"
#            h3_file = output_dir + "gal_h3.dat"
#            h4_file = output_dir + "gal_h4.dat"
#        else:
#            v_file = output_dir + "gal_" + componants[i] + "_vel.dat"
#            s_file = output_dir + "gal_" + componants[i] + "_sigma.dat"
#            h3_file = output_dir + "gal_" + componants[i] + "_h3.dat"
#            h4_file = output_dir + "gal_" + componants[i] + "_h4.dat"

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



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    man_errors()
