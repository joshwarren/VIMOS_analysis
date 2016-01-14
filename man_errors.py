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
    vel_method = 'median'
    err_method = 'median'
    bias = True

    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[8]
    wav_range = "4200-/"

# Bias not set
    glamdring_dir = "/Data/vimosindi/analysis/%s/montecarlo/" % (
        galaxy)
# Bias set
    glamdring_dir2 = "/Data/vimosindi/analysis/%s/errors_results/" % (
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
    for bin in range(n_bins):
	if bias:
            glamdring_file = glamdring_dir2 + str(bin) + ".dat"
            glamdring_file2 = glamdring_dir2 + "/errors/" + str(bin) + ".dat"
        else:
            glamdring_file = glamdring_dir + str(bin) + ".dat"
	    glamdring_file2 = glamdring_dir + "/errors/" + str(bin) + ".dat"
            if err_method is not None:
		err_method = 'std_vel'



        vel, sig, h3s, h4s =np.loadtxt(glamdring_file, unpack=True)
        vel_e, sig_e, h3s_e, h4s_e =np.loadtxt(glamdring_file2, unpack=True)

# Error thrown if all 500 reps have not been completed for this bin        
        if len(vel) != 5000:
            print("Not all reps completed")
            sys.exit()

	if vel_method !='median' and err_method != 'median':
            v_sorted = np.argsort(vel)
            vel = np.clip(vel, vel[v_sorted[1]], vel[v_sorted[-2]])
	    vel_e = np.clip(vel_e, vel_e[v_sorted[1]], vel_e[v_sorted[-2]])

            v_sorted = np.argsort(vel_e)
            vel = np.clip(vel, vel[v_sorted[1]], vel[v_sorted[-2]])
            vel_e = np.clip(vel_e, vel_e[v_sorted[1]], vel_e[v_sorted[-2]])

            v_sorted = np.argsort(sig)
            sig = np.clip(sig, sig[v_sorted[1]], sig[v_sorted[-2]])
            sig_e = np.clip(sig_e, sig_e[v_sorted[1]], sig_e[v_sorted[-2]])

            v_sorted = np.argsort(sig_e)
            sig = np.clip(sig, sig[v_sorted[1]], sig[v_sorted[-2]])
            sig_e = np.clip(sig_e, sig_e[v_sorted[1]], sig_e[v_sorted[-2]])

            v_sorted = np.argsort(h3s)
            h3s = np.clip(h3s, h3s[v_sorted[1]], h3s[v_sorted[-2]])
	    h3s_e = np.clip(h3s_e, h3s_e[v_sorted[1]], h3s_e[v_sorted[-2]])

            v_sorted = np.argsort(h3s_e)
            h3s = np.clip(h3s, h3s[v_sorted[1]], h3s[v_sorted[-2]])
	    h3s_e = np.clip(h3s_e, h3s_e[v_sorted[1]], h3s_e[v_sorted[-2]])

            v_sorted = np.argsort(h4s)
            h4s = np.clip(h4s, h4s[v_sorted[1]], h4s[v_sorted[-2]])
	    h4s_e = np.clip(h4s_e, h4s_e[v_sorted[1]], h4s_e[v_sorted[-2]])

            v_sorted = np.argsort(h4s_e)
            h4s = np.clip(h4s, h4s[v_sorted[1]], h4s[v_sorted[-2]])
	    h4s_e = np.clip(h4s_e, h4s_e[v_sorted[1]], h4s_e[v_sorted[-2]])

  

	if vel_method == 'mean':
            v[bin] = np.mean(vel)
            s[bin] = np.mean(sig)
            h3[bin] = np.mean(h3s)
            h4[bin] = np.mean(h4s)

	if vel_method == 'median':
            v[bin] = np.median(vel)
            s[bin] = np.median(sig)
            h3[bin] = np.median(h3s)
            h4[bin] = np.median(h4s)

	
	if err_method == 'mean':
            v_s[bin] = np.mean(vel_e)
            s_s[bin] = np.mean(sig_e)
      	    h3_s[bin] = np.mean(h3s_e)
            h4_s[bin] = np.mean(h4s_e)

	if err_method == 'std_vel':
            v_s[bin] = np.std(vel)
            s_s[bin] = np.std(sig)
            h3_s[bin] = np.std(h3s)
            h4_s[bin] = np.std(h4s)

	if err_method == 'median':
            v_s[bin] = np.median(vel_e)
            s_s[bin] = np.median(sig_e)
      	    h3_s[bin] = np.median(h3s_e)
            h4_s[bin] = np.median(h4s_e)

    
    dir = "/Data/vimosindi/analysis/%s/results/4200-/" % (galaxy)
    if vel_method is None:
        v = np.unpack(dir + "no_MC/gal_vel.dat", usecols=(0,), unpack=True)
    if err_method is None:
        v_s = np.unpack(dir + "no_MC/gal_vel.dat", usecols=(1,), unpack=True)




    v_file = dir + "gal_vel.dat"
    s_file = dir + "gal_sigma.dat"
    h3_file = dir + "gal_h3.dat"
    h4_file = dir + "gal_h4.dat"


    f_v = open(v_file, 'w')
    f_s = open(s_file, 'w')
    f_h3 = open(h3_file, 'w')
    f_h4 = open(h4_file, 'w')
    for i in range(len(v)):
        f_v.write(str(v[i]) + '    ' + str(v_s[i]) + '\n')
        f_s.write(str(s[i]) + '    ' + str(s_s[i]) + '\n')
        f_h3.write(str(h3[i]) + '    ' + str(h3_s[i]) + '\n')
        f_h4.write(str(h4[i]) + '    ' + str(h4_s[i]) + '\n')


    print "Bias:       " + str(bias)
    print "vel_method: " + vel_method
    print "err_method: " + err_method
##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    man_errors()
