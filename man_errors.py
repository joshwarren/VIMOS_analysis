## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20151015 Routine to manipulate the output from the 
## galamdring monte carlo simulations into a useful form.
## warrenj 20160114 Combined man_errors and man_errors2 with use of 
## keywords.
## warrenj 20160215 Removed options as they were incorrect methods
## and not in use any longer.

## *************************** KEYWORDS ************************* ##
# Galaxy
# wav_range
## ************************************************************** ##

import numpy as np # for reading files
import os
from checkcomp import checkcomp
cc = checkcomp()

def man_errors(galaxy, wav_range="4200-"):
	dir = "%s/Data/vimos/analysis/%s/gas_MC/" % (cc.base_dir, galaxy)
	output_dir = "%s/Data/vimos/analysis/%s/results/%s/" % (cc.base_dir, 
		galaxy, wav_range)
	tessellation_File = "%s/Data/vimos/analysis/%s/" %(cc.base_dir, galaxy) +\
		"voronoi_2d_binning_output.txt"
	x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
		skiprows = 1) 
	n_bins = int(max(bin_num)+1)

	componants = ["stellar"]

	componants.extend([d for d in os.listdir(dir + "gas/") if \
			   os.path.isdir(os.path.join(dir + "gas/", d))])
	componants = np.array(componants)
	
#    componants.append("stellar")#, "gas"] #,"gas_SF", "gas_Shocks"]

	dynamics = [v, s, h3, h4]
	dynamics_uncert = [v_s, s_s, h3_s, h4_s]

# ------------====== Reading Files and calc std_devs ==========----------
	for i, componant in enumerate(componants):
		v = np.zeros(n_bins)*np.nan
		v_s = np.zeros(n_bins)*np.nan
		s = np.zeros(n_bins)*np.nan
		s_s = np.zeros(n_bins)*np.nan
		h3 = np.zeros(n_bins)*np.nan
		h3_s = np.zeros(n_bins)*np.nan
		h4 = np.zeros(n_bins)*np.nan
		h4_s = np.zeros(n_bins)*np.nan
		for bin in range(n_bins):
			# Bestfit values
			glamdring_file = dir + str(bin) + ".dat"
			c_in_bin = np.loadtxt(glamdring_file, unpack=True, usecols=(0,), 
				dtype=str)
			vel, sig, h3s, h4s = np.loadtxt(glamdring_file, unpack=True, 
				usecols=(1,2,3,4))

			# check componant is in bin
			if componant in c_in_bin:
				v[bin] = vel[i]
				s[bin] = sig[i]
				h3[bin] = h3s[i]
				h4[bin] = h4s[i]

				# Calculating uncertainties
				ex_dir = ""
				if componant != "stellar": ex_dir = "gas/"
				glamdring_file = dir + ex_dir + componant + "/" + \
					str(bin) + ".dat"
				vel, sig, h3s, h4s =np.loadtxt(glamdring_file, unpack=True)

				# Error thrown if 5000 reps not completed in this bin        
				if len(vel) != 5000:
					raise ValueError("Not all reps completed")

				v_s[bin] = np.std(vel)
				s_s[bin] = np.std(sig)
				h3_s[bin] = np.std(h3s)
				h4_s[bin] = np.std(h4s)

			else:
				v[bin] = np.nan
				s[bin] = np.nan
				h3[bin] = np.nan
				h4[bin] = np.nan

				v_s[bin] = np.nan
				s_s[bin] = np.nan
				h3_s[bin] = np.nan
				h4_s[bin] = np.nan

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
	galaxy = galaxies[6]
	wav_range = "4200-"

	man_errors(galaxy, wav_range=wav_range)
