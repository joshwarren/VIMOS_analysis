## ==================================================================
## 		Load and Pickle data
## ==================================================================
## warrenj 20160810 Routine to load pPXF fits from Glamdring into 
##		data object and pickle it, ready for use by other routines.
## warrenj 20161217 Added section to replace man_errors.py routine.
##
## *************************** KEYWORDS ************************* ##
## opt 		'kin'	Option of kinamatics (kin) or stellar population
##					(pop).
## ************************************************************** ##

import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
import os
import cPickle as pickle
import warnings
from Bin import Data, emission_line
from checkcomp import checkcomp
cc = checkcomp()


vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
vin_dir_cube = '%s/Data/vimos/cubes' % (cc.base_dir)
out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)



#-----------------------------------------------------------------------------
def pickler(galaxy, discard=0, wav_range="", norm="lwv", opt="kin",	**kwargs):
	print "    Loading D"

	tessellation_File = "%s/%s/voronoi_2d_binning_output_%s.txt" % (vin_dir, galaxy, opt)
	tessellation_File2 = "%s/%s/voronoi_2d_binning_output2_%s.txt" %(vin_dir, galaxy, opt)
	dataCubeDirectory = "%s/%s.cube.combined.corr.fits" % (vin_dir_cube, galaxy)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range)
	if opt == "kin":
		vin_dir_gasMC = "%s/%s/gas_MC" % (vin_dir, galaxy)
	elif opt == "pop":
		vin_dir_gasMC = "%s/%s/pop_MC" % (vin_dir, galaxy)
	out_pickle = '%s/pickled' % (output)

	# Check tessellation file is older than pPXF outputs (checks vin_dir_gasMC/0.dat only).
	if os.path.getmtime(tessellation_File) > os.path.getmtime('%s/0.dat' % (vin_dir_gasMC)
		): 
		bin_num = np.loadtxt(tessellation_File, unpack=True, skiprows = 1, usecols=(2,), 
			dtype=int)
		if os.path.exists('%s/%i.dat' % (vin_dir_gasMC,max(bin_num))) and not \
			os.path.exists('%s/%i.dat' % (vin_dir_gasMC, max(bin_num)+1)):
			# Issue warning, but do not stop.
			warnings.warn('WANING: The tesselation file '+\
				'voronoi_2d_binning_output_%s.txt may to have been changed.' % (opt))
		else:
			# Stop and raise exception
			raise UserWarning('WANING: The tesselation file '+\
				'voronoi_2d_binning_output_%s.txt has been overwritten.' % (opt))
# ------------======== Reading the spectrum  ============----------
	D = Data(np.loadtxt(tessellation_File, unpack=True, skiprows = 1, 
			usecols=(0,1,2)))

	galaxy_data, header = pyfits.getdata(dataCubeDirectory, 0, header=True)
	
	s = galaxy_data.shape
	rows_to_remove = range(discard)
	rows_to_remove.extend([s[1]-1-i for i in range(discard)])
	cols_to_remove = range(discard)
	cols_to_remove.extend([s[2]-1-i for i in range(discard)])
	
	galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
	galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
	
	D.unbinned_flux = np.nansum(galaxy_data, axis=0)
	
	FWHM_gal = 4*0.71
	temp_wav = np.loadtxt('%s/models/miles_library/m0001V' % (cc.home_dir),
		usecols=(0,), unpack=True)
	for i in range(D.number_of_bins):
		D.bin[i].spectrum = np.loadtxt("%s/input/%d.dat" % (vin_dir_gasMC,i), 
			unpack=True)
		D.bin[i].noise = np.loadtxt("%s/noise_input/%d.dat" % 
			(vin_dir_gasMC,i), unpack=True)
		D.bin[i].lam = np.loadtxt("%s/lambda/%d.dat" % (vin_dir_gasMC, i))
		
		# Getting emission templates used
		lines = {'Hdelta':4101.76, 'Hgamma':4340.47, 'Hbeta':4861.33, 
			'[OIII]5007d':5004.0,'[NI]d':5200.0}
		matrix = np.loadtxt("%s/bestfit/matrix/%d.dat" % (vin_dir_gasMC, i),dtype=str)
		ms = matrix.shape
		for j in range(ms[0]):
			if not matrix[j,0].isdigit():
				D.bin[i].components[matrix[j,0]] = emission_line(D.bin[i],
					matrix[j,0],lines[matrix[j,0]],matrix[j,1:].astype(float))
				D.add_e_line(matrix[j,0],lines[matrix[j,0]])
		
		#Setting the weighting given to the gas templates 
		temp_name, temp_weight = np.loadtxt("%s/temp_weights/%d.dat" % (
			vin_dir_gasMC, i), unpack=True, dtype='str')
		temp_weight = temp_weight.astype(float)
		D.bin[i].set_templates(temp_name, temp_weight)

		D.bin[i].bestfit = np.loadtxt("%s/bestfit/%d.dat" %(vin_dir_gasMC,i), 
			unpack=True)
		if opt == 'kin':
			D.bin[i].apweight = np.loadtxt("%s/apweights/%d.dat" %(vin_dir_gasMC,i), 
				unpack=True)
		elif opt == 'pop':
			D.bin[i].mpweight = np.loadtxt("%s/mpweights/%d.dat" %(vin_dir_gasMC,i), 
				unpack=True)
	D.xBar, D.yBar = np.loadtxt(tessellation_File2, unpack=True, skiprows = 1)
# ------------=========== Read kinematics results ==============----------
	componants = [d for d in os.listdir(vin_dir_gasMC + "/gas") if \
		os.path.isdir(os.path.join(vin_dir_gasMC + "/gas", d))]
	if len(componants) == 0: gas =0
	elif 'gas' in componants: gas = 1
	elif 'Shocks' in componants and 'SF' in componants: gas = 2
	else: gas = 3

	

	for c in D.list_components:
		dynamics = {'vel':np.zeros(D.number_of_bins)*np.nan, 
			'sigma':np.zeros(D.number_of_bins)*np.nan, 'h3':np.zeros(D.number_of_bins)*np.nan, 
			'h4':np.zeros(D.number_of_bins)*np.nan}
		dynamics_uncert = {'vel':np.zeros(D.number_of_bins)*np.nan, 
			'sigma':np.zeros(D.number_of_bins)*np.nan, 'h3':np.zeros(D.number_of_bins)*np.nan, 
			'h4':np.zeros(D.number_of_bins)*np.nan}

		for bin in range(D.number_of_bins):
			# Bestfit values
			glamdring_file = "%s/%i.dat" % (vin_dir_gasMC, bin)
			c_in_bin = np.loadtxt(glamdring_file, unpack=True, usecols=(0,), 
				dtype=str)
			# vel, sig, h3s, h4s = np.loadtxt(glamdring_file, unpack=True, 
			# 	usecols=(1,2,3,4))

			if gas == 1 and c != 'stellar':
				c_type = 'gas'
			elif gas == 2 and c != 'stellar':
				if 'H' in c:
					c_type = 'SF'
				else:
					c_type = 'Shocks'
			else:
				c_type = c

			# check componant is in bin
			if c_type in c_in_bin:
				i = np.where(c_in_bin == c_type)[0][0]
				with open(glamdring_file, 'r') as g:
					rows = g.read().splitlines()
					row = np.array(rows[i].split('   '))
				for j, d in enumerate(['vel', 'sigma', 'h3', 'h4']):
					try:
						dynamics[d][bin] = float(row[j+1])
					except IndexError:
						pass
					# 	dynamics[d][bin] = 0
				# dynamics['vel'][bin] = vel[i]				
				# dynamics['sigma'][bin] = sig[i]
				# dynamics['h3'][bin] = h3s[i]
				# dynamics['h4'][bin] = h4s[i]

				# Calculating uncertainties
				if c_type != "stellar": MC_dir = "%s/gas" % (vin_dir_gasMC) 
				else: MC_dir=vin_dir_gasMC

				glamdring_file = '%s/%s/%i.dat' % (MC_dir, c_type, bin)
				# Ignore empty file warnings
				with warnings.catch_warnings():
					warnings.simplefilter("ignore")
					f = np.loadtxt(glamdring_file, unpack=True)
				# Check if MC has been run.
				if len(f) != 0:
					vel, sig, h3s, h4s = f[0,:], f[1,:], f[2,:], f[3,:]
					dynamics_uncert['vel'][bin] = np.std(vel)
					dynamics_uncert['sigma'][bin] = np.std(sig)
					dynamics_uncert['h3'][bin] = np.std(h3s)
					dynamics_uncert['h4'][bin] = np.std(h4s)
				else:
					dynamics_uncert = dynamics

		for kine in dynamics:
			if np.isnan(dynamics[kine]).all():
				D.components[c].unset(kine)
			else:
				if c in D.list_components:
					D.components[c].setkin(kine, dynamics[kine])
					D.components[c].setkin_uncert(kine, dynamics_uncert[kine])

	D.find_restFrame()
# ------------================ Pickling =================----------

	print "    Pickling D"
	if not os.path.exists(out_pickle):
		os.makedirs(out_pickle) 
	if opt == 'kin':
		pickleFile = open("%s/dataObj_%s.pkl" % (out_pickle, wav_range), 'wb')
	elif opt == 'pop':
		pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'wb')
	pickle.dump(D,pickleFile)
	pickleFile.close()
# ------------======== Save flux for KINEMTRY (IDL) =====----------
	if opt == 'kin':
		print "    Saving flux for KINEMETRY (IDL)"
		with open('%s/flux.dat' % (output), 'wb') as f:
			for i in range(D.number_of_bins):
				f.write(str(D.flux[i]) + '\n')

	return D
##############################################################################

# Use of pickler.py

if __name__ == '__main__':

	galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
	galaxy = galaxies[6]

	wav_range="4200-"
	discard = 2 # rows of pixels to discard- must have been the same 
			#	for all routines 

	pickler(galaxy, discard=discard, wav_range=wav_range, opt='pop')