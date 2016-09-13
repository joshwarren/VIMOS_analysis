# ==================================================================
# Rebinning using Voronoi tessellation
# ==================================================================
# warrenj 20150515 Process to rebin spaxels together in order to
# create a minimum S/N ratio. targetSN is approx 10-15 for just v and
# sigma, and 40-50 for h3 and h4. 
# warrenj 20160913 Ported to python

import numpy as np
import glob
from astropy.io import fits
import ppxf_util as util
from voronoi_2d_binning import voronoi_2d_binning
from checkcomp import checkcomp
cc=checkcomp()


# ----------===============================================---------
# ----------======== Check overwrite of target SN =========---------
# ----------===============================================---------
def check_overwrite(new, old):
	if new != old:
		A = input('Are you sure you want to overwrite the old target ' + 
			'of %d with a new target of %d? (Y/N) ' % (old, new))
	if A == "N" or A == "n": new = old
	return new




def binning_spaxels(galaxy, discard=2, targetSN=None):
# ----------===============================================---------
# ----------============ Default parameters ===============---------
# ----------===============================================---------

	dir = "%s/Data/vimos" % (cc.base_dir)
	data_file = "%s/analysis/galaxies.txt" %(dir)
	galaxy_gals = np.loadtxt(data_file, usecols=(0,), unpack=True, 
		dtype=str, skiprows=1)
	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_used_gals = np.loadtxt(
		data_file, skiprows=1, usecols=(1,2,3,4,5,6), 
		dtype='float,float,float,int,int,float')
	try:
		i_gal = np.where(galaxy_gals == galaxy)[0]
	except:
		i_gal = -1
		galaxy_gals = [galaxy_gals, galaxy]


	if targetSN is None and i_gal != -1:
		targetSN=SN_used_gals[i_gal]
	elif i_gal  != -1: 
		targetSN = check_overwrite(targetSN, SN_used_gals[i_gal])
	elif targetSN is not None:
		SN_used_gals = [SN_used_gals, targetSN]
	else:
		targetSN = 30.0
		SN_used_gals = [SN_used_gals, targetSN]




# ----------================= Save SN_used ===============---------
	np.savetxt(data_file, np.column_stack([galaxy_gals, z_gals, vel_gals[i],
		sig_gals[i], x_gals[i], y_gals[i], SN_used_gals]), 
		fmt='%s %10.9f %10.9f %10.5f %2i %2i %2.2f', 
		header ="Galaxy	  z	 velocity	 velocity dispersion	x	 y	 Target SN")
	# with open(data_file, 'w') as f:
	# 	f.write("Galaxy	  z	 velocity	 velocity dispersion	x	 y	 Target SN")
	# 	for i in range(len(galaxy_gals)):
	# 		f.write(galaxy_gals[i] + '  ' + z_gals[i] + '  ' + 
	# 			vel_gals[i] + '  ' + sig_gals[i] + '  ' + x_gals[i] + 
	# 			'  ' + y_gals[i] + '  ' + SN_used_gals[i] + '\n')


# ----------================ Find S/N ================------------
# Final wildcard notes that depending on the method used the quadrants
#may or may not have been flux calibrated. 
	dataCubeDirectory = glob.glob(dir+"cubes/%s.cube.combined.fits" % (galaxy)) 
		
	galaxy_data, header = fits.getdata(dataCubeDirectory[0], 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory[0], 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory[0], 3)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CDELT3']
	s = galaxy_data.shape

	rows_to_remove = range(discard)
	rows_to_remove.extend([s[1]-1-i for i in range(discard)])
	cols_to_remove = range(discard)
	cols_to_remove.extend([s[2]-1-i for i in range(discard)])

	galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
	galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
	galaxy_noise = np.delete(galaxy_noise, rows_to_remove, axis=1)
	galaxy_noise = np.delete(galaxy_noise, cols_to_remove, axis=2)
	galaxy_badpix = np.delete(galaxy_noise, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_noise, cols_to_remove, axis=2)

	n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])

	s = galaxy_data.shape


	signal = np.zeros(n_spaxels)
	noise = np.zeros(n_spaxels)
	x = np.zeros(n_spaxels)
	y = np.zeros(n_spaxels)




# collapsing the spectrum for each spaxel.
	for i in range(s[1]):
		for j if range(s[2]):
			gap=12
			ignore = int((5581 - CRVAL_spec)/CDELT_spec) + np.arange(
				-gap+1,gap)  
			ignore2 = int((5199 - CRVAL_spec)/CDELT_spec) + np.arange(
				-gap+1,gap) 

			## h is the spectrum with the peak enclosed by 'ignore' removed.
			h = np.delete(bin_lin_temp, ignore)
			h = np.delete(h,ignore2)


			half = s[0]/2
			a = np.delete(h,np.arange(-4,0)+len(h),None)/np.median(
				h[np.nonzero(h)]) - h[4:]/np.median(h[np.nonzero(h)])
			a = np.where(np.isfinite(a), a, 0)

			if any(np.abs(a[:0.5*half]) > 0.2):
				lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
			else: 
				lower_limit = -1
			
			# lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
			if any(np.abs(a[1.5*half:]) > 0.2):
				upper_limit = min(np.where(np.abs(a[1.5*half:]) > 0.2)[0])+int(1.5*half)
			else:
				upper_limit = -1
				
			if upper_limit > ignore2[0]: upper_limit+=gap 
			if upper_limit > ignore[0]: upper_limit+=gap

			if lower_limit < 0:
				lower_limit = min(np.where(a[:half] != 0)[0]) + 5
				if lower_limit < 0: lower_limit = 0

			else:
				lower_limit +=5

			if upper_limit > s[0]-1 or upper_limit < half:
				upper_limit = s[0]-6 
			else:
				upper_limit += -5

			signal[i*s[1] + j] = np.mean(
				galaxy_data[lower_limit:upper_limit,i,j])
			noise[i*s[1] + j] = np.mean(
				galaxy_noise[lower_limit:upper_limit,i,j])

			# Assign x and y
			x[*s[1]+j] = i
			y[i*s[2]+j] = j

			a = where(galaxy_badpix[:,i,j] == 1, count)
			if count != 0:
				signal[i*s[1] + j] = 0
				noise[i*s[1] + j] = 0


	binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, signal, noise, targetSN, plot=True, quiet=True)




	order = sorted(binNum)
	xBin = np.zeros(n_spaxels)
	yBin = np.zeros(n_spaxels)

	# spaxel number
	i = 0
	for bin in range(max(binNum)+1):
		while i < n_spaxels and bin == binNum[order[i]]:
			xBin[order[i]] = xBar[bin]
			yBin[order[i]] = yBar[bin]
			# move onto next spaxel
			i = i + 1



# Save to a text file the initial coordinates of each pixel together
# with the corresponding bin number computed by this procedure.
# binNum uniquely specifies the bins and for this reason it is the only
# number required for any subsequent calculation on the bins.
	if not os.path.exists("%s/analysis/%s" % (dir,galaxy)):
		os.makedirs("%s/analysis/%s" % (dir, galaxy))

	np.savetxt("%s/analysis/%s/voronoi_2d_binning_output.txt" % (dir,galaxy), 
		np.column_stack([x, y, binNum, xBin, yBin]), 
		fmt='%2i %2i %3i %2.5f %2.5f'
		header=' X"    Y"    BIN_NUM     XBIN   YBIN')

	np.savetxt("%s/analysis/%s/voronoi_2d_binning_output2.txt" % (dir,galaxy), 
		np.column_stack([xBar, yBar]), 
		fmt='%2.5f %2.5f', header='XBAR     YBAR')

