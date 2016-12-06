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
from errors2 import remove_anomalies
import os


# ----------===============================================---------
# ----------======== Check overwrite of target SN =========---------
# ----------===============================================---------
def check_overwrite(new, old):
	if new != old:
		A = raw_input('Are you sure you want to overwrite the old target ' + 
			'of %d with a new target of %d? (Y/N) ' % (old, new))
		if A == "N" or A == "n": new = old
	return new




def binning_spaxels(galaxy, discard=2, targetSN=None):
	print '     Voronoi Binning'
# ----------===============================================---------
# ----------============ Default parameters ===============---------
# ----------===============================================---------
	dir = "%s/Data/vimos" % (cc.base_dir)
	data_file = "%s/analysis/galaxies.txt" %(dir)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_used_gals = np.loadtxt(
		data_file, skiprows=1, usecols=(1,2,3,4,5,6), unpack=True,
		dtype='float,float,float,int,int,float')
	try:
		i_gal = np.where(galaxy_gals == galaxy)[0]
	except:
		i_gal = -1
		galaxy_gals = [galaxy_gals, galaxy]


	if targetSN is None and i_gal != -1:
		targetSN=SN_used_gals[i_gal]
	elif targetSN is not None and i_gal  != -1: 
		targetSN = check_overwrite(targetSN, SN_used_gals[i_gal])
		SN_used_gals[i_gal] = targetSN
	elif targetSN is not None and i_gal == -1:
		SN_used_gals = [SN_used_gals, targetSN]
	else:
		targetSN = 30.0
		SN_used_gals = [SN_used_gals, targetSN]




# ----------================= Save SN_used ===============---------
	temp = "{0:12}{1:11}{2:9}{3:15}{4:4}{5:4}{6:10}\n"
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "z", "velocity", "vel dispersion", "x", "y", 
			"Target SN"))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(z_gals[i],7)), 
				str(round(vel_gals[i],4)), str(round(sig_gals[i],4)), 
				str(int(x_gals[i])), str(int(y_gals[i])), str(round(SN_used_gals[i],2))))

# ----------================ Find S/N ================------------
# Final wildcard notes that depending on the method used the quadrants
#may or may not have been flux calibrated. 
	dataCubeDirectory = glob.glob("%s/cubes/%s.cube.combined.fits" % (dir,galaxy)) 

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
	galaxy_badpix = np.delete(galaxy_badpix, rows_to_remove, axis=1)
	galaxy_badpix = np.delete(galaxy_badpix, cols_to_remove, axis=2)

	# Check for nan is data set.
	galaxy_badpix[np.isnan(galaxy_data)] = 1

	n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])

	s = galaxy_data.shape


	#signal = np.zeros(n_spaxels)
	#noise = np.zeros(n_spaxels)
	x = np.zeros(n_spaxels)
	y = np.zeros(n_spaxels)

# collapsing the spectrum for each spaxel.
	signal = np.mean(galaxy_data, axis=0).flatten()
	noise = np.sqrt(np.mean(galaxy_noise**2, axis=0)).flatten()


	for i in range(s[1]):
		for j in range(s[2]):
			# Assign x and y
			x[i*s[1]+j] = i
			y[i*s[2]+j] = j

			a = np.where(galaxy_badpix[:,i,j] == 1)[0]
			if len(a) != 0:
				signal[i*s[1] + j] = 0
				noise[i*s[1] + j] = 0.0000000001
				print 'Spaxel containing badpixels: ', i, j
				print 'Number of affected pixels:   ', len(a)

	binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, signal, noise, targetSN, quiet=True, 
        saveTo='%s/analysis/%s/binning.png' %(dir,galaxy))


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

	temp = "{0:3}{1:3}{2:8}{3:9}{4:9}\n"
	temp2 = "{0:9}{1:9}\n"

	with open("%s/analysis/%s/voronoi_2d_binning_output.txt" % (dir,galaxy), 'w') as f:
		f.write(temp.format('X"', 'Y"', 'BIN_NUM', 'XBIN', 'YBIN'))
		for i in range(len(xBin)):
			f.write(temp.format(str(int(x[i])), str(int(y[i])), str(int(binNum[i])), 
				str(round(xBin[i],5)), str(round(yBin[i],5))))


	with open("%s/analysis/%s/voronoi_2d_binning_output2.txt" % (dir,galaxy), 'w') as f:
		f.write(temp2.format('XBAR','YBAR'))
		for i in range(len(xBar)):
			f.write(temp2.format(str(round(xBar[i],5)), str(round(yBar[i],5)))) 

	print 'Number of bins: ', max(binNum)+1