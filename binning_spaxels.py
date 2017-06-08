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
from errors2 import remove_anomalies, get_dataCubeDirectory
import os


# ----------===============================================---------
# ----------======== Check overwrite of target SN =========---------
# ----------===============================================---------
def check_overwrite(new, old, auto_override=False):
	if not auto_override and new != old:
		A = raw_input('Are you sure you want to overwrite the old target ' + 
			'of %d with a new target of %d? (Y/N) ' % (old, new))
		if A == "N" or A == "n": new = old
	return new




def binning_spaxels(galaxy, discard=2, targetSN=None, opt='kin', auto_override=False):
	print '     Voronoi Binning'
# ----------===============================================---------
# ----------============ Default parameters ===============---------
# ----------===============================================---------
	dir = "%s/Data/vimos" % (cc.base_dir)
	data_file = "%s/analysis/galaxies.txt" %(dir)
	# Check if file has anything in it - it does need to exsist.
	try:
		d = np.loadtxt(data_file, unpack=True, dtype=str)
		galaxy_gals = d[0][1:]
		z_gals, vel_gals, sig_gals = d[1][1:].astype(float), d[2][1:].astype(float), \
			d[3][1:].astype(float),
		x_gals, y_gals = d[4][1:].astype(int), d[5][1:].astype(int)
		SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(6,len(d))}
	except StopIteration:
		galaxy_gals = np.array([])
		z_gals = np.array([])
		vel_gals = np.array([])
		sig_gals = np.array([])
		x_gals = np.array([])
		y_gals = np.array([])
		SN_gals = {}

	try:
		SN_used_gals = SN_gals['SN_%s' % (opt)]
	except KeyError:
		SN_used_gals = np.zeros([len(galaxy_gals)])


	i_gal = np.where(galaxy_gals == galaxy)[0]
	if len(i_gal) == 0:
		i_gal = -1
		galaxy_gals = np.append(galaxy_gals, galaxy)

	if targetSN is None and i_gal != -1:
		targetSN=SN_used_gals[i_gal]
	elif targetSN is not None and i_gal  != -1: 
		targetSN = check_overwrite(float(targetSN), SN_used_gals[i_gal], auto_override)
		SN_used_gals[i_gal] = targetSN
	elif targetSN is not None and i_gal == -1:
		SN_used_gals = np.append(SN_used_gals, targetSN)
	else:
		targetSN = 30.0
		SN_used_gals = np.append(SN_used_gals, targetSN)

	if i_gal == -1:
		z_gals = np.append(z_gals, 0)
		vel_gals = np.append(vel_gals, 0)
		sig_gals = np.append(sig_gals, 0)
		x_gals = np.append(x_gals, 0)
		y_gals = np.append(y_gals, 0)

		SN_gals = {t:np.append(v, 0) for t, v in SN_gals.iteritems()}

# ----------================= Save SN_used ===============---------
	SN_gals['SN_%s' % (opt)] = SN_used_gals

	temp = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}" + ''.join(['{%i:%i}'%(i+6,len(t)+1) 
		for i, t in enumerate(SN_gals.keys())]) + "\n"

	SN_titles = list(SN_gals.keys())
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "z", "velocity", "sigma", "x", "y", 
			*(s for s in SN_titles)))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(z_gals[i],7)), 
				str(round(vel_gals[i],4)), str(round(sig_gals[i],4)), 
				str(int(x_gals[i])), str(int(y_gals[i])), 
				*(str(round(SN_gals[s][i],2)) for s in SN_titles)))

# ----------================ Find S/N ================------------
# Final wildcard notes that depending on the method used the quadrants
#may or may not have been flux calibrated. 
	dataCubeDirectory = get_dataCubeDirectory(galaxy) 

	galaxy_data, header = fits.getdata(dataCubeDirectory, 0, header=True)
	galaxy_noise = fits.getdata(dataCubeDirectory, 1)
	galaxy_badpix = fits.getdata(dataCubeDirectory, 3)

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

	s = galaxy_data.shape

	
	x = np.zeros(s[1]*s[2], dtype=int)
	y = np.zeros(s[1]*s[2], dtype=int)


# collapsing the spectrum for each spaxel.
	signal = np.nanmedian(galaxy_data, axis=0).flatten()
	noise = np.nanmedian(galaxy_noise,axis=0).flatten()

	for i in range(s[1]):
		for j in range(s[2]):
			# Assign x and y
			x[i*s[1]+j] = i
			y[i*s[2]+j] = j

	mask = (np.isfinite(signal)) * (np.isfinite(noise))
	signal = signal[mask]
	noise = noise[mask]
	x = x[mask]
	y = y[mask]
	n_spaxels = np.sum(mask)

	if not os.path.exists("%s/analysis/%s/%s/setup" % (dir, galaxy, opt)):
		os.makedirs("%s/analysis/%s/%s/setup" % (dir, galaxy, opt))

	binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, signal, noise, targetSN, quiet=True, plot=False,
        saveTo='%s/analysis/%s/%s/setup/binning.png' %(dir, galaxy, opt))

	order = np.argsort(binNum)
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

# ------------================ Saving Results ===============---------------			

	temp = "{0:3}{1:3}{2:8}{3:9}{4:9}\n"
	temp2 = "{0:9}{1:9}\n"

	with open("%s/analysis/%s/%s/setup/voronoi_2d_binning_output.txt" % (dir,galaxy,opt), 
		'w') as f:
		f.write(temp.format('X"', 'Y"', 'BIN_NUM', 'XBIN', 'YBIN'))
		for i in range(len(xBin)):
			f.write(temp.format(str(int(x[i])), str(int(y[i])), str(int(binNum[i])), 
				str(round(xBin[i],5)), str(round(yBin[i],5))))


	with open("%s/analysis/%s/%s/setupvoronoi_2d_binning_output2.txt" % (dir,galaxy,opt), 
		'w') as f:
		f.write(temp2.format('XBAR','YBAR'))
		for i in range(len(xBar)):
			f.write(temp2.format(str(round(xBar[i],5)), str(round(yBar[i],5)))) 

	print 'Number of bins: ', max(binNum)+1