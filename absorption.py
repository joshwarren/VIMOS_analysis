## ==================================================================
## 		Absorption line indices
## ==================================================================
## warrenj 20161218 A routine to measure absorption line indices using Ryan 
##	Hougton's spectools.py
##
##############################################################################
# Keywords:
# line_name: 			Name of line to be found
# lam:					Wavelngth array of observed spectrum and conv_spec if supplied 
# spec:					Observed spectrum
# unc_lam:		None	Wavelength of unconvolved spectrum - from the templates used by pPXF
# unc_spec:		None	Unconcolved spectum
# conv_spec:	None	Convolved spectrum i.e. pPXF bestfit
# noise:		None	Observed noise from reduction pipeline
# lick:			False	Return corrected indices to LICK resolution
##############################################################################
from spectools import *
from tools import length as len
import numpy as np 
from scipy.ndimage.filters import gaussian_filter1d

c = 299792.458 # speed of light in km/s


def absorption(line_name, lam, spec, unc_lam=None, unc_spec=None, conv_spec=None, 
	noise=None, lick=False):
	if line_name=='H_beta' or line_name=='Hbeta':
		line_name='Hb'
	elif line_name=='Mg_b':
		line_name='Mgb'

	if len(unc_lam) != len(unc_spec):
		raise ValueError('Length of unconvoled spectrum and corresponding '+\
			'wavelength should be the same.')
	if conv_spec is not None:
		if len(lam) != len(conv_spec) != len(spec):
			raise ValueError('Length of spectrum (bestfit and observed) and '+\
				'corresponding wavelength should be the same.')
	if noise is not None:
		if len(noise) != len(spec):
			raise ValueError('Length of noise and observed spectrum '+\
				'should be the same.')

	# Catch if index is not in wavelenght range.
	try:
		# Calculate correction for velocity dispersion spreading
		if unc_lam is not None and unc_spec is not None and conv_spec is not None:
			# Bring unc_spec and conv_spec to the same resolution
			sig_pix = np.sqrt(np.abs((unc_lam[1]-unc_lam[0])**2 - (lam[1]-lam[0])**2)
				)/(unc_lam[1]-unc_lam[0])
			if unc_lam[1]-unc_lam[0] > lam[1]-lam[0]:
				unc_spec = gaussian_filter1d(unc_spec, sig_pix)
			else:
				conv_spec = gaussian_filter1d(conv_spec, sig_pix)
			# Line strength of unconvolved (/convolved to LICK resolution) spectrum.
			if lick:
				sig_pix = np.sqrt(200**2 - 64**2)*np.median(lam)/c/(unc_lam[1]-unc_lam[0])
				# unc_spec becomes convolved to 200km/s dispersion (as required for LICK 
				#	indices)
				unc_spec = gaussian_filter1d(unc_spec, sig_pix)
			spectra = spectrum(unc_lam, unc_spec)
			unc_index_value, unc_index_va, _, _ = spectra.irindex(unc_lam[1]-unc_lam[0], 
				line_name)
			# Line strength of convolved spectrum (bestfit - emission lines)
			spectra = spectrum(lam, conv_spec)
			conv_index_value, conv_index_va, _, _ = spectra.irindex(lam[1] - lam[0],
				line_name)
			# LOSVD correction (From SAURON VI: Section 4.2.2)
			corr = unc_index_value/conv_index_value
			corr_var = np.sqrt((conv_index_va/conv_index_value)**2 + 
				(unc_index_va/unc_index_value)**2)*corr
		else:
			corr = 1
			corr_var = 0

		# Reduce spec to Lick resolution
		sig_pix = np.sqrt(8.**2 - (lam[1]-lam[0])**2)/(unc_lam[1] - unc_lam[0])
		spec = gaussian_filter1d(spec, sig_pix)

		spectra = spectrum(lam=lam, lamspec=spec)
		if noise is not None: 
			variance = spectrum(lam=lam, lamspec=noise**2)
		else: variance = None
		index_value, index_va, con_band, index_band = spectra.irindex(lam[1] - lam[0], 
			line_name, varSED=variance, verbose=False)

		# Apply correction
		index_va = np.sqrt(index_va)
		index_va = np.sqrt((index_va/index_value)**2 + (corr_var/corr)**2)*(
			index_value*corr)
		index_value *= corr
	
	except IndexError:
		index_value=np.nan
		index_va = np.nan

	if noise is not None:
		return index_value, index_va
	else:
		return index_value