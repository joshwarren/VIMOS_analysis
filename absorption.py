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
##############################################################################
from spectools import *
from tools import length as len
#from Thomas_models_index_model_variation import *
from scipy.ndimage.filters import gaussian_filter1d
import numpy as np 

def absorption(line_name, lam, spec, unc_lam=None, unc_spec=None, conv_spec=None, 
	noise=None):
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
			# Line strength of unconvolved (/convolved to LICK resolution) spectrum.
			sig_pix = 200*np.median(lam)/c/(unc_lam[1]-unc_lam[0])
			lick_spec = gaussian_filter1d(unc_spec, sig_pix)
			spectra = spectrum(unc_lam, unc_spec)
			unc_index_value, unc_index_va, _, _ = spectra.irindex(0.71, line_name)
			# Line strength of convolved spectrum (bestfit - emission lines)
			spectra = spectrum(lam, conv_spec)
			conv_index_value, conv_index_va, _, _ = spectra.irindex(0.71, line_name)
			# LOSVD correction (From SAURON VI: Section 4.2.2)
			corr = unc_index_value/conv_index_value
			corr_var = np.sqrt((conv_index_va/conv_index_value)**2 + 
				(unc_index_va/unc_index_value)**2)*corr
		else:
			corr = 1
			corr_var = 0

		spectra = spectrum(lam, spec)
		if noise is not None: variance = spectrum(lam, noise)
		else: variance = None
		index_value, index_va, _, _ = spectra.irindex(0.71, line_name, varSED=variance)

		# Apply correction
		index_va = np.sqrt((index_va/index_value)**2 + (corr_var/corr)**2)*(index_value*corr)
		index_value *= corr
	
	except IndexError:
		index_value=np.nan
		index_va = np.nan

	if noise is not None:
		return index_value, index_va
	else:
		return index_value