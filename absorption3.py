from spectools import *
from Thomas_models_index_model_variation import *
import numpy as np 

def absorption(line_name, lam, spec, unc_lam=None, unc_spec=None, conv_spec=None, 
	noise=None):
	if line_name=='H_beta' or line_name=='Hbeta':
		line_name='Hb'
	elif line_name=='Mg_b':
		line_name='Mgb'

	spectra = spectrum(lam, spec)
	variance = spectrum(lam, noise)
	try:
		index_value, index_va, continuum_def,feature_def = spectra.irindex(0.71,
			line_name,varSED=variance)
	except ValueError:
		index_value=np.nan
		index_va = np.nan

	if noise is not None:
		return index_value, index_va
	else:
		return index_value