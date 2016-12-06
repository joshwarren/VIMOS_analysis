from spectools import *
from Thomas_models_index_model_variation import *
import numpy as np 

def absorption(line_name, D, uncert=False):
	if line_name=='H_beta' or line_name=='Hbeta':
		line_name='Hb'
	elif line_name=='Mg_b':
		line_name='Mgb'

	spectra=[]
	variance=[]
	index1_ob=[]
	index1_ob_var=[]
	for i in range(D.number_of_bins):
		spectra.append(spectrum(D.bin[i].lam, D.bin[i].continuum))
		variance.append(spectrum(D.bin[i].lam, D.bin[i].noise))
		try:
			index_value, index_va, continuum_def,feature_def = spectra[i].irindex(0.71,line_name,varSED=variance[i])
		except ValueError:
			index_value=np.nan
			index_va = np.nan
		index1_ob.append(index_value)
		index1_ob_var.append(index_va)

	if uncert:
		return np.array(index1_ob).flatten(), np.array(index1_ob_var).flatten()
	else:
		return np.array(index1_ob).flatten()