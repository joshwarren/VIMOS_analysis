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
# unc_lam:		None	Wavelength of unconvolved spectrum - from the templates used by 
#						pPXF
# unc_spec:		None	Unconcolved spectum
# conv_spec:	None	Convolved spectrum i.e. pPXF bestfit
# noise:		None	Observed noise from reduction pipeline
# lick:			False	Return corrected indices to LICK resolution
##############################################################################
from spectools import *
# from tools import length as len
import numpy as np 

c = 299792.458 # speed of light in km/s


def absorption(line_name, lam, spec, unc_lam=None, unc_spec=None, 
	conv_spec=None, conv_lam=None, noise=None):
	mag = False
	if line_name=='H_beta' or line_name=='Hbeta':
		line_name='Hb'
		line_name2='hb'
	elif line_name=='Mg_b' or line_name=='Mgb':
		line_name='Mgb'
		line_name2='mgb'
	elif line_name=='Mg1':
		line_name2='mg1'
		mag = True
	elif line_name=='Mg2':
		line_name2='mg2'
		mag = True
	elif line_name=='NaD':
		line_name2='nad'
	elif line_name=='TiO1':
		line_name2='tio1'
		mag = True
	elif line_name=='TiO2':
		line_name2='tio2'
		mag = True
	elif line_name=='Fe5270':
		line_name2='fe52'
	elif line_name=='Fe5335':
		line_name2='fe53'
	elif line_name=='HdeltaA' or line_name=='HdA':
		line_name2='hdA'
	else:
		line_name2=line_name

	if unc_lam is not None and unc_spec is not None:
		if len(unc_lam) != len(unc_spec):
			raise ValueError('Length of unconvoled spectrum and corresponding '+
				'wavelength should be the same.')
	if conv_spec is not None and conv_lam is None:
		if len(lam) != len(conv_spec) != len(spec):
			raise ValueError('Length of spectrum (bestfit and observed) and '+
				'corresponding wavelength should be the same.')
		else:
			conv_lam = np.array(lam)
	elif conv_spec is not None and conv_lam is not None:
		if len(conv_lam) != len(conv_spec):
			raise ValueError('Convolved wavelength and spectrum must be the same'+
				' length')
	if noise is not None:
		if len(noise) != len(spec):
			raise ValueError('Length of noise and observed spectrum '+
				'should be the same.')

	# spec /= spec[np.argmin(np.abs(lam-4300))]

	# Catch if index is not in wavelenght range.
	try:
		# if disp is not None:
			# sig_pix = np.sqrt(8.4**2 - disp**2)/(disp)
			# spec = gaussian_filter1d(spec, sig_pix)
		spectra = spectrum(lam=lam, lamspec=spec)
		disp = np.diff(lam)[np.argmin(np.abs(lam 
			- np.mean(getattr(spectra,line_name2))))]
		if noise is not None: 
			variance = spectrum(lam=lam, lamspec=noise**2)
		else: 
			variance = None

		if mag:
			index_value, index_va = spectra.irmagindex(disp, line_name, 
				varSED=variance, verbose=False)
		else:
			index_value, index_va, cont, band = spectra.irindex(disp, 
				line_name, varSED=variance, verbose=False)

		# Calculate correction for velocity dispersion spreading
		if unc_lam is not None and unc_spec is not None and \
			conv_spec is not None:

			# Line strength of unconvolved spectrum from templates
			spectra = spectrum(unc_lam, unc_spec)
			unc_disp = np.diff(unc_lam)[np.argmin(np.abs(unc_lam - 
				np.mean(getattr(spectra,line_name2))))]
			if mag:
				unc_index_value, unc_index_va = spectra.irmagindex(
					unc_disp, line_name)
			else:
				unc_index_value, unc_index_va, _, _ = spectra.irindex(
					unc_disp, line_name)

			# Line strength of convolved spectrum
			spectra = spectrum(conv_lam, conv_spec)
			conv_disp = np.diff(conv_lam)[np.argmin(np.abs(conv_lam - 
				np.mean(getattr(spectra,line_name2))))]
			if mag:
				conv_index_value, conv_index_va = spectra.irmagindex(
					conv_disp, line_name)
			else:
				conv_index_value, conv_index_va, _, _ = spectra.irindex(
					conv_disp, line_name)
			# LOSVD correction (From SAURON VI: Section 4.2.2)
			corr = unc_index_value/conv_index_value
			# corr_var = np.sqrt((conv_index_va/conv_index_value)**2 + 
			# 	(unc_index_va/unc_index_value)**2)*corr
		else:
			corr = 1
			# corr_var = 0


		# Apply correction
		index_va = np.sqrt(index_va) * corr
		# index_va = np.sqrt((index_va/index_value)**2 + (corr_var/corr)**2)*(
		# 	index_value*corr)
		index_value *= corr
	
	except IndexError:
		index_value = np.array([np.nan])
		index_va = np.array([np.nan])

	if noise is not None:
		return index_value, index_va
	else:
		return index_value

# Variable Lick resolution
def get_Lick_res(index):
	# Find band ; dummy spectra
	s = spectrum(np.array([100,101,102]),np.array([1,1,1]))	
	if index=='Mgb' or index=='Mg_b':
		cont = s.mgbcont; band = s.mgb
	elif index=='Hb' or index=='Hbeta' or index=='H_beta':
		cont = s.hbcont; band = s.hb
	elif index=='G4300':
		cont = s.G4300cont; band = s.G4300
	elif index=='NaD':
		cont = s.nadcont; band = s.nad
	elif index=='TiO1':
		cont = s.tio1cont; band = s.tio1
	elif index=='TiO2':
		cont = s.tio2cont; band = s.tio2
	elif index=='Fe4383':
		cont = s.Fe4383cont; band = s.Fe4383
	elif index=='Ca4455':
		cont = s.Ca4455cont; band = s.Ca4455
	elif index=='Fe4531':
		cont = s.Fe4531cont; band = s.Fe4531
	elif index=='Fe4668':
		cont = s.Fe4668cont; band = s.Fe4668
	elif index=='Fe5015':
		cont = s.Fe5015cont; band = s.Fe5015
	elif index=='Fe52' or index=='Fe5270':
		cont = s.fe52cont; band = s.fe52
	elif index=='Fe53' or index=='Fe5335':
		cont = s.fe53cont; band = s.fe53
	elif index=='Fe5406':
		cont = s.Fe5406cont; band = s.Fe5406
	elif index=='Fe5709':
		cont = s.Fe5709cont; band = s.Fe5709
	elif index=='Fe5782':
		cont = s.Fe5782cont; band = s.Fe5782


	# Wavelength of center of band
	wav = np.median(band)

	from scipy.interpolate import interp1d
	# From table 8 Worthey, Ottaviani 1997 ApJS 111 (2) 377
	from scipy import __version__ as scipy_v
	if scipy_v == '0.19.0':
		res_wav_dep = interp1d(np.array([4000, 4400, 4900, 5400, 6000]), 
			np.array([11.5, 9.2, 8.4, 8.4, 9.8]), bounds_error=False, 
			fill_value=(11.5, 9.8))
		out = res_wav_dep(wav)
	else: 
		if wav < 4000:
			out = 11.5
		elif wav > 6000:
			out = 9.8
		else:
			res_wav_dep = interp1d(np.array([4000, 4400, 4900, 5400, 6000]), 
				np.array([11.5, 9.2, 8.4, 8.4, 9.8]))
			out = res_wav_dep(wav)
	return out