# ==================================================================
# Run full IDL analysis
# ==================================================================
# warrenj 20150918 Routine to call a neccessary wrapper scripts for
# binning, finding best inital guesses, finding templates to use, and
# actually running the pPXF and Gandalf code. 
# By default all routine will analyse NGC3557

import numpy as np
from binning_spaxels import binning_spaxels
from find_template import find_template
from mcmc import mcmc

def full_analysis(galaxy=None, opt='kin'):

	galaxies = np.array(['ngc3557',
		'ic1459', 
		'ic1531', 
		'ic4296', 
		'ngc0612',
		'ngc1399',
		'ngc3100',
		'ngc7075', 
		'pks0718-34', 
		'eso443-g024'])

	# an inital guess from quick internet search of redshift.
	z_gals = [0.01, 0.005, 0.025, 0.01, 0.028, 0.005, 0.01, 0.02, 0.03, 0.015] 
	gal=6
	if galaxy is None:
		galaxy = galaxies[gal]
	else:
		gal = np.where(galaxies == galaxy)[0][0]
	
	print galaxy
	
	z = z_gals[gal]
	discard = 0
	targetSN = 30.
	set_range = [4200, 10000]

	binning_spaxels(galaxy, discard=discard, targetSN=targetSN, opt=opt, 
		auto_override=True)

	# find_template(galaxy, z=z, discard=discard, set_range=set_range)

	# mcmc(galaxy, z=z, discard=discard, set_range=set_range)


if __name__=="__main__":
	galaxies = [#'ngc3557',
		# 'ic1459'#,
		#'ic1531', 
		# 'ic4296', 
		# 'ngc0612',
		# 'ngc1399',
		'ngc3100',
		# 'ngc7075', 
		# 'pks0718-34', 
		# 'eso443-g024'
		]
	for g in galaxies: full_analysis(galaxy=g, opt='kin')
	# full_analysis(galaxy='ic1531', opt='pop')