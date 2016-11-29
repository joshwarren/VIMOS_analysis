## ==================================================================
## Run full python analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 

from pickler import pickler
from plot_results import plot_results
from kinematics import kinematics
from GH_plots import GH_plots
from man_errors2 import man_errors
from plot_absorption import plot_absorption
import matplotlib.pyplot as plt # used for plotting
from stellar_pop import stellar_pop
from use_kinemetry import use_kinemetry
from classify import classify

galaxies = ['ic1459',
			'ic1531',
			'ic4296',
			'ngc0612',# missing 216
			'ngc1399',
			'ngc3100',
			'ngc7075',
			'pks0718-34',
			'ngc3557', # missing 686, 688
			'eso443-g024'] # missing 453, 455
#galaxies = ['ngc3557']
#galaxies = ['ic1459']
#galaxies = ['ic1531']
galaxies = ['ic4296']
#galaxies = ['ngc0612']
#galaxies = ['ngc1399']
galaxies = ['ngc3100']
#galaxies = ['ngc7075']
#galaxies = ['pks0718-34']
#galaxies = ['eso443-g024']


discard = 2
wav_range = '4200-'
vLimit = 2
norm='lwv'

# Arrays for error catching
gal_err=[]
err = []
for galaxy in galaxies:
	D = None
	print galaxy
	#try:
		#man_errors(galaxy, wav_range=wav_range)
		#D = pickler(galaxy, discard=discard, wav_range=wav_range, norm=norm)
		#D = plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit, 
		#	nointerp=True, CO = True, residual="median", norm=norm, D=D)
	# 	plt.close("all")
	# 	#GH_plots(galaxy, wav_range=wav_range)
	# 	#plt.close("all")
	#	kinematics(galaxy, discard=discard, wav_range=wav_range)
	# 	#plt.close("all")
	#	use_kinemetry(galaxy)
	#	classify(galaxy)
		
	# D = pickler(galaxy, discard=discard, wav_range=wav_range, norm=norm, opt='pop')
	# D = plot_absorption(galaxy, wav_range=wav_range, vLimit=vLimit, D=D)#, uncert=False)
	D = stellar_pop(galaxy, wav_range=wav_range, vLimit=vLimit, D=D)
	#except Exception as e:
	#	gal_err.append(galaxy)
	#	err.append(e)
#v_vd_ellip(wav_range=wav_range)

# Display errors
for i in range(len(gal_err)):
	print ''
	print gal_err[i], ' FAILED'
	print err[i]
	print ''