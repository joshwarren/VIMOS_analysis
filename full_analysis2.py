## ==================================================================
## Run full python analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
from pickler import pickler
from sav_for_kinemetry import sav_for_kinemetry
from plot_results import plot_results, mapping
from kinematics import kinematics
# from GH_plots import GH_plots
from rotation_curve import rotation_curve
from plot_absorption import plot_absorption
import matplotlib.pyplot as plt # used for plotting
from plot_stellar_pop import plot_stellar_pop
from use_kinemetry import use_kinemetry
from classify import classify
from fit_disk import fit_disk
import traceback, sys

galaxies = [
			# 'eso443-g024',
			'ic1459',
			# 'ic1531', 
			# 'ic4296',
			# 'ngc0612',
			# 'ngc1399',
			# 'ngc3100',
			'ngc3557',
			# 'ngc7075',
			# 'pks0718-34'
			]
# galaxies = ['eso443-g024']
# galaxies = ['ic1459']
# galaxies = ['ic1531']
# galaxies = ['ic4296']
# galaxies = ['ngc0612']
# galaxies = ['ngc1399']
# galaxies = ['ngc3100']
# galaxies = ['ngc3557']
# galaxies = ['ngc7075']
# galaxies = ['pks0718-34']


discard = 0
norm= ''#'lwv' # 'lws'
opt_dir=''

m=mapping()
# m.SNR = False
# m.image = False
# m.equivalent_width = False
# m.amp_noise = False
# m.kinematics = False
# m.plot_resid = False
# m.line_ratios = False

# Arrays for error catching
gal_err=[]
err = []
trace =[]
for galaxy in galaxies:
	D = None
	print galaxy
	try:
		# D = pickler(galaxy, discard=discard, norm=norm, opt='kin'+opt_dir)
		# D = sav_for_kinemetry(galaxy, opt='kin'+opt_dir)
		# D = plot_results(galaxy, discard=discard, overplot={'CO':'c', 'radio':'r'}, 
		# 	residual="median", norm=norm, D=D, mapping=m, opt='kin'+opt_dir,
		# 	show_bin_num=True)
		# plt.close("all")
		# # GH_plots(galaxy)
		# plt.close("all")
		# kinematics(galaxy, discard=discard, D=D, opt='kin'+opt_dir)
		# rotation_curve(galaxy, D=D, opt='kin'+opt_dir)
		# plt.close("all")
		# fit_disk(galaxy, opt='kin'+opt_dir, D=D)
		# plt.close('all')

		# Requires the IDL kinemetry routine to have been run. 
		# classify(galaxy, opt='kin'+opt_dir)
		# use_kinemetry(galaxy, opt='kin'+opt_dir)

		D = None
		D = pickler(galaxy, discard=discard, norm=norm, opt='pop'+opt_dir)
		D = plot_results(galaxy, discard=discard, overplot={'CO':'c', 'radio':'r'}, 
			residual="median", norm=norm, D=D, mapping=m, opt='pop'+opt_dir,
			show_bin_num=True)
		D = plot_absorption(galaxy, D=D, opt='pop'+opt_dir, uncert=True, 
			overplot={'CO':'c', 'radio':'r'})
		D = plot_stellar_pop(galaxy, method='mostlikely', D=D, opt='pop'+opt_dir, 
			overplot={'CO':'c', 'radio':'r'})
	except Exception as e:
		gal_err.append(galaxy)
		err.append(e)
		trace.append(sys.exc_info())
		 
#v_vd_ellip(wav_range=wav_range)

# Display errors
for i in range(len(gal_err)):
	print ''
	print gal_err[i], ' FAILED'
	print err[i]
	exc_type, exc_value, exc_traceback = trace[i]
	traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
	print ''