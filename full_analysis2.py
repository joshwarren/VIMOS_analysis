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

galaxies = ['ic1459',
            'ic1531',
            'ic4296',
            'ngc0612',
            'ngc3100',
            'ngc7075',
            'pks0718-34',
            'ngc1399',
            'ngc3557',
            'eso443-g024']
#galaxies = ['ngc3557']
#galaxies = ['ic1459']
#galaxies = ['ic1531']
#galaxies = ['ic4296']
#galaxies = ['ngc0612']
#galaxies = ['ngc1399']
#galaxies = ['ngc3100']
#galaxies = ['ngc7075']
#galaxies = ['pks0718-34']
#galaxies = ['eso443-g024']


discard = 2
wav_range = '4200-'
vLimit = 2
norm='lwv'
for galaxy in galaxies:

    print galaxy
    man_errors(galaxy, wav_range=wav_range)
    pickler(galaxy, discard=discard, wav_range=wav_range, norm=norm)
    #plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit, 
    #    nointerp=True, CO = True, residual="median", norm=norm)
    #plt.close("all")
    #GH_plots(galaxy, wav_range=wav_range)
    #plt.close("all")
    #kinematics(galaxy, discard=discard, wav_range=wav_range)
    #plt.close("all")
    #plot_absorption(galaxy, wav_range=wav_range, vLimit=vLimit)
    stellar_pop(galaxy, wav_range=wav_range, vLimit=vLimit)
#v_vd_ellip(wav_range=wav_range)
