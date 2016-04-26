## ==================================================================
## Run full python analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 


from kinematics import kinematics
#from plot_results_CO import plot_results
from plot_results import plot_results
from GH_plots import GH_plots
from man_errors2 import man_errors
import matplotlib.pyplot as plt # used for plotting
import os # for creating directory

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024', 'ngc1399']
#galaxies = ['ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
galaxies = ['ngc1399']
#galaxies = ['ngc7075']
#galaxies = ['eso443-g024']
#galaxies = ['pks0718-34']
#galaxies = ['ic1459']
galaxies = ['ngc3100']
#galaxies = ['ngc0612']

discard = 2
wav_range = '4200-'
vLimit = 1
for galaxy in galaxies:

    path = "/Data/vimosindi/analysis/%s/results/%s/plots" % (galaxy, wav_range)
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(path+"/notinterpolated"):
        os.makedirs(path+"/notinterpolated")

    print galaxy
#    man_errors(galaxy, wav_range=wav_range)
#    plt.close("all")
    plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit, 
        nointerp=True, residual=False)#residual="median")#, norm='sig')
    plt.close("all")
#    plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit, 
#        nointerp=False)
#    plt.close("all")
#    GH_plots(galaxy, wav_range=wav_range)
#    plt.close("all")
#    kinematics(galaxy, discard=discard, wav_range=wav_range)
#    plt.close("all")
#v_vd_ellip(wav_range=wav_range)
