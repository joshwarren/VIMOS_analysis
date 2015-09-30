## ==================================================================
## Run full IDL analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 


from kinematics import kinematics
from plot_results import plot_results
from GH_plots import GH_plots
import matplotlib.pyplot as plt # used for plotting
import os # for creating directory


galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']

galaxy = galaxies[-2]
discard = 2
wav_range = '4200-'
vLimit = 2

path = "/Data/vimosindi/analysis/%s/results/%s/plots" % (galaxy, wav_range)
if not os.path.exists(path):
    os.makedirs(path)

plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit)
plt.close("all")
kinematics(galaxy, discard=discard, wav_range=wav_range)
plt.close("all")
GH_plots(galaxy, wav_range=wav_range)
