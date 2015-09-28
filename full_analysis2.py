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


galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612']

galaxy = galaxies[-1]
discard = 2
wav_range = '4200-'
vLimit = 2

os.makedirs("/Data/vimosindi/analysis/%s/results/%s/plots" % (galaxy, wav_range))

plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit)
plt.close("all")
kinematics(galaxy, discard=discard, wav_range=wav_range)
plt.close("all")
GH_plots(galaxy, wav_range=wav_range)
