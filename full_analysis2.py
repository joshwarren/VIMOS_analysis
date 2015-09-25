## ==================================================================
## Run full IDL analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 


from kinematics import kinematics
from plot_results import plot_results
from GH_plots import GH_plots
import matplotlib.pyplot as plt # used for plotting



galaxy = 'ic1459'
galaxy = 'ic1531'
discard = 2
wav_range = '4200-'
vLimit = 2

plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit)
plt.close("all")
kinematics(galaxy, discard=discard, wav_range=wav_range)
plt.close("all")
GH_plots(galaxy, wav_range=wav_range)
