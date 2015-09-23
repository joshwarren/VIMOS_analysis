## ==================================================================
## Run full IDL analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 


from kinematics import kinematics
from plot_results import plot_results
from GH_plots import GH_plots


galaxy = 'ic1459'
discard = 2
wav_range = '4200-'
vLimit = 2

plot_results(galaxy, discard=discard, wav_range=wav_range, vLimit=vLimit)

kinematics(galaxy, discard=discard, wav_range=wav_range)

GH_plots(galaxy, wav_range=wav_range)
