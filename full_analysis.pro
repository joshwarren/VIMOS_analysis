;; ==================================================================
;; Run full analysis
;; ==================================================================
;; warrenj 20150918 Routine to call a neccessary wrapper scripts for
;; binning, finding best inital guesses, finding templates to use, and
;; actually running the pPXF and Gandalf code. 
;; By default all routine will analyse NGC3557

pro full_analysis

galaxy = 'ic1459'
z = 0.01
discard = 2
SN = 30.0
range = [4200, 10000]
vel = 0.0d
sig = 200.0d


binning_spaxels, galaxy, discard, SN

find_template, galaxy, z, discard, range

mcmc, galaxy, z, discard, range, vel, sig

gandalf_VIMOS, galaxy, z, discard, range, vel, sig



return
end
