;; ==================================================================
;; Run full IDL analysis
;; ==================================================================
;; warrenj 20150918 Routine to call a neccessary wrapper scripts for
;; binning, finding best inital guesses, finding templates to use, and
;; actually running the pPXF and Gandalf code. 
;; By default all routine will analyse NGC3557

pro full_analysis

galaxy = 'ic1459'
galaxy = 'ngc3557'
z = 0.01
discard = 2
SN = 30.0
range = [4200, 10000]


;binning_spaxels, galaxy, discard, SN

find_template, galaxy, z=z, discard=discard, range=range

mcmc, galaxy, z=z, discard=discard, range=range

;gandalf_VIMOS, galaxy, discard=discard, range=range



return
end
