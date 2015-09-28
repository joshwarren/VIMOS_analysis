;; ==================================================================
;; Run full IDL analysis
;; ==================================================================
;; warrenj 20150918 Routine to call a neccessary wrapper scripts for
;; binning, finding best inital guesses, finding templates to use, and
;; actually running the pPXF and Gandalf code. 
;; By default all routine will analyse NGC3557

pro full_analysis

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612']
; an inital guess from quick internet search of redshift.
z_gals = [0.01, 0.005, 0.025, 0.01, 0.03] 
gal = -1
galaxy = galaxies[gal]
z = z_gals[gal]
discard = 2
targetSN = 30.0
range = [4200, 10000]


;binning_spaxels, galaxy, discard=discard, targetSN=targetSN

;find_template, galaxy, z=z, discard=discard, range=range

;mcmc, galaxy, z=z, discard=discard, range=range

gandalf_VIMOS, galaxy, discard=discard, range=range



return
end
