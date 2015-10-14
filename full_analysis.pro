;; ==================================================================
;; Run full IDL analysis
;; ==================================================================
;; warrenj 20150918 Routine to call a neccessary wrapper scripts for
;; binning, finding best inital guesses, finding templates to use, and
;; actually running the pPXF and Gandalf code. 
;; By default all routine will analyse NGC3557

pro full_analysis

galaxies = ['ngc3557', $
	'ic1459', $
	'ic1531', $
	'ic4296', $
	'ngc0612', $
	'ngc1399', $
	'ngc3100', $
	'ngc7075', $
	'pks0718-34', $
	'eso443-g024']
; an inital guess from quick internet search of redshift.
z_gals = [0.01, 0.005, 0.02, 0.01, 0.03, 0.005, 0.01, 0.02, 0.03, 0.015] 
;for gal=0, n_elements(z_gals)-1 do begin
gal=1
galaxy = galaxies[gal]
print, galaxy
z = z_gals[gal]
discard = 2
targetSN = 30.0
range = [4200, 10000]


binning_spaxels, galaxy, discard=discard, targetSN=targetSN

;find_template, galaxy, z=z, discard=discard, range=range

;mcmc, galaxy, z=z, discard=discard, range=range

;gandalf_VIMOS, galaxy, discard=discard, range=range
;endfor
return
end
