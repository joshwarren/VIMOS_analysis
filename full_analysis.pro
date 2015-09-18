;; ==================================================================
;; Run full analysis
;; ==================================================================
;; warrenj 20150918 Routine to call a neccessary wrapper scripts for
;; binning, finding best inital guesses, finding templates to use, and
;; actually running the pPXF and Gandalf code. 
;; By default all routine will analyse NGC3557

pro full_analysis

galaxy = 'ic1459'




find_template, galaxy

binning_spaxels, galaxy

mcmc, galaxy

gandalf_VIMOS, galaxy



return
end
