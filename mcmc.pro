;; ==================================================================
;;               		MCMC
;; ==================================================================
;; warrenj 20150811 Rountine to run an MCMC code on the vel and sigma
;; inital guesses for the ppxf wrapper scripts. 

pro mcmc
	v=0
	sigma=200
	chi=1
repeats = 100
results = MAKE_ARRAY(2,repeats)
for i = 0, repeats-1 do begin
print, i
;print, v, sigma
	v_sav = v
	sigma_sav = sigma
	chi_sav = chi

	mcmc_fit_bin, v, sigma, chi
	results[*,i]=[v,sigma]

IF (chi GT chi_sav) AND (i LT 499) THEN BEGIN
	v=v_sav + (randomu(seed)*10)
	sigma=sigma_sav + (randomu(seed)*10)
ENDIF

endfor 

print, MEAN(results[0,*]), MEAN(results[1,*])
show_rst = SCATTERPLOT(results[0,*],results[1,*], XTITLE="velocity", YTITLE="velocity dispersion")


return
end
