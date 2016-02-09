;; ==================================================================
;;               		MCMC
;; ==================================================================
;; warrenj 20150811 Rountine to run an MCMC code on the vel and sigma
;; inital guesses for the ppxf wrapper scripts. 

pro mcmc, galaxy, z=z, discard=discard, range=range, vel=vel, sig=sig

;; ----------===============================================---------
;; ----------============ Default parameters ===============---------
;; ----------===============================================---------

if not keyword_set(z) then z=0.01
if not keyword_set(discard) then discard=2
if not keyword_set(range) then range=[4200,10000]
if not keyword_set(vel) then vel=0.0d
if not keyword_set(sig) then sig=200.0d



;	galaxy='ic1459'
;	z = 0.01
;	vel=0
;	sig=200


	chi=1
	repeats = 100
	c = 299792.458d


results = MAKE_ARRAY(2,repeats-2)
for i = 0, repeats-1 do begin
print, i
	v_sav = vel
	sigma_sav = sig
	chi_sav = chi

	mcmc_fit_bin, galaxy, z, vel, sig, chi

if i le 1 then z = (z + 1)*((1 + vel/c)/(1 - vel/c)) - 1 else begin

	results[*,i-2]=[vel,sig]
	local_min=randomu(seed0)
IF ((chi GT chi_sav) AND (i LT 499)) OR local_min GT 0.85 THEN BEGIN
	vel=v_sav + (randomu(seed1)*10)
	sig=sigma_sav + (randomu(seed2)*10)
ENDIF
endelse
endfor 

print, MEAN(results[0,*]), MEAN(results[1,*])
show_rst = SCATTERPLOT(results[0,*],results[1,*], XTITLE="velocity", YTITLE="velocity dispersion")


vel = MEAN(results[0,*])
sig = MEAN(results[1,*])



;; ----------===============================================---------
;; ----------================= Save Result =================---------
;; ----------===============================================---------
data_file = "/Data/vimosindi/analysis/galaxies.txt"
readcol, data_file, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, $
    y_gals, SN_used, skipline=1, format='A,D,D,D,I,I,I', /SILENT

i_gal = where(galaxy_gals eq galaxy)
if i_gal eq -1 then begin
    galaxy_gals = [galaxy_gals, galaxy]
    z_gals = [z_gals, z]
    vel_gals = [vel_gals, vel]
    sig_gals = [sig_gals, sig]
endif else begin
    z_gals[i_gal] = z
    vel_gals[i_gal] = vel
    sig_gals[i_gal] = sig
endelse


forprint, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_used, $
    textout=data_file, /SILENT, $
    Comment = "Galaxy      z     velocity     velocity dispersion    x    y    Target SN"



return
end
