;; ==================================================================
;; 		Plot given spaxels on top of each other
;; ==================================================================
;; warrenj 20150618 Compare spaxels



pro compare_spaxels

galaxy = 'ngc3557'

;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data, header


spaxels = [[10,37],[10,38],[10,39],[11,36],[11,37],[11,38],[11,39]]

s = size(spaxels)

	spectrum_lin = MAKE_ARRAY(s[2],n_elements(galaxy_data[0,0,*]))


for i = 0, s[2]-1 do begin
for k = 0, n_elements(spectrum_lin[0,*])-1 do begin
	spectrum_lin[i,k] = galaxy_data[spaxels[0,i], spaxels[1,i], k]
endfor
if (MEDIAN(spectrum_lin[i,*]) NE 0) then begin
spectrum_lin[i,*] = spectrum_lin[i,*]/MEDIAN(spectrum_lin[i,*])
endif
endfor



;; ----------=========== Setting the x axis  ==============---------
x = MAKE_ARRAY(n_elements(spectrum_lin[0,*]))
for i = 0, n_elements(x)-1 do begin
x[i] =  sxpar(header,'CRVAL3') + sxpar(header,'CD3_3')*i
endfor



plot, x, spectrum_lin[0,*]
for i=1,s[2]-1 do begin
oplot, x, spectrum_lin[i,*], color = 10^(spaxels[0,i]+spaxels[1,i]-35) 
endfor


return
end

