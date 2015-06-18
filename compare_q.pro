;; ==================================================================
;; 		Plot 4 central spaxels on top of each other
;; ==================================================================
;; warrenj 20150618 Compare the quadrants by plotting the 4 central
;; spaxels on top of each other. 



pro compare_q

galaxy = 'ngc3557'

;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data, header



	spectrum_lin = MAKE_ARRAY(4,n_elements(galaxy_data[19,19,*]))

;; ----------============ Spaxel [19,19]  ================---------
for i = 0, n_elements(spectrum_lin[0,*])-1 do begin
	spectrum_lin[0,i] = galaxy_data[19, 19, i]
     endfor
spectrum_lin[0,*] = spectrum_lin[0,*]/MEDIAN(spectrum_lin[0,*])

;; ----------============ Spaxel [20,19]  ================---------
for i = 0, n_elements(spectrum_lin[0,*])-1 do begin
	spectrum_lin[1,i] = galaxy_data[20, 19, i]
endfor
spectrum_lin[1,*] = spectrum_lin[1,*]/MEDIAN(spectrum_lin[1,*])

;; ----------============ Spaxel [20,20]  ================---------
for i = 0, n_elements(spectrum_lin[0,*])-1 do begin
	spectrum_lin[2,i] = galaxy_data[20, 20, i]
endfor
spectrum_lin[2,*] = spectrum_lin[2,*]/MEDIAN(spectrum_lin[2,*])

;; ----------============ Spaxel [19,20]  ================---------
for i = 0, n_elements(spectrum_lin[0,*])-1 do begin
	spectrum_lin[3,i] = galaxy_data[19, 20, i]
endfor
spectrum_lin[3,*] = spectrum_lin[3,*]/MEDIAN(spectrum_lin[3,*])




plot, spectrum_lin[0,*]
oplot, spectrum_lin[1,*], color = 500
oplot, spectrum_lin[2,*], color = 1000000
oplot, spectrum_lin[3,*], color = 60000000


return
end

