;; ==================================================================
;; Analyse reduced VIMOS data using pPFX
;; ==================================================================
;; warrenj 20150216 Process to analyse the reduced VIMOS data.

pro run_analysis;, galaxy

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'ngc3557'
	c = 299792.458d
  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
	vel = 2000.0d ; Initial estimate of the galaxy velocity in km/s
	sig = 270.0d ; Initial estimate of the galaxy dispersion in km/s 
		     ; within its rest frame
        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
        FWHM_gal = FWHM_gal/(1+z) ; Adjust resolution in Angstrom
	moments = 4 ; number of comonants to calc with ppxf (see 
                    ; keyword moments in ppxf.pro for more details)
;; File for output: an array containing the calculated dynamics of the
;; galaxy. 
	output_v = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_vel.dat'
	output_temp_weighting = '/Data/vimosindi/analysis/' + $
		galaxy + '/results/template_weighting.dat'
	output_sigma = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_sigma.dat'
	output_h3 = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_h3.dat'
	output_h4 = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_h4.dat'
	output_h5 = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_h5.dat'
	output_h6 = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_h6.dat'
	output_Chi = '/Data/vimosindi/analysis/' + galaxy + $
		'/results/gal_Chi.dat'
	
;; Tessellation input
;	binning_spaxels, galaxy
	tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'


;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------
	CLOSE, 1, 2, 3, 4, 5, 6, 7, 8
	OPENW, 1, output_temp_weighting
	OPENW, 2, output_v
	OPENW, 3, output_sigma
	OPENW, 4, output_h3
	OPENW, 5, output_h4
	OPENW, 6, output_h5
	OPENW, 7, output_h6
	OPENW, 8, output_Chi




;+
;;; ----------=============== Elodie library ================---------
; Finding the template files
;	templatesDirectory = '/Data/ppxf/LL_ELODIE_31_23/'
;	templateFiles = FILE_SEARCH(templatesDirectory + '*.fits', $
;		COUNT=nfiles)
;
;
;;; Creating the templates array with correct dimension
;	FITS_READ, templateFiles[0], temp_template, h2
;	lamRange2 = sxpar(h2,'CRVAL1') + $
;		[0d, sxpar(h2,'CDELT1') * (sxpar(h2, 'NAXIS1') - 1d)]
;	log_rebin, lamRange2, temp_template, log_temp_template, $
;                logLam_template, velscale=velscale
;	templates = MAKE_ARRAY(n_elements(log_temp_template), nfiles)
;
;
;       
;
;        
;	FWHM_tem = 0.5     ; Elodie spectra have a resolution
;	                   ; FWHM of 0.5A.
;         
;;; Reading the contents of the fits files into the array
;;; templates. Including rebinning them.
;for i = 0, nfiles - 1 do begin
;
;	FITS_READ, templateFiles[i], temp_template, h2
;
;;; Calibrating resolutions
;	FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
;	sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference
;                                		  ; in pixels
;	temp_template = gauss_smooth(temp_template, sigma)
;
;;; Rebinning templates logarthmically
;	lamRange_template = sxpar(h2,'CRVAL1') + $
;		[0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
;	log_rebin, lamRange_template, temp_template, log_temp_template, $
;		velscale=velscale
;
;;; Normalizing templates
;        templates[*,i] = log_temp_template/median(log_temp_template)
;endfor
;-




;+
;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/ppxf/MILES_library/'
	templateFiles = FILE_SEARCH(templatesDirectory + $
		'm0[0-9][0-9][0-9]V', COUNT=nfiles)

;v1 is wavelength, v2 is spectrum
	READCOL, templateFiles[0], v1,v2, FORMAT = 'D,D', /SILENT

; Using same keywords as fits headers
	CRVAL1 = v1[0]		; starting wavelength
	NAXIS1 = size(v2, /N_ELEMENTS) ; Number of entries
	; wavelength increments (resolution?)
	CDELT1 = (v1[NAXIS1-1]-v1[0])/(NAXIS1-1)

; Creating the templates array with correct dimension
;	temp_template = MAKE_ARRAY(NAXIS1, nfiles)

	lamRange_template = CRVAL1 + [0d, CDELT1 * (NAXIS1 - 1d)]
        log_rebin, lamRange_template, v1, log_temp_template, $
		logLam_template, velscale=velscale

;; ****************************************************************
;; NB: shouldn't this be 0.9A as this is resolution?
        FWHM_tem = 2.5     ; Miles spectra have a resolution
                           ; FWHM of 2.5A.


;; Which templates to use are given in use_templates.pro. This is
;; transfered to the array templatesToUse.
	use_templates, templatesToUse
	nfiles = N_ELEMENTS(templatesToUse)
	templates = MAKE_ARRAY(n_elements(log_temp_template), nfiles)

         
;; Reading the contents of the files into the array templates. 
;; Including rebinning them.
for i = 0, nfiles - 1 do begin


	READCOL, templateFiles[templatesToUse[i]-1], v1,v2, $
		FORMAT = 'D,D', /SILENT


;	READCOL, templateFiles[i], v1,v2, FORMAT = 'D,D', /SILENT

;; Rebinning templates logarthmically
	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

;; Normalizing templates
	templates[*,i] = log_temp_template/median(log_temp_template)
endfor
;-




;; ----------========= Reading Tessellation  =============---------

;; Reads the txt file containing the output of the binning_spaxels
;; routine. 
	RDFLOAT, tessellation_File, x, y, bin_num, COLUMNS = [1,2,3], $
		SKIPLINE = 1, /SILENT 
	
	n_bins = max(bin_num) + 1
;; Contains the order of the bin numbers in terms of index number.
	order = sort(bin_num)




;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data, header

	n_spaxels = n_elements(galaxy_data[*,0,0]) * $
		n_elements(galaxy_data[0,*,0])


;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/sxpar(header,'CD3_3') ; Sigma difference 
						     ; in pixels

;; For rebinning logarthmically
	lamRange = sxpar(header,'CRVAL3') + $
                   [0,sxpar(header,'CD3_3')*(sxpar(header,'NAXIS3')-1)]
        lamRange = lamRange/(1+z) ; Compute approximate restframe
        			    ; wavelength range



;; ----------========== Spatially Binning =============---------

i=0
for bin = 0, n_bins-1 do begin
;for bin = 0, 5 do begin

;; Need to create a new spectrum for a new bin.
bin_lin = MAKE_ARRAY(n_elements(galaxy_data[0,0,*]))

while (i LT n_spaxels && bin EQ bin_num[order[i]]) do begin

for k = 0 , n_elements(galaxy_data[x[order[i]], y[order[i]],*]) - 1 $
	do begin 
;; add spectrums together with the bin
	bin_lin[k] = bin_lin[k] + $
		galaxy_data[x[order[i]], y[order[i]],k] 
endfor

i = i + 1
endwhile
;; bin_lin now contains linearly binned spectrum of the spatial bin. 

;; smooth spectrum to fit with templates resolution
	bin_lin = gauss_smooth(bin_lin, sigma)



;; rebin spectrum logarthmically
	log_rebin, lamrange, bin_lin, bin_log, logLam_bin, $
		velscale=velscale



;; ----------========= Assigning noise variable =============---------
;;   NOISE: vector containing the 1*sigma error (per pixel) in the
;;   galaxy spectrum, or covariance matrix describing the correlated
;;   errors in the galaxy spectrum. Of course this vector/matrix must
;;   have the same units as the galaxy spectrum.
;;   - If GALAXY is a Nx2 array, NOISE has to be an array with the
;;     same dimensions.
;;   - When NOISE has dimensions NxN it is assumed to contain the
;;     covariance matrix with elements sigma(i,j). When the errors in
;;     the spectrum are uncorrelated it is mathematically equivalent
;;     to input in PPXF an error vector NOISE=errvec or a NxN diagonal
;;     matrix NOISE=DIAG_MATRIX(errvec^2) (note squared!).
;;   - IMPORTANT: the penalty term of the pPXF method is based on the
;;     *relative* change of the fit residuals. For this reason the
;;     penalty will work as expected even if no reliable estimate of
;;     the NOISE is available (see Cappellari & Emsellem [2004] for
;;     details).
;;     If no reliable noise is available this keyword can just be set
;;     to:
	noise = MAKE_ARRAY(n_elements(bin_log), $
		VALUE = 1d)
		;galaxy*0+1 ; Same weight for all pixels




; The galaxy and the template spectra do not have the same starting
; wavelength. For this reason an extra velocity shift DV has to be
; applied to the template to fit the galaxy spectrum. We remove this
; artificial shift by using the keyword VSYST in the call to PPXF
; below, so that all velocities are measured with respect to DV. This
; assume the redshift is negligible. In the case of a high-redshift
; galaxy one should de-redshift its wavelength to the rest frame
; before using the line below (see above). 

dv = (logLam_template[0]-logLam_bin[0])*c ; km/s


; Find the pixels to ignore to avoid being distracted by gas emission
; lines or atmospheric absorbsion line.  
goodPixels = ppxf_determine_goodpixels(logLam_bin,lamRange_template,vel) 





        


	start = [vel, sig] ; starting guess

print, bin
	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics, GOODPIXELS=goodPixels, $
		MOMENTS = moments, DEGREE = 2, VSYST = dv, $
		WEIGHTS = weights, /QUIET
;;		ERROR = error


;	print, 'Best-fitting redshift z:', (z + 1)*((1 + $
;		bin_dynamics[0]/c)/(1 - bin_dynamcics[0]/c)) - 1

; CALLING SEQUENCE:
;  PPXF, templates, galaxy, noise, velScale, start, sol, BESTFIT=bestFit, $
;	BIAS=bias, CHI2DOF=chi2dof, /CLEAN, COMPONENT=component, $
;	DEGREE=degree, ERROR=error, GOODPIXELS=goodPixels, LAMBDA=lambda, $
;	MDEGREE=mdegree, MOMENTS=moments, MPOLYWEIGHTS=mpolyweights, $
;	/OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, $
;	REDDENING=reddening, REGUL=regul, REG_DIM=reg_dim, SKY=sky, $
;	VSYST=vsyst, WEIGHTS=weights




	PRINTF, 2, bin_dynamics[0]
	PRINTF, 3, bin_dynamics[1]
;	PRINTF, 4, bin_dynamics[3]
;	PRINTF, 5, bin_dynamics[4]
;	PRINTF, 6, bin_dynamics[5]
;	PRINTF, 7, bin_dynamics[6]
;	PRINTF, 8, bin_dynamics[7]



;+
;;; Write weightings for each template used to file 1.
;for k = 0, nfiles-1 do begin
;if weights[k] ne 0 then begin
;;; Use (uncomment) this line if using limited template files.
;	PRINTF, 1, string(templatesToUse[k]) + ' ' + string(weights[k])
;;; Use (uncomment) this line if using full library.
;;	PRINTF, 1, string(k+1) + ' ' + string(weights[k])
;endif
;endfor
;-


endfor 



CLOSE, 1, 2, 3, 4, 5, 6, 7, 8



;; Error check - making sure all spaxels have been read into some
;; bin. 
if (i EQ n_spaxels-1) THEN BEGIN
	print, 'ERROR: not all spaxels have been read'
endif 



return
end

























;; ==================================================================
;; 		Check the fit of a central spaxel
;; ==================================================================
;; warrenj 20150602 Routine to plot the best fit against the spectrum
;; and show the residuals, for spaxels 20,20 (zero weighted).

pro check_fit

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'ngc3557'
	c = 299792.458d
;  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
	vel = 2000.0d ; Initial estimate of the galaxy velocity in km/s
	sig = 270.0d ; Initial estimate of the galaxy dispersion in km/s 
		     ; within its rest frame
;	spaxel = [20,20] ; spaxel to read off
        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
 ;       FWHM_gal = FWHM_gal/(1+z) ; Adjust resolution in Angstrom
	moments = 4 ; number of comonants to calc with ppxf (see 
                    ; keyword moments in ppxf.pro for more details)
;; File for output: an array containing the calculated dynamics of the
;; galaxy. 




;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------



;+
;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/ppxf/MILES_library/'
	templateFiles = FILE_SEARCH(templatesDirectory + $
		'm0[0-9][0-9][0-9]V', COUNT=nfiles)

;v1 is wavelength, v2 is spectrum
	READCOL, templateFiles[0], v1,v2, FORMAT = 'D,D', /SILENT

; Using same keywords as fits headers
	CRVAL1 = v1[0]		; starting wavelength
	NAXIS1 = size(v2, /N_ELEMENTS) ; Number of entries
	; wavelength increments (resolution?)
	CDELT1 = (v1[NAXIS1-1]-v1[0])/(NAXIS1-1)

; Creating the templates array with correct dimension
;	temp_template = MAKE_ARRAY(NAXIS1, nfiles)

	lamRange_template = CRVAL1 + [0d, CDELT1 * (NAXIS1 - 1d)]
        log_rebin, lamRange_template, v1, log_temp_template, $
		logLam_template, velscale=velscale

;; ****************************************************************
;; NB: shouldn't this be 0.9A as this is resolution?
        FWHM_tem = 2.5     ; Miles spectra have a resolution
                           ; FWHM of 2.5A.


;; Which templates to use are given in use_templates.pro. This is
;; transfered to the array templatesToUse.
	use_templates, templatesToUse
	nfiles = N_ELEMENTS(templatesToUse)
	templates = MAKE_ARRAY(n_elements(log_temp_template), nfiles)

         
;; Reading the contents of the files into the array templates. 
;; Including rebinning them.
for i = 0, nfiles - 1 do begin


	READCOL, templateFiles[templatesToUse[i]-1], v1,v2, $
		FORMAT = 'D,D', /SILENT


;	READCOL, templateFiles[i], v1,v2, FORMAT = 'D,D', /SILENT

;; Rebinning templates logarthmically
	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

;; Normalizing templates
	templates[*,i] = log_temp_template/median(log_temp_template)
endfor
;-





;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data, header

	spectrum_lin = MAKE_ARRAY(n_elements(galaxy_data[20, 20, *]))
for i = 0, n_elements(spectrum_lin)-1 do begin
	spectrum_lin[i] = galaxy_data[20, 20, i]
endfor

;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/sxpar(header,'CD3_3') ; Sigma difference 
						     ; in pixels

;; For rebinning logarthmically
	lamRange = sxpar(header,'CRVAL3') + $
                   [0,sxpar(header,'CD3_3')*(sxpar(header,'NAXIS3')-1)]
;        lamRange = lamRange/(1+z) ; Compute approximate restframe
        			    ; wavelength range


;; smooth spectrum to fit with templates resolution
	spectrum_lin = gauss_smooth(spectrum_lin, sigma)



;; rebin spectrum logarthmically
	log_rebin, lamrange, spectrum_lin, spectrum_log, $
		logLam_spectrum, velscale=velscale



;; ----------========= Assigning noise variable =============---------
;;   NOISE: vector containing the 1*sigma error (per pixel) in the
;;   galaxy spectrum, or covariance matrix describing the correlated
;;   errors in the galaxy spectrum. Of course this vector/matrix must
;;   have the same units as the galaxy spectrum.
;;   - If GALAXY is a Nx2 array, NOISE has to be an array with the
;;     same dimensions.
;;   - When NOISE has dimensions NxN it is assumed to contain the
;;     covariance matrix with elements sigma(i,j). When the errors in
;;     the spectrum are uncorrelated it is mathematically equivalent
;;     to input in PPXF an error vector NOISE=errvec or a NxN diagonal
;;     matrix NOISE=DIAG_MATRIX(errvec^2) (note squared!).
;;   - IMPORTANT: the penalty term of the pPXF method is based on the
;;     *relative* change of the fit residuals. For this reason the
;;     penalty will work as expected even if no reliable estimate of
;;     the NOISE is available (see Cappellari & Emsellem [2004] for
;;     details).
;;     If no reliable noise is available this keyword can just be set
;;     to:
	noise = MAKE_ARRAY(n_elements(spectrum_log), $
		VALUE = 1d)
		;galaxy*0+1 ; Same weight for all pixels




; The galaxy and the template spectra do not have the same starting
; wavelength. For this reason an extra velocity shift DV has to be
; applied to the template to fit the galaxy spectrum. We remove this
; artificial shift by using the keyword VSYST in the call to PPXF
; below, so that all velocities are measured with respect to DV. This
; assume the redshift is negligible. In the case of a high-redshift
; galaxy one should de-redshift its wavelength to the rest frame
; before using the line below (see above). 

dv = (logLam_template[0]-logLam_spectrum[0])*c ; km/s


; Find the pixels to ignore to avoid being distracted by gas emission
; lines or atmospheric absorbsion line.  
goodPixels = ppxf_determine_goodpixels(logLam_spectrum,$
	lamRange_template,vel) 





        


	start = [vel, sig] ; starting guess

	PPXF, templates, spectrum_log, noise, velscale, start, $
		spaxel_dynamics, BESTFIT = bestfit, $
		GOODPIXELS=goodPixels, MOMENTS = moments, $
		DEGREE = 2, VSYST = dv, WEIGHTS = weights, /QUIET
;;		ERROR = error


;	print, 'Best-fitting redshift z:', (z + 1)*(1 + dynamics[0]/c)- 1

; CALLING SEQUENCE:
;  PPXF, templates, galaxy, noise, velScale, start, sol, BESTFIT=bestFit, $
;	BIAS=bias, CHI2DOF=chi2dof, /CLEAN, COMPONENT=component, $
;	DEGREE=degree, ERROR=error, GOODPIXELS=goodPixels, LAMBDA=lambda, $
;	MDEGREE=mdegree, MOMENTS=moments, MPOLYWEIGHTS=mpolyweights, $
;	/OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, $
;	REDDENING=reddening, REGUL=regul, REG_DIM=reg_dim, SKY=sky, $
;	VSYST=vsyst, WEIGHTS=weights




	residuals = MAKE_ARRAY(n_elements(spectrum_log))
for i = 0, n_elements(spectrum_log)-1 do begin
	residuals[i] = spectrum_log[i] - bestfit[i]
endfor



plot, spectrum_log
oplot, bestfit, COLOR = 100160
oplot, residuals, COLOR = 5090


print, 'v = ', spaxel_dynamics[0]
print, 'sigma = ', spaxel_dynamics[1]
print, 'h3 = ', spaxel_dynamics[2]
print, 'Best-fitting redshift z:', ((1 + spaxel_dynamics[0]/c)/(1 - spaxel_dynamics[0]/c)) - 1
return
end



























;; ==================================================================
;; 		Print the data cube to a table format
;; ==================================================================
;; warrenj 20150330 Process to read the cube format and print it in
;; table form into a text file.

pro print_to_file


	galaxy = 'ngc3557'
	OB = '1'


	galaxyDirectoryArray = FILE_SEARCH('/Data/vimosindi/' + Galaxy + $
		'-' + OB + $
;'/combined/combined_exposures/*.fits')
		'/Final/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits')
	dataCube = galaxyDirectoryArray[0]
	FITS_READ, dataCube, galaxy_data, header

	OPENW, 1, '~/VIMOS_project/analysis/testIO.dat'
;	print, size(galaxy_data)
	printf, 1, 'x       y       lambda       flux'

for i = 0, 39 do begin
for j = 0, 39 do begin
for k = 0, 2799 do begin
	PRINTF, 1, STRING(i) + STRING(j) + STRING(k) + $
		STRING(galaxy_data[i,j,k])
endfor
endfor
endfor
	CLOSE, 1


return
end





