;; ==================================================================
;; 		Check the fit of a given spaxel
;; ==================================================================
;; warrenj 20150602 Routine to plot the best fit against the spectrum
;; and show the residuals, for spaxels 20,20 (zero weighted).
;; warrenj 20150604 Altered to fit for any given spaxel.
;; warrenj 20150624 The routine now automatically cuts the spectrum to
;; size such that ppxf does not try to fit non-real data.

pro fit_spaxel

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'ngc3557'
	discard = 2
	spaxel = [19, 21]
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
	degree = 4 ; order of addative Legendre polynomial used to 
		   ; correct the template continuum shape during the fit  
;; File for output: an array containing the calculated dynamics of the
;; galaxy. 



;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------



;+
;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/idl_libraries/ppxf/MILES_library/'
	templateFiles = FILE_SEARCH(templatesDirectory + $
		'm0[0-9][0-9][0-9]V', COUNT=nfiles)

;v1 is wavelength, v2 is spectrum/flux
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
	log_rebin, lamRange_template, v2, log_temp_template;, $
;		velscale=velscale

;; Normalizing templates
	templates[*,i] = log_temp_template/median(log_temp_template)
endfor
;-




;; For rebinning logarthmically
;	lamRange = sxpar(header,'CRVAL3') + $
;                  [0,sxpar(header,'CD3_3')*(sxpar(header,'NAXIS3')-1)]


;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	s = size(galaxy_data_temp)
	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_data = galaxy_data_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]

;; --------======== Finding limits of the spectrum ========--------
;; limits are the cuts in pixel units, while lamRange is the cuts in
;; wavelength unis.
ignore = FIX((5581 - sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')) + [-1,+1]*12
ignore2 =FIX((5199 - sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')) + [-1,+1]*12

	h =[bin_lin_temp[0:ignore[0]],bin_lin_temp[ignore[1]:*]]
;	lower_limit = MIN(WHERE(h/MEDIAN(h) GT 0.6), MAX=upper_limit)
	h =[h[0:ignore2[0]],h[ignore2[1]:*]]

;; --------======= Finding limits of the spectrum 2 =======--------
;; a is an array containing the difference between the ith element and
;; the (i+1)th element of h
a = MAKE_ARRAY(n_elements(h)-2)
half = s[3]/2
for i=0, n_elements(h)-3 do begin
	a[i]=h[i]/MEDIAN(h)-h[i+2]/MEDIAN(h)
	if (FINITE(a[i]) NE 1) THEN a[i]=0
endfor
print, "here", WHERE(ABS(a) GT 0.3)
;oplot, a[100:200]+1, color = 10000

;	lower_limit = MIN(WHERE(ABS(a) GT 0.3), MAX=upper_limit)
	lower_limit = MAX(WHERE(ABS(a[0:half]) GT 0.3))
	upper_limit = MIN(WHERE(ABS(a[half :*]) GT 0.3)+half)




IF (upper_limit GT ignore2[0]) then upper_limit += (ignore2[1]-ignore2[0])
IF (upper_limit GT ignore[0]) then upper_limit += (ignore[1]-ignore[0])
	lower_limit = lower_limit + 5
	upper_limit = upper_limit - 5

IF (lower_limit LT 0) THEN lower_limit = 0
IF (upper_limit GT half*2-1) OR (upper_limit LT 0) THEN upper_limit=half*2-1


;lower_limit = 0 
;upper_limit = n_elements(galaxy_data[0,0,*])-1
print, "lower limit:", lower_limit, $
	lower_limit*sxpar(header,'CD3_3') + sxpar(header,'CRVAL3')
print, "upper limit:", upper_limit, $
	upper_limit*sxpar(header,'CD3_3') + sxpar(header,'CRVAL3')


	lamRange = MAKE_ARRAY(2)
	lamRange[0] = lower_limit*sxpar(header,'CD3_3') + $
		sxpar(header,'CRVAL3')
	lamRange[1] = upper_limit*sxpar(header,'CD3_3') + $
		sxpar(header,'CRVAL3')



;; ----------========= Writing the spectrum  =============---------
	spectrum_lin = MAKE_ARRAY(upper_limit-lower_limit)
for i = 0, n_elements(spectrum_lin)-1 do begin
	spectrum_lin[i] = galaxy_data[spaxel[0], spaxel[1], lower_limit+i]
endfor




;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/sxpar(header,'CD3_3') ; Sigma difference 
						     ; in pixels
;; For rebinning logarthmically
;	lamRange = sxpar(header,'CRVAL3') + $
;                   [0,sxpar(header,'CD3_3')*(sxpar(header,'NAXIS3')-1)]
;        lamRange = lamRange/(1+z) ; Compute approximate restframe
        			    ; wavelength range


;; smooth spectrum to fit with templates resolution
	spectrum_lin = gauss_smooth(spectrum_lin, sigma)

;; rebin spectrum logarthmically
	log_rebin, lamrange, spectrum_lin, spectrum_log, $
		logLam_spectrum, velscale=velscale

;; normalise the spectrum
	spectrum_log = spectrum_log/MEDIAN(spectrum_log)




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


;; Attempt at a more sophisticated approach to handling errors/noise 
;+
;	galaxy_noise = MRDFITS(dataCubeDirectory[0], 2, /SILENT)
;	noise_lin = MAKE_ARRAY(upper_limit-lower_limit)
;for i = 0, n_elements(noise_lin)-1 do begin
;	noise_lin[i] = galaxy_noise[spaxel[0], spaxel[1], lower_limit+i]*$
;		galaxy_noise[spaxel[0], spaxel[1], lower_limit+i]
;endfor
;
;;; smooth spectrum to fit with templates resolution
; 	noise_lin = gauss_smooth(spectrum_lin, sigma)
;
;;; rebin spectrum logarthmically
;	log_rebin, lamrange, noise_lin, noise
;
;;; normalise the spectrum
;	noise = (SQRT(noise))/MEDIAN(spectrum_log)
;-


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







	lambda =EXP(logLam_spectrum)


	start = [vel, sig] ; starting guess

print, 'spaxel: [', spaxel[0], ',', spaxel[1], ']'
	PPXF, templates, spectrum_log, noise, velscale, start, $
		spaxel_dynamics, BESTFIT = bestfit, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /PLOT
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



;+
;	residuals = MAKE_ARRAY(n_elements(spectrum_log))
;for i = 0, n_elements(spectrum_log)-1 do begin
;	residuals[i] = spectrum_log[i] - bestfit[i]
;endfor
;
;
;
;;plot, spectrum_log
;;oplot, bestfit, COLOR = 100160
;;oplot, residuals, COLOR = 5090
;
;
;print, 'v = ', spaxel_dynamics[0]
;print, 'sigma = ', spaxel_dynamics[1]
;print, 'h3 = ', spaxel_dynamics[2]
;print, 'Best-fitting redshift z:', ((1 + spaxel_dynamics[0]/c)/(1 - spaxel_dynamics[0]/c)) - 1
;-
return
end




