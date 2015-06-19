;; ==================================================================
;; 		Check the fit of a given spaxel
;; ==================================================================
;; warrenj 20150602 Routine to plot the best fit against the spectrum
;; and show the residuals, for spaxels 20,20 (zero weighted).
;; warrenj 20150604 Altered to fit for any given spaxel.

pro fit_spaxel

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'ngc3557'
	spaxel = [20, 20]
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

lower_cut = 0 	; in pixel units
upper_cut = 0 	; in pixel units


;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------



;+
;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/ppxf/MILES_library/'
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


;; Applying a 
if (lower_cut GT CRVAL1 + CDELT1*(NAXIS1-1d)) then begin
	lower_cut = CRVAL1 + CDELT1*(NAXIS1-1d)
endif
if (lower_cut GT CRVAL1) then begin
for k=0, (lower_cut-CRVAL1)/CDELT1 do begin
	v2[k] = 0
endfor
endif

if (upper_cut LT CRVAL1) then begin
	upper_cut = CRVAL1
endif
if (upper_cut NE 0) then begin
;for k=(upper_cut-CRVAL1)/CDELT1, n_elements(v2)-1 do begin
;	v2[k] = 0
;endfor
endif
;; Rebinning templates logarthmically
	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template;, $
;		velscale=velscale

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

	spectrum_lin = MAKE_ARRAY(n_elements(galaxy_data[spaxel[0], $
		spaxel[1], *]))
for i = 0, n_elements(spectrum_lin)-1 do begin
	spectrum_lin[i] = galaxy_data[spaxel[0], spaxel[1], i]
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


for k=0, (lower_cut-sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3') $
	do begin
spectrum_lin[k] = 0
endfor
if (upper_cut NE 0) then begin
;for k=(upper_cut-sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3'), $
;	n_elements(spectrum_lin)-1 do begin
;spectrum_lin[k] = 0
;endfor
endif

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




