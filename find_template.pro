;; ==================================================================
;; Analyse reduced VIMOS data using pPFX
;; ==================================================================
;; warrenj 20150216 Process to analyse the reduced VIMOS data.

pro find_template, galaxy, z=z, discard=discard, range=range 

;; ----------===============================================---------
;; ----------============ Default parameters ===============---------
;; ----------===============================================---------

if not keyword_set(z) then z=0.01
if not keyword_set(discard) then discard=2
if not keyword_set(range) then range=[4200,10000]

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
;  	galaxy = 'ngc3557'
;	galaxy = 'ic1459'
;	discard = 2
;	range = [4200,10000]
	c = 299792.458d
; 	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
	vel = 0.0d ; Initial estimate of the galaxy velocity and
	sig = 200.0d ;velocity dispersion in km/s in the rest frame
        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
        FWHM_gal = FWHM_gal/(1+z) ; Adjust resolution in Angstrom
	moments = 4 ; number of componants to calc with ppxf (see 
                    ; keyword moments in ppxf.pro for more details)
	degree = 4 ; order of addative Legendre polynomial used to 
		   ; correct the template continuum shape during the fit 
;; File for output: an array containing the calculated dynamics of the
;; galaxy. 
	output_temp_weighting = '/Data/vimos/analysis/' + $
		galaxy + '/templates.txt'

	CLOSE, 1
	OPENW, 1, output_temp_weighting

;; Tessellation input
;	binning_spaxels, galaxy
	tessellation_File = '/Data/vimos/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'


;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------


;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/idl_libraries/ppxf/MILES_library/'
	templateFiles = FILE_SEARCH(templatesDirectory + $
		'm0[0-9][0-9][0-9]V', COUNT=nfiles)

;v1 is wavelength, v2 is spectrum
	READCOL, templateFiles[0], v1,v2, FORMAT = 'D,D', /SILENT

; Using same keywords as fits headers
	CRVAL_temp = v1[0]		; starting wavelength
	NAXIS_temp = size(v2, /N_ELEMENTS) ; Number of entries
	; wavelength increments (resolution?)
	CDELT_temp = (v1[NAXIS_temp-1]-v1[0])/(NAXIS_temp-1)

; Creating the templates array with correct dimension
;	temp_template = MAKE_ARRAY(NAXIS1, nfiles)

	lamRange_template = CRVAL_temp + [0d, CDELT_temp*(NAXIS_temp-1d)]
        log_rebin, lamRange_template, v1, log_temp_template, $
		logLam_template, velscale=velscale

;; ****************************************************************
;; NB: shouldn't this be 0.9A as this is resolution?
        FWHM_tem = 2.5     ; Miles spectra have a resolution
                           ; FWHM of 2.5A.

	templates = MAKE_ARRAY(n_elements(log_temp_template), nfiles)

         
;; Reading the contents of the files into the array templates. 
;; Including rebinning them.
for i = 0, nfiles - 1 do begin

	READCOL, templateFiles[i], v1,v2, FORMAT = 'D,D', /SILENT

;; Rebinning templates logarthmically
;	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

;; Normalizing templates
	templates[*,i] = log_temp_template/median(log_temp_template)
endfor



;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
; Final wildcard reflects the fact that depending on reduction method
; quadrants may or may not have beenflux calibrated.
	dataCubeDirectory = FILE_SEARCH('/Data/vimos/cubes/' + $
		Galaxy + '.cube.combined.fits') 


	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header

;; write key parameters from header - can then be altered in future	
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CDELT3')
	s = size(galaxy_data_temp)

;; Change to pixel units
IF keyword_set(range) THEN range = FIX((range - CRVAL_spec)/CDELT_spec)

	gal_temp = total(total(galaxy_data_temp[discard:s[1]-discard-1, discard:s[2]-discard-1,*], 1), 1)


;; --------======== Finding limits of the spectrum ========--------
;; limits are the cuts in pixel units, while lamRange is the cuts in
;; wavelength unis.
	gap=12
	ignore = FIX((5581 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap  
	ignore2 =FIX((5199 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap 


;; h is the spectrum with the peak enclosed by 'ignore' removed.
	if 5581 lt CRVAL_spec+s[3]*CDELT_spec then begin
		h =[gal_temp[0:ignore[0]],gal_temp[ignore[1]:*]]
        endif else h = gal_temp

	h =[h[0:ignore2[0]],h[ignore2[1]:*]]


	half = s[3]/2
	a = h/MEDIAN(h) - h[4:*]/MEDIAN(h)
	a[WHERE(~FINITE(a))] = 0
	
;	lower_limit = MIN(WHERE(ABS(a) GT 0.2), MAX=upper_limit)
	lower_limit = MAX(WHERE(ABS(a[0:0.5*half]) GT 0.2))
	upper_limit = MIN(WHERE(ABS(a[1.5*half:*]) GT 0.2))+1.5*half




IF (upper_limit GT ignore2[0]) then upper_limit += gap
IF (upper_limit GT ignore[0]) then upper_limit += gap

IF (lower_limit LT 0) THEN BEGIN
	lower_limit = MIN(WHERE(a[0:half] NE 0)) + 5
	IF (lower_limit LT 0) THEN lower_limit = 0 
ENDIF ELSE lower_limit += 5
IF (upper_limit GT s[3]-1) OR (upper_limit LT half) THEN upper_limit=s[3]-6 $
	ELSE upper_limit += - 5

;; --------=========== Using range variable ===========--------
IF keyword_set(range) THEN BEGIN
IF range[0] GT lower_limit THEN lower_limit = range[0]
IF range[1] LT upper_limit THEN upper_limit = range[1]
ENDIF

;lower_limit = 0 
;upper_limit = n_elements(galaxy_data[0,0,*])-1
	lamRange = MAKE_ARRAY(2)
	lamRange[0] = lower_limit*CDELT_spec + CRVAL_spec
	lamRange[1] = upper_limit*CDELT_spec + CRVAL_spec

;; ----------========= Writing the spectrum  =============---------

	gal = gal_temp[lower_limit:upper_limit-1]

;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/CDELT_temp ; Sigma difference 
						     ; in pixels

;; smooth spectrum to fit with templates resolution
	gal = gauss_smooth(gal, sigma)


	lamRange = lamRange/(1+z)
;; rebin spectrum logarthmically
	log_rebin, lamrange, gal, bin_log, logLam_bin, $
		velscale=velscale

;; normalise the spectrum
	bin_log = bin_log/MEDIAN(bin_log)

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
goodPixels = ppxf_determine_goodpixels(logLam_bin,lamRange_template,vel, z) 




        
	lambda = EXP(logLam_bin)

	start = [vel, sig] ; starting guess

print, "lower limit:", lower_limit, lower_limit*CDELT_spec + CRVAL_spec
print, "upper limit:", upper_limit, upper_limit*CDELT_spec + CRVAL_spec

	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics_temp, BESTFIT = bestfit, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /PLOT;, $
;;		/QUIET, ERROR = error
print, ""
print, ""

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

 
 


;+
;; Write weightings for each template used to file 1.
for k = 0, nfiles-1 do begin
if weights[k] ne 0 then begin
;; Use (uncomment) this line if using limited template files.
;	PRINTF, 1, string(templatesToUse[k]) + ' ' + string(weights[k])
;; Use (uncomment) this line if using full library.
	PRINTF, 1, string(k+1) + ' ' + string(weights[k])
endif
endfor
;-




CLOSE, 1




return
end









