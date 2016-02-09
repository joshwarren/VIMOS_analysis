;; ==================================================================
;; 		Check the fit of a given bin
;; ==================================================================
;; warrenj 20150604 Routine to plot the best fit against the spectrum
;; for a given bin.

pro mcmc_fit_bin, galaxy, z, vel, sig, chi

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
;  	galaxy = 'ngc3557'
	discard = 2
	range = [4000,4500]
	c = 299792.458d
;  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
;	vel = 100.0d ; Initial estimate of the galaxy velocity in km/s
;	sig = 200.0d ; Initial estimate of the galaxy dispersion in km/s 
		     ; within its rest frame
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




;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------



;+
;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = '/Data/idl_libraries/ppxf/MILES_library/'
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
	use_templates, galaxy, templatesToUse
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
	tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'
	RDFLOAT, tessellation_File, x, y, bin_num, COLUMNS = [1,2,3], $
		SKIPLINE = 1, /SILENT 
	
	n_bins = max(bin_num) + 1
;; Contains the order of the bin numbers in terms of index number.
	order = sort(bin_num)



	fit_bin_num = FIX(n_bins/2)


;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	s = size(galaxy_data_temp)
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CD3_3')
	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_data = galaxy_data_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]

;+
;; normalise the entire cube
;	galaxy_data = galaxy_data/MEDIAN(galaxy_data)
;-
	n_spaxels = n_elements(galaxy_data[*,0,0]) * $
		n_elements(galaxy_data[0,*,0])


;; ----------========== Spatially Binning =============---------


; spaxels in the given bin
spaxels_in_bin = WHERE(bin_num EQ fit_bin_num, n_spaxels)


;; Need to create a new spectrum for a new bin.
bin_lin_temp = MAKE_ARRAY(n_elements(galaxy_data[0,0,*]), VALUE = 0d)

for i = 0, n_spaxels-1 do begin
	x_i = x[spaxels_in_bin[i]]
	y_i = y[spaxels_in_bin[i]]
for k = 0, s[3]-1 do begin 
	bin_lin_temp[k] += galaxy_data[x_i, y_i, k]
endfor
endfor
;; bin_lin now contains linearly binned spectrum of the spatial bin. 


;; --------======== Finding limits of the spectrum ========--------
;; limits are the cuts in pixel units, while lamRange is the cuts in
;; wavelength unis.
	gap=12
	ignore = FIX((5581 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap  
	ignore2 =FIX((5199 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap 


;; h is the spectrum with the peak enclosed by 'ignore' removed.
	h =[bin_lin_temp[0:ignore[0]],bin_lin_temp[ignore[1]:*]]

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
;; Change to pixel units
	range = FIX((range - CRVAL_spec)/CDELT_spec)
IF range[0] GT lower_limit THEN lower_limit = range[0]
IF range[1] LT upper_limit THEN upper_limit = range[1]
ENDIF



;print, "lower limit:", lower_limit, $
;	lower_limit*sxpar(header,'CD3_3') + sxpar(header,'CRVAL3')
;print, "upper limit:", upper_limit, $
;	upper_limit*sxpar(header,'CD3_3') + sxpar(header,'CRVAL3')


	lamRange = MAKE_ARRAY(2)
	lamRange[0] = lower_limit*sxpar(header,'CD3_3') + $
		sxpar(header,'CRVAL3')
	lamRange[1] = upper_limit*sxpar(header,'CD3_3') + $
		sxpar(header,'CRVAL3')

;print, bin_lin_temp[upper_limit-20:upper_limit+20]/MEDIAN(bin_lin_temp)







;; ----------========= Writing the spectrum  =============---------
	bin_lin = MAKE_ARRAY(upper_limit-lower_limit)
for i = 0, n_elements(bin_lin)-1 do begin
	bin_lin[i] = bin_lin_temp[lower_limit+i]
endfor
;for i=0,n_elements(bin_lin)-1 do $
;print, i, i*sxpar(header,'CD3_3')+lamRange[0], bin_lin[i]/MEDIAN(bin_lin)


;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/sxpar(header,'CD3_3') ; Sigma difference 
						     ; in pixels
;; smooth spectrum to fit with templates resolution
	bin_lin = gauss_smooth(bin_lin, sigma)


lamRange = lamRange/(1+z)
;; rebin spectrum logarthmically
	log_rebin, lamRange, bin_lin, bin_log, logLam_bin, $
		velscale=velscale

;; normalise the bin to the medium value
	bin_log = bin_log/MEDIAN(bin_log)
;for i=0,n_elements(bin_log)-1 do print, i, exp(logLam_bin[i]), bin_log[i]

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
goodPixels = ppxf_determine_goodpixels(logLam_bin,$
	lamRange_template,vel, z) 



        
	lambda =EXP(logLam_bin)

	start = [vel, sig] ; starting guess

	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics, BESTFIT = bestfit, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /QUIET, /PLOT
;;		ERROR = error


; CALLING SEQUENCE:
;  PPXF, templates, galaxy, noise, velScale, start, sol, BESTFIT=bestFit, $
;	BIAS=bias, CHI2DOF=chi2dof, /CLEAN, COMPONENT=component, $
;	DEGREE=degree, ERROR=error, GOODPIXELS=goodPixels, LAMBDA=lambda, $
;	MDEGREE=mdegree, MOMENTS=moments, MPOLYWEIGHTS=mpolyweights, $
;	/OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, $
;	REDDENING=reddening, REGUL=regul, REG_DIM=reg_dim, SKY=sky, $
;	VSYST=vsyst, WEIGHTS=weights

	vel = bin_dynamics[0]
	sig = bin_dynamics[1]
	chi = bin_dynamics[6]


return
end


