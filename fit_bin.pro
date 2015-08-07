;; ==================================================================
;; 		Check the fit of a given bin
;; ==================================================================
;; warrenj 20150604 Routine to plot the best fit against the spectrum
;; for a given bin.

pro fit_bin

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'ngc3557'
	discard = 2
	fit_bin_num = 50
	c = 299792.458d
;  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
	vel = 3000.0d ; Initial estimate of the galaxy velocity in km/s
	sig = 270.0d ; Initial estimate of the galaxy dispersion in km/s 
		     ; within its rest frame
;	spaxel = [20,20] ; spaxel to read off
        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
 ;       FWHM_gal = FWHM_gal/(1+z) ; Adjust resolution in Angstrom
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
	tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'
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

	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	s = size(galaxy_data_temp)
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
for k = 0, n_elements(galaxy_data[x_i,y_i,*])-1 do begin 
	bin_lin_temp[k] = bin_lin_temp[k] + galaxy_data[x_i, y_i, k]
endfor
endfor
;; bin_lin now contains linearly binned spectrum of the spatial bin. 


;; --------======== Finding limits of the spectrum ========--------
;; limits are the cuts in pixel units, while lamRange is the cuts in
;; wavelength unis.
ignore = FIX((5581 - sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')) + [-1,+1]*12

ignore2 =FIX((5199 - sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')) + [-1,+1]*12


;; h is the spectrum with the peak enclosed by 'ignore' and ignore2
;; removed. 
	h =[bin_lin_temp[0:ignore[0]],bin_lin_temp[ignore[1]:*]]

;	lower_limit = MIN(WHERE(h/MEDIAN(h) GT 1.2), MAX=upper_limit)
	h =[h[0:ignore2[0]],h[ignore2[1]:*]]



;; --------======= Finding limits of the spectrum 2 =======--------
;; a is an array containing the difference between the ith element and
;; the (i+#)th element of h
;a = MAKE_ARRAY(n_elements(h)-4)
half = s[3]/2
;plot, h[0:half]/MEDIAN(h)
a = h/MEDIAN(h) - h[4:*]/MEDIAN(h)
a[WHERE(~FINITE(a))] = 0
;for i=0, n_elements(h)-5 do begin
;	a[i]=h[i]/MEDIAN(h)-h[i+4]/MEDIAN(h)
;	if (FINITE(a[i]) NE 1) THEN a[i]=0
;endfor

;print, "here", WHERE(ABS(a) GT 0.2)
;oplot, ABS(a[0:half])+0.5, color = 10000

;	lower_limit = MIN(WHERE(ABS(a) GT 0.3), MAX=upper_limit)
	lower_limit = MAX(WHERE(ABS(a[0:0.5*half]) GT 0.2))
	upper_limit = MIN(WHERE(ABS(a[1.5*half:*]) GT 0.2))+1.5*half

;print, MAX(ABS(a[0:half]))
;print, ABS(a[0:half])

;print, "a", lower_limit, upper_limit

IF (upper_limit GT ignore2[0]) then upper_limit += (ignore2[1]-ignore2[0])
;print, "b", lower_limit, upper_limit
IF (upper_limit GT ignore[0]) then upper_limit += (ignore[1]-ignore[0])

;print, "c", lower_limit, upper_limit

IF (lower_limit LT 0) THEN BEGIN
	lower_limit = MIN(WHERE(a[0:half] NE 0)) + 5
	IF (lower_limit LT 0) THEN lower_limit = 0 
ENDIF ELSE lower_limit += 5
;print, "d", lower_limit, upper_limit
IF (upper_limit GT s[3]-1) OR (upper_limit LT half) THEN upper_limit=s[3]-6 $
	ELSE upper_limit += - 5
;print, "e", lower_limit, upper_limit
	

;upper_limit +=200

;lower_limit = FIX(4600-sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')
;upper_limit = FIX(4700-sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3')
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

;print, bin_lin_temp[upper_limit-20:upper_limit+20]/MEDIAN(bin_lin_temp)









;; ----------========= Writing the spectrum  =============---------
	bin_lin = MAKE_ARRAY(upper_limit-lower_limit)
for i = 0, n_elements(bin_lin)-1 do begin
	bin_lin[i] = bin_lin_temp[lower_limit+i]
endfor


;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/sxpar(header,'CD3_3') ; Sigma difference 
						     ; in pixels
;; smooth spectrum to fit with templates resolution
	bin_lin = gauss_smooth(bin_lin, sigma)



;; rebin spectrum logarthmically
	log_rebin, lamRange, bin_lin, bin_log, logLam_bin, $
		velscale=velscale

;; normalise the bin to the medium value
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
goodPixels = ppxf_determine_goodpixels(logLam_bin,$
	lamRange_template,vel) 





        
	lambda =EXP(logLam_bin)

	start = [vel, sig] ; starting guess

print, 'bin:', fit_bin_num
print, 'spaxels:'
print, 'x = ', x[spaxels_in_bin]
print, 'y = ', y[spaxels_in_bin]

	PPXF, templates, bin_log, noise, velscale, start, $
		spaxel_dynamics, BESTFIT = bestfit, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /PLOT
;;		ERROR = error


	print, 'Best-fitting redshift z:', ((1 + $
		spaxel_dynamics[0]/c)/ (1 - spaxel_dynamics[0]/c)) - 1

; CALLING SEQUENCE:
;  PPXF, templates, galaxy, noise, velScale, start, sol, BESTFIT=bestFit, $
;	BIAS=bias, CHI2DOF=chi2dof, /CLEAN, COMPONENT=component, $
;	DEGREE=degree, ERROR=error, GOODPIXELS=goodPixels, LAMBDA=lambda, $
;	MDEGREE=mdegree, MOMENTS=moments, MPOLYWEIGHTS=mpolyweights, $
;	/OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, $
;	REDDENING=reddening, REGUL=regul, REG_DIM=reg_dim, SKY=sky, $
;	VSYST=vsyst, WEIGHTS=weights



return
end


