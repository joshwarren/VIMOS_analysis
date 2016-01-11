;; ==================================================================
;; Propergate uncertainty
;; ==================================================================
;; warrenj 20150216 Process to progerate the uncertainty using Monty
;; Carlo methods to get uncertainty in velocity space.


pro errors, i_gal, bin
resolve_routine, ['log_rebin', 'ppxf'];, 'ppxf_determine_goodpixels']
;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
; 	galaxy = galaxies[1]
galaxy = galaxies[i_gal]
	reps = 5000 ;; number of monte carlo reps per bin.
	discard = 2
	range = [4200,10000]
	c = 299792.458d
;  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
;	vel = 114.0d ; Initial estimate of the galaxy velocity and
;	sig = 269.0d ;velocity dispersion in km/s in the rest frame
        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
 
	moments = 4 ; number of componants to calc with ppxf (see 
                    ; keyword moments in ppxf.pro for more details)
	degree = 4 ; order of addative Legendre polynomial used to 
		   ; correct the template continuum shape during the fit 
;; File for output: an array containing the calculated dynamics of the
;; galaxy. 




;dir = '~/'
dir = '/Data/vimosindi/'
;dir2 = '~/'
dir2 = '/Data/idl_libraries/'
	
;; Tessellation input
;	binning_spaxels, galaxy
	tessellation_File = dir + 'analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'



data_file = dir + "analysis/galaxies.txt"
readcol, data_file, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, $
    y_gals, SN_used, skipline=1, format='A,D,D,D,D,D,D', /SILENT

i_gal = where(galaxy_gals eq galaxy)
index=i_gal[0]
vel = vel_gals[index]
sig = sig_gals[index]
z = z_gals[index]



       FWHM_gal = FWHM_gal/(1+z) ; Adjust resolution in Angstrom

;; ----------===============================================---------
;; ----------=============== Run analysis  =================---------
;; ----------===============================================---------


;; ----------=============== Miles library ================---------
; Finding the template files
	templatesDirectory = dir2 + 'ppxf/MILES_library/'
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
;	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

;; Normalizing templates
	templates[*,i] = log_temp_template
endfor
TEMPLATES /= median(log_temp_template)



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
	dataCubeDirectory = FILE_SEARCH(dir + 'reduced/' + Galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 
        
;; For analysis of just one quadrant - mst have used rss2cube_quadrant
;;                                     and have binned the quadrant.
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		galaxy + $
;		'-3/Q2/calibrated/cube/*_fluxcal_cube.fits')



	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 2, /SILENT)

;; write key parameters from header - can then be altered in future	
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CD3_3')
	s = size(galaxy_data_temp)

;; Change to pixel units
IF keyword_set(range) THEN range = FIX((range - CRVAL_spec)/CDELT_spec)

	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
        galaxy_noise = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_data = galaxy_data_temp[discard:s[1]-discard-1, $
		discard:s[2]-discard-1,*]
	galaxy_noise = galaxy_noise_temp[discard:s[1]-discard-1, $
		discard:s[2]-discard-1,*]
;; array to hold results
;	bin_dynamics = MAKE_ARRAY(7, n_bins);, reps)

	n_spaxels = n_elements(galaxy_data[*,0,0]) * $
		n_elements(galaxy_data[0,*,0])

;; ----------========== Spatially Binning =============---------

;; endfor is near the end - after ppxf has been run on this bin.
;for bin=0, n_bins-1 do begin
	spaxels_in_bin = WHERE(bin_num EQ bin, n_spaxels_in_bin)

;; Creates a new spectrum for a new bin.
        bin_lin_temp = MAKE_ARRAY(n_elements(galaxy_data[0,0,*]), $
		VALUE = 0d) 
        bin_lin_noise_temp = MAKE_ARRAY(n_elements(galaxy_noise[0,0,*]), $
		VALUE = 0d) 

for i = 0, n_spaxels_in_bin-1 do begin
	x_i = x[spaxels_in_bin[i]]
	y_i = y[spaxels_in_bin[i]]
for k = 0, s[3]-1 do begin
	bin_lin_temp[k] += galaxy_data[x_i, y_i, k]
        bin_lin_noise_temp[k] += galaxy_noise[x_i, y_i, k]^2
endfor
endfor

	bin_lin_noise_temp = sqrt(bin_lin_noise_temp)
;; bin_lin now contains linearly binned spectrum of the spatial bin.
;; bin_lin_noise contains the errors combined in quadrature. 

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
IF range[0] GT lower_limit THEN lower_limit = range[0]
IF range[1] LT upper_limit THEN upper_limit = range[1]
ENDIF

;lower_limit = 0 
;upper_limit = n_elements(galaxy_data[0,0,*])-1
	lamRange = MAKE_ARRAY(2)
	lamRange[0] = lower_limit*CDELT_spec + CRVAL_spec
	lamRange[1] = upper_limit*CDELT_spec + CRVAL_spec


;; ----------========= Writing the spectrum  =============---------
	bin_lin = MAKE_ARRAY(upper_limit-lower_limit)
	bin_lin_noise = MAKE_ARRAY(upper_limit-lower_limit)
for i = 0, n_elements(bin_lin)-1 do begin
	bin_lin[i] = bin_lin_temp[lower_limit+i]
	bin_lin_noise[i] = bin_lin_noise_temp[lower_limit+i]
endfor

;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/CDELT_temp ; Sigma difference 
						     ; in pixels
;; smooth spectrum to fit with templates resolution
	bin_lin = gauss_smooth(bin_lin, sigma)
        bin_lin_noise = gauss_smooth(bin_lin_noise, sigma)
;;;;;**************should there be any scaling here???*********;;;;;;;;;;;;;;;

	lamRange = lamRange/(1+z)
;; rebin spectrum logarthmically
	log_rebin, lamrange, bin_lin, bin_log, logLam_bin, $
		velscale=velscale
	log_rebin, lamrange, bin_lin_noise^2, bin_log_noise, $
		velscale=velscale
	bin_log_noise = sqrt(bin_log_noise) ;; from log_rebin.pro notes

;; normalise the spectrum
        med_bin = MEDIAN(bin_log)
	bin_log = bin_log/med_bin
        bin_log_noise = bin_log_noise/med_bin

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
;	noise = MAKE_ARRAY(n_elements(bin_log), $
;		VALUE = 1d)
		;galaxy*0+1 ; Same weight for all pixels
noise = bin_log_noise+0.0000000000001

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

;; Run once to get bestfit. Add noise to bestfit so that noise
;; is not 'added twice'.
	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics_temp, BESTFIT = bestfit_sav, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /QUIET, $
		ERROR = errors


bin_output=MAKE_ARRAY(reps,4, /FLOAT)
bin_errors=MAKE_ARRAY(reps,4, /FLOAT)
seed = !NULL
TIC
for rep=0,reps-1 do begin
print, 'rep ', rep
random = randomu(seed, n_elements(noise), /NORMAL)
gaussian = gaussian(random, [1/sqrt(2*!pi),0,1])
add_noise = (random/abs(random))*sqrt((-2*noise^2)*alog(gaussian*noise))
bin_log = bestfit_sav + add_noise

	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics_temp, BIAS = 0.001, $
		GOODPIXELS=goodPixels, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, /QUIET, ERROR = errors


bin_output[rep,0] = bin_dynamics_temp[0]
bin_output[rep,1] = bin_dynamics_temp[1]
bin_output[rep,2] = bin_dynamics_temp[2]
bin_output[rep,3] = bin_dynamics_temp[3]

bin_errors[rep,0] = errors[0]
bin_errors[rep,1] = errors[1]
bin_errors[rep,2] = errors[2]
bin_errors[rep,3] = errors[3]


endfor
TOC
;endfor

FILE_MKDIR, dir + "analysis/" + galaxy + "/errors_results/errors"
bin_file =  dir + "analysis/" + galaxy + "/errors_results/" + $
	STRTRIM(STRING(bin),2) + ".dat"
errors_file =  dir + "analysis/" + galaxy + "/errors_results/errors/" + $
	STRTRIM(STRING(bin),2) + ".dat"

CLOSE,1
OPENW, 1, bin_file
CLOSE, 1

forprint, bin_output[*,0], bin_output[*,1], bin_output[*,2], $
          bin_output[*,3], TEXTOUT = bin_file, /SILENT, /NOCOMMENT

forprint, bin_errors[*,0], bin_errors[*,1], bin_errors[*,2], $
          bin_errors[*,3], TEXTOUT = errors_file, /SILENT, /NOCOMMENT




return
end










