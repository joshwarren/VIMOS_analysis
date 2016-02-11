;; ==================================================================
;; 		MCMC to find systematic v and sigma
;; ==================================================================
;; warrenj 20160210 Changed to fit the whole galaxy spectrum and have
;; combined mcmc and mcmc_fit_bin into one routine. 

pro mcmc, galaxy, z=z, vel=vel, sig=sig, discard=discard, range=range

chi = 1000000
repeats = 100

results = MAKE_ARRAY(2,repeats-2)

;; ----------===============================================---------
;; ----------============ Default parameters ===============---------
;; ----------===============================================---------
if not keyword_set(z) then z=0.01
if not keyword_set(discard) then discard=2
if not keyword_set(range) then range=[4200,10000]
if not keyword_set(vel) then vel=0.0d
if not keyword_set(sigma) then sig=200.0d

;; ----------===============================================---------
;; ----------============= Input parameters ================---------
;; ----------===============================================---------
	quiet = boolean(1) ; 1 = yes = true
	plot = boolean(1)
	c = 299792.458d

        FWHM_gal = 4*0.571 ; The fibre FWHM on VIMOS is
                           ; about 4px with a dispersion of
                           ; 0.571A/px. (From: http://www.eso.org
                           ; /sci/facilities/paranal/instruments
                           ; /vimos/inst/ifu.html)
	moments = 4 ; number of componants to calc with ppxf (see 
                    ; keyword moments in ppxf.pro for more details)
	degree = 4 ; order of addative Legendre polynomial used to 
		   ; correct the template continuum shape during the fit 


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

	templates[*,i] = log_temp_template
endfor
;; Normalizing templates
templates /= median(templates)


;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
; Final wildcard reflects the fact that depending on reduction method
; quadrants may or may not have beenflux calibrated.
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 
        


	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 2, /SILENT)

;; write key parameters from header - can then be altered in future	
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CD3_3')
	s = size(galaxy_data_temp)

;; Change to pixel units
IF keyword_set(range) THEN range = FIX((range - CRVAL_spec)/CDELT_spec)

galaxy_data = total(total(galaxy_data_temp[discard:s[1]-discard-1, $
	discard:s[2]-discard-1,*], 1), 1)
;; Summed in quadrature
galaxy_noise = SQRT(total(total(galaxy_noise_temp[discard:s[1]-discard-1, $ 
	discard:s[2]-discard-1,*]^2, 1), 1))


;; --------======== Finding limits of the spectrum ========--------
;; limits are the cuts in pixel units, while lamRange is the cuts in
;; wavelength unis.
	gap=12
	ignore = FIX((5581 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap  
	ignore2 =FIX((5199 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap 


;; h is the spectrum with the peak enclosed by 'ignore' removed.
	h =[galaxy_data[0:ignore[0]],galaxy_data[ignore[1]:*]]

	h =[h[0:ignore2[0]],h[ignore2[1]:*]]

	half = s[3]/2
	a = h/MEDIAN(h) - h[4:*]/MEDIAN(h)
	a[WHERE(~FINITE(a))] = 0
	
;	lower_limit = MIN(WHERE(ABS(a) GT 0.2), MAX=upper_limit)
	lower_limit = MAX(WHERE(ABS(a[0:0.5*half]) GT 0.2))
	upper_limit = MIN(WHERE(ABS(a[1.5*half:*]) GT 0.15))+1.5*half


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
print, lamRange
;; ----------========= Writing the spectrum  =============---------
gal_lin = galaxy_data[lower_limit:upper_limit-1]
gal_lin_noise = galaxy_noise[lower_limit:upper_limit-1]









for i = 0, repeats-1 do begin
print, "i=",i
    v_sav = vel
    sigma_sav = sig
    chi_sav = chi
print, "input v: ",v_sav
print, "input sig: ",sigma_sav
print, "input chi2: ",chi_sav
;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
    FWHM_gal_MC = FWHM_gal/(1+z) ; Adjust resolution in Angstrom
    FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal_MC^2)
    sigma = FWHM_dif/2.355/CDELT_temp ; Sigma difference in pixels

;; smooth spectrum to fit with templates resolution
    gal_lin_MC = gauss_smooth(gal_lin, sigma)
    gal_lin_noise_MC = gauss_smooth(gal_lin_noise, sigma)
    lamRange_MC = lamRange/(1+z)
;; rebin spectrum logarthmically
    log_rebin, lamrange_MC, gal_lin_MC, gal_log, logLam_bin, $
	velscale=velscale

    log_rebin, lamrange_MC, gal_lin_noise_MC^2, gal_log_noise, $
	velscale=velscale
    gal_log_noise = sqrt(gal_log_noise) ;; from log_rebin.pro notes
    lambda = EXP(logLam_bin)
    log_gal_start = logLam_bin[0]
    log_gal_step = logLam_bin[1]-logLam_bin[0]

;; normalise the spectrum
    med_bin = MEDIAN(gal_log)
    gal_log = gal_log/med_bin
    gal_log_noise = gal_log_noise/med_bin

;; ----------========= Assigning noise variable =============---------
    noise = gal_log_noise+0.0000000000001



     
; The galaxy and the template spectra do not have the same starting
; wavelength. For this reason an extra velocity shift DV has to be
; applied to the template to fit the galaxy spectrum. We remove this
; artificial shift by using the keyword VSYST in the call to PPXF
; below, so that all velocities are measured with respect to DV. This
; assume the redshift is negligible. 

    dv = (logLam_template[0]-logLam_bin[0])*c ; km/s
    
;; ----------======= Stellar spectrum fitting (pPXF) ========--------- 
; Initial mask for pPXF with all emission lines masked.
    goodPixels = ppxf_determine_goodPixels(logLam_bin, lamRange_template, $
	vel, z)    

    start = [vel, sig] ; starting guess

    PPXF, templates, gal_log, noise, velscale, start, bin_dynamics, $
	BESTFIT = Pbestfit, GOODPIXELS=goodPixels, LAMBDA=lambda, $
	MOMENTS = moments, DEGREE = degree, VSYST = dv, WEIGHTS = weights, $
	PLOT=plot, QUIET = quiet, ERROR = perror

    vel = bin_dynamics[0]
    sig = bin_dynamics[1]
    chi = bin_dynamics[6]
print, "chi2 from this fit: ",chi

    if i le 1 then begin
        z = (z + 1)*sqrt((1 + vel/c)/(1 - vel/c)) - 1 
        vel = v_sav
        sig = sigma_sav
        chi = chi_sav
    endif else begin
        results[*,i-2]=[vel,sig]
        local_min=randomu(seed0)
        IF ((abs(chi-1) GT abs(chi_sav-1)) AND (i NE repeats-1)) OR $
	    local_min GT 0.85 THEN BEGIN
            step_range=30
	    vel=v_sav + (randomu(seed1)*2*step_range)-step_range
            sig=sigma_sav + (randomu(seed2)*2*step_range)-step_range
	    chi=chi_sav
        ENDIF
    endelse
endfor

vel = MEAN(results[0,*])
sig = MEAN(results[1,*])

print, "Mean vel: :",vel
print, "MEAN vel dispersion: ",sig
show_rst = SCATTERPLOT(results[0,*],results[1,*], XTITLE="velocity", $
    YTITLE="velocity dispersion", TITLE="MCMC for initial conditions")
show_mean = SCATTERPLOT(vel,sig, SYMBOL='star', SYM_SIZE=2.0, /OVERPLOT, $
    /SYM_FILLED)

show_mean.Save, "/Data/vimosindi/analysis/"+galaxy+"/MCMC_inital fit.png", $
    BORDER=10, RESOLUTION=300;, TRANSPARENT=[255, 255, 255] 







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


forprint2, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_used, $
    textout=data_file, /SILENT, WIDTH=90,$
    Comment = "Galaxy      z     velocity     velocity dispersion    x    y    Target SN"


return
end


