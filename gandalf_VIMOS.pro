;; ==================================================================
;; Analyse reduced VIMOS data using pPFX and Gandalf
;; ==================================================================
;; warrenj 20150216 Process to analyse the reduced VIMOS data.
;; warrenj 20150915 Gandalf code added to analyse hot gas componant.



;-----------------------------------------------------------------------------
function determine_goodPixels, emission, logLam, lamRangeTemp, vel, z
; warrenj 20150905 Copied from ppxf_determine_goodPixels.pro
;
; PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of
;	goodPixels to be used as input keyword for the routine
;	PPXF. This is useful to mask gas emission lines or atmospheric
;	absorptions. It can be trivially adapted to mask different
;	lines. 
; 
; INPUT PARAMETERS:
; - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
;     of each pixel of the log rebinned *galaxy* spectrum.
; - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the
;     minimum and maximum wavelength in Angstrom in the stellar
;     *template* used in PPXF. 
; - VEL: Estimate of the galaxy velocity in km/s.
; 
; V1.0: Michele Cappellari, Leiden, 9 September 2005
; V1.01: Made a separate routine and included additional common
;   emission lines. MC, Oxford 12 January 2012
; V1.02: Included more lines. MC, Oxford, 7 Januray 2014

c = 299792.458d ; speed of light in km/s

;; 20150617 warrenj Added Telluric lines (tell) at 5199 (is a blended sky
;; line)

 
;dv = lines*0+800d ; width/2 of masked gas emission region in km/s
dv = 800d ; width/2 of masked gas emission region in km/s

flag = bytarr(n_elements(logLam))

; Marks telluric line
tell = 5199
flag or= logLam gt alog(tell) - z - dv/c $
     and logLam lt alog(tell) - z + dv/c 

; Mask emission lines passed to routine
FOR i=0, n_elements(emission.i)-1 DO BEGIN
    IF (emission.action[i] EQ 'm') THEN BEGIN
        flag or= logLam gt alog(emission.lambda[i]) + (vel - dv)/c $
             and logLam lt alog(emission.lambda[i]) + (vel + dv)/c
    ENDIF
ENDFOR

flag or= logLam lt alog(lamRangeTemp[0]) + (vel + 900d)/c ; Mask edges of
flag or= logLam gt alog(lamRangeTemp[1]) + (vel - 900d)/c ; stellar library


flag[0:3] = 1 ; Mask edge of data
flag[-4:*]= 1 ; to remove edge effects
return, where(flag eq 0)
end
;-----------------------------------------------------------------------------








pro pause
in=' '
READ,"Press enter",in
end












pro gandalf_VIMOS;, galaxy, discard, limits

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
	quiet = boolean(1) ; 1 = yes = true
;  	galaxy = 'ngc3557'
	galaxy = 'ic1459'
	discard = 2
	range = [4200,10000]
	c = 299792.458d
  	z = 0.01 ; redshift to move galaxy spectrum to its rest frame 
	vel = 114.0d ; Initial estimate of the galaxy velocity and
	sig = 269.0d ;velocity dispersion in km/s in the rest frame
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
output_OIII = '/Data/vimosindi/analysis/' + galaxy + '/results/' + $
    'gal_OIII.dat'
output_NI = '/Data/vimosindi/analysis/' + galaxy + '/results/' + $
    'gal_NI.dat'
output_Hb = '/Data/vimosindi/analysis/' + galaxy + '/results/' + $
    'gal_Hb.dat'
output_Hd = '/Data/vimosindi/analysis/' + galaxy + '/results/' + $
    'gal_Hd.dat'
	
;; Tessellation input
;	binning_spaxels, galaxy
	tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'
;; Emission line file
	emission_File = "/Data/vimosindi/analysis/emission_line.dat"


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
;	lamRange_template = CRVAL1 + [0d, CDELT1*(NAXIS1 - 1d)]
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

	templates[*,i] = log_temp_template
endfor
;; Normalizing templates
templates /= median(templates)



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
; Final wildcard reflects the fact that depending on reduction method
; quadrants may or may not have beenflux calibrated.
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 
        
;; For analysis of just one quadrant - mst have used rss2cube_quadrant
;;                                     and have binned the quadrant.
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		galaxy + $
;		'-3/Q2/calibrated/cube/*_fluxcal_cube.fits')



	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header

;; write key parameters from header - can then be altered in future	
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CD3_3')
	s = size(galaxy_data_temp)

;; Change to pixel units
IF keyword_set(range) THEN range = FIX((range - CRVAL_spec)/CDELT_spec)

	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_data = galaxy_data_temp[discard:s[1]-discard-1, $
		discard:s[2]-discard-1,*]


	n_spaxels = n_elements(galaxy_data[*,0,0]) * $
		n_elements(galaxy_data[0,*,0])

;; ----------========= Arrays to hold results =============---------
stellar_bin_dynamics = MAKE_ARRAY(7, n_bins)
OIII_dynamics = MAKE_ARRAY(n_bins, VALUE = 0.0)
NI_dynamics = MAKE_ARRAY(n_bins, VALUE = 0.0)
Hb_dynamics = MAKE_ARRAY(n_bins, VALUE = 0.0)
Hd_dynamics = MAKE_ARRAY(n_bins, VALUE = 0.0)

;; ----------===== Masking structure for emission lines =====---------
;; Find the pixels to ignore to avoid being distracted by gas emission
;; lines or atmospheric absorbsion line.  
;goodPixels = ppxf_determine_goodpixels(logLam_bin,lamRange_template,vel, z) 


; Read in emission-line setup file and fill in emission-line structure
	readcol, emission_File, eml_i, eml_name, eml_lambda, $
                eml_action, eml_kind, eml_a, eml_v, eml_s, eml_fit, $
		f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent  
; creating emission setup structure, to be used for masking the
; emission-line contaminated regions and fitting the emission lines in
; ppxf_gas
	emission = create_struct('i',eml_i,'name',eml_name,$
		'lambda', eml_lambda, 'action', eml_action, $
		'kind', eml_kind, 'a', eml_a, 'v', eml_v, 's', eml_s, $
		'fit', eml_fit)  
; setting the galaxy systemic velocity and 50km/s as initial guesses
; for the gas velocities and velocity dispersions, over-writing the
; emission-line setup values
	emission.v = vel & emission.s = 50.0d


;; ----------========== Spatially Binning =============---------

;; endfor is near the end - after ppxf has been run on this bin.
for bin=0, n_bins-1 do begin
;bin = 1 ;;;;;************************************************
	spaxels_in_bin = WHERE(bin_num EQ bin, n_spaxels_in_bin)


;; Creates a new spectrum for a new bin.
        bin_lin_temp = MAKE_ARRAY(n_elements(galaxy_data[0,0,*]), $
		VALUE = 0d) 

for i = 0, n_spaxels_in_bin-1 do begin
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
for i = 0, n_elements(bin_lin)-1 do begin
	bin_lin[i] = bin_lin_temp[lower_limit+i]
endfor

;; ----------======== Calibrating the spectrum  ===========---------
;; For calibrating the resolutions between templates and observations
;; using the gauss_smooth command
	FWHM_dif = SQRT(FWHM_tem^2 - FWHM_gal^2)
	sigma = FWHM_dif/2.355/CDELT_temp ; Sigma difference 
						     ; in pixels

;; smooth spectrum to fit with templates resolution
	bin_lin = gauss_smooth(bin_lin, sigma)


	lamRange = lamRange/(1+z)
;; rebin spectrum logarthmically
	log_rebin, lamrange, bin_lin, bin_log, logLam_bin, $
		velscale=velscale
	lambda = EXP(logLam_bin)
	log_gal_start = logLam_bin[0]
	log_gal_step = logLam_bin[1]-logLam_bin[0]

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
; assume the redshift is negligible. 

	dv = (logLam_template[0]-logLam_bin[0])*c ; km/s

; Initial mask for pPXF with all emission lines masked.
	goodPixels = determine_goodPixels(emission, logLam_bin, $
		lamRange_template, vel, z)    

	start = [vel, sig] ; starting guess

;; ----------======= Stellar spectrum fitting (pPXF) ========--------- 

print, "bin:", bin, "/", FIX(n_bins-1)
print, "lower limit:", lower_limit, lower_limit*CDELT_spec + CRVAL_spec
print, "upper limit:", upper_limit, upper_limit*CDELT_spec + CRVAL_spec
print, "spaxels:"
print, 'x = ', x[spaxels_in_bin]
print, 'y = ', y[spaxels_in_bin]

	PPXF, templates, bin_log, noise, velscale, start, $
		bin_dynamics_temp, BESTFIT = Pbestfit, $
		GOODPIXELS=goodPixels, LAMBDA=lambda, MOMENTS = moments, $
		DEGREE = degree, VSYST = dv, WEIGHTS = weights, /PLOT, $
		QUIET = quiet;, ERROR = error

;; Add systematic velocity to the dynamics array.
bin_dynamics_temp[0] += dv
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


;; ----------======== Lift mask on [OIII] and Hb ===========--------- 
i_HbO3 = where(emission.name eq 'Hb' or emission.name eq '[OIII]')
emission.action[i_HbO3] = 'f'
goodPixels = determine_goodPixels(emission, logLam_bin, lamRange_template, $
    vel, z) 

if not quiet then PAUSE

;; ----------======= Gas spectrum fitting 1 (Gandalf) ======--------- 
sol = bin_dynamics_temp
GANDALF, templates, bin_log, noise, velscale, sol, $;bin_dynamics_temp, $
    emission, log_gal_start, log_gal_step, GOODPIXELS = goodPixels, $
    DEGREE = degree,  BESTFIT = Gbestfit, WEIGHTS = weights, /PLOT, $
    QUIET = quiet 

; CALLING SEQUENCE:
; PRO GANDALF, templates, galaxy, noise, velScale, sol, $
; 	emission_setup, l0_gal, lstep_gal, GOODPIXELS=goodPixels, $
;	DEGREE=degree, MDEGREE=mdegree, INT_DISP=int_disp, $
;	BESTFIT=bestFit, EMISSION_TEMPLATES=emission_templates, $
;	WEIGHTS=weights, ERROR=esol, PLOT=plot, QUIET=quiet, $
;	LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ,$
;	FOR_ERRORS=for_errors 
;

if not quiet then PAUSE

;; ----------============ Lift remaining masks ============--------- 
i_mlines = where(emission.action eq 'm')
emission.action[i_mlines] = 'f'
goodPixels = determine_goodPixels(emission, logLam_bin, lamRange_template, $
    vel, z) 
; fix Hb and [OIII] kinematics
Hb = where(emission.name eq 'Hb')
O3 = where(emission.name eq '[OIII]')
emission.v[Hb] = sol[2]
emission.s[Hb] = sol[3]
emission.v[O3] = sol[6]
emission.s[O3] = sol[7]
emission.fit[i_Hbo3]  = 'h'
; fix the [NI] and Hg kinematics to that of Hb
i_n1Hg = where(emission.name eq '[NI]' and emission.name eq 'Hg')
emission.v[i_n1Hg]       = sol[2]
emission.s[i_n1Hg]       = sol[3]
;emission.action[i_n1Hg] = 'f'
goodPixels = determine_goodPixels(emission, logLam_bin, lamRange_template, $
    vel, z)

emission.fit[i_n1Hg] = 'f'

;; ----------======= Gas spectrum fitting 2 (Gandalf) ======--------- 
sol = bin_dynamics_temp
GANDALF, templates, bin_log, noise, velscale, sol, $;bin_dynamics_temp, $
    emission, log_gal_start, log_gal_step, GOODPIXELS = goodPixels, $
    DEGREE = degree,  BESTFIT = Gbestfit, WEIGHTS = weights, /PLOT, $
    QUIET = quiet 



;; ----------================ Saving results ===============--------- 
stellar_bin_dynamics[0,bin]=bin_dynamics_temp[0]
stellar_bin_dynamics[1,bin]=bin_dynamics_temp[1]
stellar_bin_dynamics[2,bin]=bin_dynamics_temp[2]
stellar_bin_dynamics[3,bin]=bin_dynamics_temp[3]
stellar_bin_dynamics[4,bin]=bin_dynamics_temp[4]
stellar_bin_dynamics[5,bin]=bin_dynamics_temp[5]
stellar_bin_dynamics[6,bin]=bin_dynamics_temp[6]
 
if sol[20] ne 0 then OIII_dynamics[bin] = sol[22]
if sol[16] ne 0 then Hb_dynamics[bin] = sol[18]
if sol[24] ne 0 then NI_dynamics[bin] = sol[26]
if sol[8] ne 0 then Hd_dynamics[bin] = sol[10]





;; ----------======= Reset emission line structure =====--------- 
    emission.action = 'm'
    emission.v      = vel 
    emission.s      = 50.0d
    emission.fit[i_Hbo3] = 'f'
    emission.fit[i_n1Hg] = 'h'
endfor ;;;;;************************************************


;; Error check - making sure all spaxels have been read into some
;; bin. 
if (i EQ n_spaxels-1) THEN BEGIN
	print, 'ERROR: not all spaxels have been read'
endif 

;; ----------========== Write results to file ==========--------- 
;; Open and print to files
	CLOSE, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
	OPENW, 1, output_temp_weighting
	OPENW, 2, output_v
	OPENW, 3, output_sigma
	OPENW, 4, output_h3
	OPENW, 5, output_h4
	OPENW, 6, output_h5
	OPENW, 7, output_h6
	OPENW, 8, output_Chi
	OPENW, 9, output_OIII
	OPENW, 10, output_Hb
	OPENW, 11, output_NI
	OPENW, 12, output_Hd
for bin=0, n_bins-1 do begin
	PRINTF, 2, stellar_bin_dynamics[0,bin]
	PRINTF, 3, stellar_bin_dynamics[1,bin]
	PRINTF, 4, stellar_bin_dynamics[2,bin]
	PRINTF, 5, stellar_bin_dynamics[3,bin]
;	PRINTF, 6, stellar_bin_dynamics[4,bin]
;	PRINTF, 7, stellar_bin_dynamics[5,bin]
	PRINTF, 8, stellar_bin_dynamics[6,bin]
	PRINTF, 9, OIII_dynamics[bin]
	PRINTF, 10, Hb_dynamics[bin]
	PRINTF, 11, NI_dynamics[bin]
	PRINTF, 12, Hd_dynamics[bin]
endfor


CLOSE, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12













return
end
