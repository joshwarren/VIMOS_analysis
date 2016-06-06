;; ==================================================================
;; Propergate uncertainty
;; ==================================================================
;; warrenj 20150216 Process to progerate the uncertainty using Monty
;; Carlo methods to get uncertainty in velocity space.




function set_lines, lines, logLam_temp, FWHM_gal
;
; In this routine all lines are free to have independent intensities.
; One can fix the intensity ratio of different lines (e.g. the [OIII] doublet)
; by placing them in the same emission template
;
lam = exp(logLam_temp)
lines = lines[where((lines gt min(lam)) and (lines lt max(lam)))]
sigma = FWHM_gal/2.355 ; Assumes instrumental sigma is constant in Angstrom
emission_lines = dblarr(n_elements(logLam_temp),n_elements(lines))
for j=0,n_elements(lines)-1 do $
    emission_lines[*,j] = exp(-0.5d*((lam - lines[j])/sigma)^2)
return, emission_lines
end
;------------------------------------------------------------------------------







;-----------------------------------------------------------------------------
function determine_goodPixels, logLam, lamRangeTemp, vel, z
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


flag or= logLam lt alog(lamRangeTemp[0]) + (vel + 900d)/c ; Mask edges of
flag or= logLam gt alog(lamRangeTemp[1]) + (vel - 900d)/c ; stellar library


flag[0:3] = 1 ; Mask edge of data
flag[-4:*]= 1 ; to remove edge effects
return, where(flag eq 0)
end
;-----------------------------------------------------------------------------







;-----------------------------------------------------------------------------
pro pause
in=' '
READ,"Press enter",in
end
;-----------------------------------------------------------------------------










;-----------------------------------------------------------------------------
pro errors_gas, i_gal, bin
resolve_routine, ['log_rebin', 'ppxf']
;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
galaxy = galaxies[i_gal]
reps = 2 ;; number of monte carlo reps per bin.

if not keyword_set(discard) then discard=2
if not keyword_set(range) then range=[4200,10000]

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

;; ----------======== Files/Directories List ===========---------	
dir = "/Data/vimosindi/"
;dir = "~/"
dir2 = "/Data/idl_libraries/"
;dir2 = "~/"

tessellation_File = dir + 'analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt'

data_file = dir + "analysis/galaxies.txt"

templatesDirectory = dir2 + 'ppxf/MILES_library/'

templateFiles = FILE_SEARCH(templatesDirectory + 'm0[0-9][0-9][0-9]V', $
	COUNT=nstemplates)

emission_File = dir + "analysis/emission_line.dat"

dataCubeDirectory = FILE_SEARCH(dir + 'reduced/' + Galaxy + $
	'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 

;; output files given at the end.




;; ----------======== More input parameters ===========---------	

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

;v1 is wavelength, v2 is spectrum
READCOL, templateFiles[0], v1,v2, FORMAT = 'D,D', /SILENT

; Using same keywords as fits headers
CRVAL_temp = v1[0]		; starting wavelength
NAXIS_temp = size(v2, /N_ELEMENTS) ; Number of entries
; wavelength increments (resolution?)
CDELT_temp = (v1[NAXIS_temp-1]-v1[0])/(NAXIS_temp-1)

; Creating the templates array with correct dimension
;	temp_template = MAKE_ARRAY(NAXIS1, nstemplates)

lamRange_template = CRVAL_temp + [0d, CDELT_temp*(NAXIS_temp-1d)]
log_rebin, lamRange_template, v1, log_temp_template, logLam_template, $
	velscale=velscale

;; ****************************************************************
;; NB: shouldn't this be 0.9A as this is resolution?
FWHM_tem = 2.5     ; Miles spectra have a resolution
                   ; FWHM of 2.5A.


;; Which templates to use are given in use_templates.pro. This is
;; transfered to the array templatesToUse.
use_templates, galaxy, templatesToUse
nstemplates = N_ELEMENTS(templatesToUse)
stellar_templates = MAKE_ARRAY(n_elements(log_temp_template), nstemplates)

         
;; Reading the contents of the files into the array templates. 
;; Including rebinning them.
for i = 0, nstemplates - 1 do begin


	READCOL, templateFiles[templatesToUse[i]-1], v1,v2, $
		FORMAT = 'D,D', /SILENT


;	READCOL, templateFiles[i], v1,v2, FORMAT = 'D,D', /SILENT

;; Rebinning templates logarthmically
	log_rebin, lamRange_template, v2, log_temp_template, $
		velscale=velscale

;; Normalizing templates
	stellar_templates[*,i] = log_temp_template
endfor
stellar_templates /= median(log_temp_template)



;; ----------========= Reading Tessellation  =============---------

;; Reads the txt file containing the output of the binning_spaxels
;; routine. 
	RDFLOAT, tessellation_File, x, y, bin_num, COLUMNS = [1,2,3], $
		SKIPLINE = 1, /SILENT 
	
	n_bins = max(bin_num) + 1
;; Contains the order of the bin numbers in terms of index number.
	order = sort(bin_num)




;; ----------========= Reading the spectrum  =============---------
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
goodPixels = determine_goodpixels(logLam_bin, lamRange_template, vel, z) 

lambda = EXP(logLam_bin)


;; ----------====== Setting up gas templates  ==========---------
; Read in emission-line setup file and fill in emission-line structure
readcol, emission_File, eml_i, eml_name, eml_lambda, eml_action, $
	eml_kind, eml_a, eml_v, eml_s, eml_fit, f='(i,a,f,a,a,f,f,f,a)', $
	skipline=2,comment='#',/silent  

in_range = where((eml_lambda ge lamrange[0]) AND (eml_lambda le lamrange[1]))
eml_i = eml_i[in_range]
eml_name = eml_name[in_range]
eml_lambda = eml_lambda[in_range]



emission_lines = set_lines(eml_lambda, logLam_template, FWHM_gal)

;; ----------========= Combining componants ==============---------
templates = [[stellar_templates], [emission_lines]]

uniq_lines = UNIQ(eml_name, SORT(eml_name))

component = [make_array(nstemplates, value=0),make_array(n_elements(eml_i), $
	/INDEX, /INTEGER)+1]

start1 = [vel, sig] ; starting guess
start = [[start1],[start1]]

for line = 1, n_elements(eml_i)-1 do begin
if eml_name[line] eq eml_name[line-1] then $
	component[line+nstemplates:-1] -= 1 $
else begin
start = [[start],[start1]]
endelse
endfor 


n_components = max(component)+1



;; ----------====== Run to get noiseless bestfit ========---------
PPXF, templates, bin_log, noise, velscale, start, bin_dynamics_sav, $
	BESTFIT = bestfit_sav, GOODPIXELS=goodPixels, LAMBDA=lambda, $
	MOMENTS = moments, DEGREE = degree, VSYST = dv, WEIGHTS = weights, $
	COMPONENT = component, /QUIET



bin_output = MAKE_ARRAY(reps, 5, n_components, /FLOAT)

seed = !NULL
for rep=0,reps-1 do begin
random = randomu(seed, n_elements(noise), /NORMAL)
gaussian = gaussian(random, [1/sqrt(2*!pi),0,1])
add_noise = (random/abs(random))*sqrt((-2*noise^2)*alog(gaussian*noise))
bin_log = bestfit_sav + add_noise

PPXF, templates, bin_log, noise, velscale, start, bin_dynamics_temp, $
;	BESTFIT = bestfit, 
	GOODPIXELS=goodPixels, $; LAMBDA=lambda, $
	MOMENTS = moments,$; DEGREE = degree, 
	VSYST = dv, $
	COMPONENT = component, WEIGHTS = weights, /QUIET

bin_output[rep,1:4,*] = bin_dynamics_temp[0:3,*]

;; weightings:
for comp=0,n_components-1 do begin
	if max(weights[where(component eq comp)]) ne 0 then $
		bin_output[rep,0,comp] = $
			total(weights[where(component eq comp)] gt 0)
endfor ;comp
endfor ;rep

;; ----------=========== Saving the outputs =============---------

;; Creating output files/directories.
output_dir = dir + "analysis/" + galaxy + "/montecarlo/"
FILE_MKDIR, output_dir + "/stellar"
FILE_MKDIR, output_dir + eml_name[uniq_lines]

bin_files = [output_dir + "/stellar/" + STRTRIM(STRING(bin),2) + ".dat"]
bin_files = [bin_files, output_dir + "/" + eml_name[uniq_lines] + "/" + $
	STRTRIM(STRING(bin),2) + ".dat"]

;; Save component dynamics
for output=0,n_elements(bin_files)-1 do begin
;CLOSE, output+1 ; closes if already open
;OPENW, output+1, bin_files[output] ; creates file
;CLOSE, output+1 ; closes such that forprint can reopen


forprint, bin_output[*,0,output], bin_output[*,1,output], $
	bin_output[*,2,output], bin_output[*,3,output], $
	bin_output[*,4,output], TEXTOUT = bin_files[output],  $
	/SILENT, /NOCOMMENT

endfor ;output

;; Save (noiseless) bestfit
FILE_MKDIR, output_dir + "/bestfit/dynamics/"
forprint, bestfit_sav, $
	TEXTOUT= output_dir + "/bestfit/" + STRTRIM(STRING(bin),2) + ".dat", $
	/SILENT, /NOCOMMENT

forprint, bin_dynamics_sav[0:3], TEXTOUT = output_dir + "/bestfit/dynamics/" $
	+ STRTRIM(STRING(bin),2) + ".dat", /SILENT, /NOCOMMENT




return
end










