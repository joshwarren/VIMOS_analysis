;; ==================================================================
;; 		Check the fit of a given spaxel
;; ==================================================================
;; warrenj 20150602 Routine to plot the best fit against the spectrum
;; and show the residuals, for spaxels 20,20 (zero weighted).
;; warrenj 20150604 Altered to fit for any given spaxel.
;; warrenj 20150624 The routine now automatically cuts the spectrum to
;; size such that ppxf does not try to fit non-real data.

pro find_skyline

;; ----------===============================================---------
;; ----------============= Input parameters  ===============---------
;; ----------===============================================---------
  	galaxy = 'eso443-g024'
	discard = 2
;	spaxel = [19, 21]
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




;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
	s = size(galaxy_data_temp)
	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_data = galaxy_data_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]


;galaxy_data *=10^15


sky = 5199
sky_px = FIX((sky-sxpar(header,'CRVAL3'))/sxpar(header,'CD3_3'))

l = where_XYZ(galaxy_data[*,*,sky_px] gt 9/7 * galaxy_data[*,*,sky_px-6], xind=xind, yind=yind)
;for i=0, n_elements(xind)-1 do print, xind[i], yind[i]
print, size(xind)

a = [[xind],[yind]]
b = galaxy_data[*,*,sky_px]
t= b[a]
print, size(t)
m = max(t,mind)
print, mind, m
print, t[mind], xind[mind], yind[mind]
print, size(xind)
print, size(galaxy_data)


e = where_XYZ(b eq m, xind = xind, yind=yind)
print, "*****************************sdfsfsdfsd****************"
for i=0, n_elements(xind)-1 do print, xind[i], yind[i]


sky_px = [-10:10]+sky_px
spaxel = make_array(n_elements(sky_px))
for i=0, n_elements(sky_px)-1 do begin
spaxel[i] = galaxy_data[xind[0],yind[0],sky_px[i]]
endfor


thing = [[sky_px*sxpar(header,'CD3_3')+sxpar(header,'CRVAL3')],[spaxel]]
print, Transpose(thing)
plot, sky_px*sxpar(header,'CD3_3')+sxpar(header,'CRVAL3'), spaxel

return
end
















