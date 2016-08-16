;; ==================================================================
;; Rebinning using Voronoi tessellation
;; ==================================================================
;; warrenj 20150515 Process to rebin spaxels together in order to
;; create a minimum S/N ratio. targetSN is approx 10-15 for just v and
;; sigma, and 40-50 for h3 and h4. 
 


;; ----------===============================================---------
;; ----------======== Check overwrite of target SN =========---------
;; ----------===============================================---------
pro check_overwrite, new, old
if new ne old then begin
A=''
read, A, prompt='Are you sure you want to overwrite the old target of ' + $
    strtrim(string(old),2) + ' with a new target of ' + $
    strtrim(string(new),2) + '? (Y/N) '
if (A eq "N") or (a eq "n") then new = old
endif
return
end







pro binning_spaxels, galaxy, discard=discard, targetSN=targetSN

;; ----------===============================================---------
;; ----------============ Default parameters ===============---------
;; ----------===============================================---------

data_file = "/Data/vimos/analysis/galaxies.txt"
readcol, data_file, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, $
    y_gals, SN_used_gals, skipline=1, format='A,D,D,D,I,I,I', /SILENT
i_gal = where(galaxy_gals eq galaxy)

if i_gal ne -1 then begin
if not keyword_set(targetSN) then targetSN=SN_used_gals[i_gal]
if keyword_set(targetSN) then check_overwrite, targetSN, SN_used_gals[i_gal]
endif else if not keyword_set(targetSN) then targetSN = 30.0

if not keyword_set(discard) then discard=2

;	galaxy = 'ngc3557'
;	galaxy = 'ic1459'
;	discard = 2
;	targetSN = 30.0

;; ----------================= Save SN_used ===============---------
if i_gal eq -1 then begin
    galaxy_gals = [galaxy_gals, galaxy]
    SN_used_gals = [SN_used_gals, targetSN]
endif else begin
    SN_used_gals[i_gal] = targetSN
endelse


forprint2, galaxy_gals, z_gals, vel_gals, sig_gals, x_gals, y_gals, $
    SN_used_gals, textout=data_file, /SILENT, WIDTH=90, $
    Comment = "Galaxy      z     velocity     velocity dispersion    x     y     Target SN"


;; ----------================ Find S/N ================------------
; Final wildcard notes that depending on the method used the quadrants
; may or may not have been flux calibrated. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimos/cubes/' + $
		galaxy + '.cube.combined.fits')

;;Bins just one quadrant - must have used rss2cube_quadrant.pro script
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		galaxy + $
;		'-3/Q2/calibrated/cube/*_fluxcal_cube.fits')

	;galaxy_data = MRDFITS(dataCubeDirectory[0], 1, header, /SILENT)
;	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 1, /SILENT)
	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header, EXTEN_NO=0
	FITS_READ, dataCubeDirectory[0], galaxy_noise_temp, EXTEN_NO=1
	FITS_READ, dataCubeDirectory[0], galaxy_badpix_temp, EXTEN_NO=3

	s = size(galaxy_data_temp)

	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_noise = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_badpix = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])

	galaxy_data = galaxy_data_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]
	galaxy_noise = galaxy_noise_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]
	galaxy_badpix = galaxy_badpix_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]
;	galaxy_noise[where(galaxy_noise<0)]=0
;print, n_elements(where(galaxy_noise<0))/1600
	
	s = size(galaxy_data)
	n_spaxels = s[1]*s[2]
	CRVAL_spec = sxpar(header,'CRVAL3')
	CDELT_spec = sxpar(header,'CDELT3')


	signal = MAKE_ARRAY(n_spaxels)
	noise = MAKE_ARRAY(n_spaxels)
	x = MAKE_ARRAY(n_spaxels)
	y = MAKE_ARRAY(n_spaxels)




;; collapsing the spectrum for each spaxel.
for i = 0, s[1]-1 do begin
for j = 0, s[2]-1 do begin

	gap=12
	ignore = FIX((5581 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap  
	ignore2 =FIX((5199 - CRVAL_spec)/CDELT_spec) + [-1,+1]*gap 
;; h is the spectrum with the peak enclosed by 'ignore' removed.
        if 5581 lt CRVAL_spec+s[3]*CDELT_spec then begin
		h =[reform(galaxy_data[i,j,0:ignore[0]]),reform(galaxy_data[i,j,ignore[1]:*])]
        endif else h = reform(galaxy_data)

	h =[h[0:ignore2[0]],h[ignore2[1]:*]]


	half = s[3]/2
	a = h/MEDIAN(h) - h[4:*]/MEDIAN(h)
	a[WHERE(~FINITE(a))] = 0
	
	lower_limit = MIN(WHERE(ABS(a) GT 0.2), MAX=upper_limit)
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

upper_limit = (5300-CRVAL_spec)/CDELT_spec
;print, i,j,i*s[1] + j, n_spaxels
;	signal[i*s[1] + j] = STDDEV(galaxy_data[i, j,
;	lower_limit:upper_limit])
	signal[i*s[1] + j] = MEAN(galaxy_data[i, j, lower_limit:upper_limit])
	noise[i*s[1] + j] = MEAN(galaxy_noise[i, j, lower_limit:upper_limit])

;; Assign x and y
	x(i*s[1]+j) = i
	y(i*s[2]+j) = j

	a =where(galaxy_badpix[i,j,*] eq 1, count)
	if count ne 0 then begin
		signal[i*s[1] + j] = 0
		noise[i*s[1] + j] = 0
	endif

;if i*s[1] + j eq 1287 then forprint, galaxy_data[i,j,*], galaxy_noise[i,j,*], textout=2
endfor
endfor

n_spaxels = n_elements(signal)

;data_file = '~/VIMOS_project/analysis_v2/rebinning/voronoi_2d_binning_input.txt'
;rdfloat, data_file, x, y, signal, noise, SKIPLINE=3, NUMLINE=3107, /DOUBLE

; Load a colortable and open a graphic window
;
loadct, 5, /SILENT
r = GET_SCREEN_SIZE()
window, xsize=r[0]*0.4, ysize=r[1]*0.8

; Perform the actual computation. The vectors
; (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
; are all generated in *output*
med = median(signal)

;; ****** Fudging!!! ******
;l=where(noise/med gt 100000)
;print, l
;signal[l]=make_array(n_elements(l), value=0)
;noise[l]=make_array(n_elements(l), value=0)
signal[where(~finite(noise))] = 0.000
noise[where(~finite(noise))] = 0.000
;print, noise[where(noise gt signal)]

;print, total(signal)/sqrt(total(noise^2))


signal[where(signal le 0)]=0.00001
voronoi_2d_binning, x, y, signal, noise, targetSN, $
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /PLOT, /QUIET




order = sort(binNum)
xBin = MAKE_ARRAY(n_elements(x))
yBin = MAKE_ARRAY(n_elements(y))
;; spaxel number
i = 0
for bin = 0, MAX(binNum) do begin
while(i LT n_spaxels && bin EQ binNum[order[i]]) do begin
	xBin[order[i]] = xBar[bin]
	yBin[order[i]] = yBar[bin]
;; move onto next spaxel
i = i + 1
endwhile
endfor	








; Save to a text file the initial coordinates of each pixel together
; with the corresponding bin number computed by this procedure.
; binNum uniquely specifies the bins and for this reason it is the only
; number required for any subsequent calculation on the bins.
;
FILE_MKDIR, '/Data/vimos/analysis/' + galaxy
OPENW, 1, '/Data/vimos/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt'
CLOSE,1
astrolib

forprint, x, y, binNum, xBin, yBin, $
	TEXTOUT='/Data/vimos/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt', $
	COMMENT='          X"              Y"           BIN_NUM           XBIN           YBIN', /SILENT



forprint, xBar, yBar, $
	TEXTOUT='/Data/vimos/analysis/' + galaxy + $
	'/voronoi_2d_binning_output2.txt', $
	COMMENT='XBAR           YBAR', /SILENT



END
;----------------------------------------------------------------------------





