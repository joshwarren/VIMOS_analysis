;; ==================================================================
;; Rebinning using Voronoi tessellation
;; ==================================================================
;; warrenj 20150515 Process to rebin spaxels together in order to
;; create a minimum S/N ratio. targetSN is approx 10-15 for just v and
;; sigma, and 40-50 for h3 and h4. 
 

pro binning_spaxels;, galaxy = "ngc3557", dicard = 2, targetSN = 30.0

;	galaxy = 'ngc3557'
	galaxy = 'ic1459'
	discard = 2
	targetSN = 30.0


;; ----------================ Find S/N ================------------
; Final wildcard notes that depending on the method used the quadrants
; may or may not have been flux calibrated. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		galaxy + $
		'/cube/*crcl_oextr1*vmcmb_darc_cexp_cube.fits')

;;Bins just one quadrant - must have used rss2cube_quadrant.pro script
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		galaxy + $
;		'-3/Q2/calibrated/cube/*_fluxcal_cube.fits')


	;galaxy_data = MRDFITS(dataCubeDirectory[0], 1, header, /SILENT)
	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 2, /SILENT)
	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header
;print, size(galaxy_data)

;plot, galaxy_data(10,10,*)



	s = size(galaxy_data_temp)
	galaxy_data = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])
	galaxy_noise = MAKE_ARRAY(s[1]-2*discard,s[2]-2*discard,s[3])

	galaxy_data = galaxy_data_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]
	galaxy_noise = galaxy_noise_temp[[discard:s[1]-discard-1], $
		[discard:s[2]-discard-1],*]
;	galaxy_noise[where(galaxy_noise<0)]=0
;print, n_elements(where(galaxy_noise<0))/1600

	
	s = size(galaxy_data)
	n_spaxels = s[1]*s[2]


	signal = MAKE_ARRAY(n_spaxels)
	noise = MAKE_ARRAY(n_spaxels)
	x = MAKE_ARRAY(n_spaxels)
	y = MAKE_ARRAY(n_spaxels)



;; collapsing the spectrum for each spaxel.
for i = 0, s[1]-1 do begin
for j = 0, s[2]-1 do begin

	signal[i*s[1] + j] = MEAN(galaxy_data[i, j, *])
	noise[i*s[1] + j] = MEAN(galaxy_noise[i, j, *])
	


;; Assign x and y
	x(i*s[1]+j) = i
	y(i*s[2]+j) = j
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
FILE_MKDIR, '/Data/vimosindi/analysis/' + galaxy
astrolib
forprint, x, y, binNum, xBin, yBin, $
	TEXTOUT='/Data/vimosindi/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt', $
	COMMENT='          X"              Y"           BIN_NUM           XBIN           YBIN'



forprint, xBar, yBar, $
	TEXTOUT='/Data/vimosindi/analysis/' + galaxy + $
	'/voronoi_2d_binning_output2.txt', $
	COMMENT='XBAR           YBAR'

END
;----------------------------------------------------------------------------





