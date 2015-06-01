;; ==================================================================
;; Finding signal and noise for each spaxel
;; ==================================================================
;; warrenj 20150515 process to create a txt file for the routine
;; binning_spaxels to then rebin to create a desired S/N ratio by
;; binning some spaxels together. 



pro find_SN, galaxy, signal, noise
	
;	galaxy = 'ngc3557'
	OB = '1'


	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits')

	FITS_READ, dataCubeDirectory[0], galaxy_data, header
;print, size(galaxy_data)

;plot, galaxy_data(10,10,*)



	spaxel_spectrum = MAKE_ARRAY(2800)
for i=0, 39 do begin
for j=0, 39 do begin
;print, 'i = ' + string(i) + '  j = ' + string(j)
;print, median(galaxy_data(i,j,*))

for k = 0, 2799 do begin
	spaxel_spectrum(k) = galaxy_data(i,j,k)
endfor
	boxplot = CREATEBOXPLOTDATA(spaxel_spectrum)

	signal(i*40+j) = MEDIAN(spaxel_spectrum)
;; Estimate the noise to be half of the difference between the upper
;; and lower quartiles. Spectrum is relatively flat so may be fair. 
	noise(i*40+j) = (boxplot(3)-boxplot(1))/2

endfor
endfor


return
end





;; ==================================================================
;; Rebinning using Voronoi tessellation
;; ==================================================================
;; warrenj 20150515 Process to rebin spaxels together in order to
;; create a minimum S/N ratio. targetSN is approx 10-15 for just v and
;; sigma, and 40-50 for h3 and h4. 


pro binning_spaxels;, galaxy
;
; Usage example for the procedure VORONOI_2D_BINNING.
; Type the name VORONOI_2D_BINNING_EXAMPLE at the IDL
; prompt to run this example.
;
; Here columns 1-4 of the text file contain respectively the x, y
; coordinates of each SAURON lens and the corresponding Signal and
; Noise. 
;
	galaxy = 'ngc3557'


	signal = MAKE_ARRAY(40*40)
	noise = MAKE_ARRAY(40*40)
	find_SN, galaxy, signal, noise

;data_file = '~/VIMOS_project/analysis_v2/rebinning/voronoi_2d_binning_input.txt'

;rdfloat, data_file, x, y, signal, noise, SKIPLINE=3, NUMLINE=3107, /DOUBLE
	targetSN = 10.0

; Load a colortable and open a graphic window
;
loadct, 5, /SILENT
r = GET_SCREEN_SIZE()
window, xsize=r[0]*0.4, ysize=r[1]*0.8

; Perform the actual computation. The vectors
; (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
; are all generated in *output*
;

;; Asigning x and y
	x = MAKE_ARRAY(40*40)
	y = MAKE_ARRAY(40*40)
for i = 0, 39 do begin
for j = 0 , 39 do begin
x(i*40+j) = i
y(i*40+j) = j
endfor
endfor

voronoi_2d_binning, x, y, signal, noise, targetSN, $
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /PLOT, /QUIET


; Save to a text file the initial coordinates of each pixel together
; with the corresponding bin number computed by this procedure.
; binNum uniquely specifies the bins and for this reason it is the only
; number required for any subsequent calculation on the bins.
;
FILE_MKDIR, '/Data/vimosindi/analysis/' + galaxy
astrolib
forprint, x, y, binNum, TEXTOUT='/Data/vimosindi/analysis/' + galaxy + $
'/voronoi_2d_binning_output.txt', $
    COMMENT='          X"              Y"           BIN_NUM'

END
;----------------------------------------------------------------------------





