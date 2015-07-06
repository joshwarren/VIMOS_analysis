;; ==================================================================
;; 		Plot the outputted dynamics maps
;; ==================================================================
;; warrenj 20150330 Process to read the cube format and print it in
;; table form into a text file.

pro plot_results
 	galaxy = 'ngc3557'

	tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
		'/voronoi_2d_binning_output.txt'
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

;; Read tessellation file
	RDFLOAT, tessellation_File, x, y, bin_num, xBin, yBin, $
		SKIPLINE = 1, /SILENT 
        n_spaxels = (MAX(x) + 1) * (MAX(y) + 1)
        number_of_bins = MAX(bin_num)
	order = sort(bin_num)

;; Read results files - each entry in array corresponds to a bin (not
;; a spaxel)
	RDFLOAT, output_v, v_binned, /SILENT


;; ------------========== Total flux per bin ===========----------
;; ----------========= Reading the spectrum  =============---------

;; FILE_SEARCH returns an array even in cases where it only returns
;; one result. This is NOT equivalent to a scalar. 
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data, header

;	n_spaxels = n_elements(galaxy_data[*,0,0]) * $
;		n_elements(galaxy_data[0,*,0])



;; ----------========== Spatially Binning =============---------
;; endfor is near the end - after ppxf has been run on this bin.
for bin=0, n_bins-1 do begin
	spaxels_in_bin = WHERE(bin_num EQ bin, n_spaxels_in_bin)
        bin_flux = MAKE_ARRAY(number_of_bins, VALUE = 0d) 


for i = 0, n_spaxels_in_bin-1 do begin
	x_i = x[spaxels_in_bin[i]]
	y_i = y[spaxels_in_bin[i]]
for k = 0, n_elements(galaxy_data[x_i,y_i,*])-1 do begin
	bin_flux[i] = bin_flux[i] + galaxy_data[x_i, y_i, k]
endfor
endfor
;; bin_lin now contains linearly binned spectrum of the spatial bin.









        

;v_binned = v_binned - MEDIAN(v_binned)



b = UNIQ(bin_num, SORT(bin_num))

xNode = xBin[b]
yNode = yBin[b]






;sauron_colormap

plot_velfield, xNode, yNode, v_binned, FLUX = bin_flux





;;example use of contour
; Create a simple dataset:
;data = RANDOMU(seed, 9, 9)

; Plot the unsmoothed data:
;unsmooth = CONTOUR(data, TITLE='Unsmoothed', $
;   LAYOUT=[2,1,1], RGB_TABLE=13, /FILL, N_LEVELS=10)
; Draw the outline of the 10 levels
;outline1 = CONTOUR(data, N_LEVELS=10, /OVERPLOT)
 
; Plot the smoothed data:
;smooth = CONTOUR(MIN_CURVE_SURF(data), TITLE='Smoothed', $
;   /CURRENT, LAYOUT=[2,1,2], RGB_TABLE=13, $
;   /FILL, N_LEVELS=10)
; Draw the outline of the 10 levels
;outline2 = CONTOUR(MIN_CURVE_SURF(data), $
;   N_LEVELS=10, /OVERPLOT)


return
end
