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
	RDFLOAT, tessellation_File, x, y, bin_num, SKIPLINE = 1, $
		/SILENT 
	n_spaxels = (MAX(x) + 1) * (MAX(y) + 1)
	order = sort(bin_num)

;; Read results files - each entry in array corresponds to a bin (not
;; a spaxel)
	RDFLOAT, output_v, v_binned, /SILENT

;; 2D array to hold the results in spaxel form.
	v_map = MAKE_ARRAY(MAX(x) + 1, MAX(y) + 1)

i = 0
;; defines the working bin
;for bin = 0, MAX(bin_num) -1 do begin
for bin = 0, MAX(bin_num) do begin
;; loops over all spaxels within that bin
while (i LT n_spaxels && bin EQ bin_num[order[i]]) do begin
	
	v_map[x[order[i]], y[order[i]]] = v_binned[bin]
;; moves onto the next spaxel
i = i + 1
endwhile
endfor

CONTOUR, v_map, /FILL

return
end
