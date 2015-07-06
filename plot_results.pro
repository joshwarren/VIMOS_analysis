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
	order = sort(bin_num)

;; Read results files - each entry in array corresponds to a bin (not
;; a spaxel)
	RDFLOAT, output_v, v_binned, /SILENT



;v_binned = v_binned - MEDIAN(v_binned)

b = UNIQ(bin_num, SORT(bin_num))

xNode = xBin[b]
yNode = yBin[b]

;sauron_colormap

;n_levels = n_elements(v_binned)-1
;colors = MAKE_ARRAY(n_levels)
;for i=0, n_levels-1 do colors(i)=10000*i

;loadct, 8
;tvscl, DIST(300)
plot_velfield, xNode, yNode, v_binned;, NCOLORS = n_levels;, C_COLORS =colors 






;  data = cgDemoData(2)
;    x = cgScaleVector(Findgen(41), 0, 100)
;    y = cgScaleVector(Findgen(41), 0, 50)
;
;print,size(data)
;print,size(x)
;print,size(y)
;
;   levels = 12
;    cgLoadCT, 33, NColors=12, Bottom=3
;
    ; Draw the first plot. Let IDL select contour intervals by
    ; using the NLEVELS keyword.
;    Window, 0, Title='IDL Selected Contour Intervals', XSize=300, YSize=400
;    SetDecomposedState, 0, CurrentState=state
;    Contour, v_binned, xNode, yNode, /Fill, C_Colors=Indgen(levels)+3, Background=cgColor('white'), $
;       NLevels=levels, Position=[0.1, 0.1, 0.9, 0.80], Color=cgColor('black')
;    Contour, v_binned, xNode, yNode, /Overplot, NLevels=levels, /Follow, Color=cgColor('black')
;    SetDecomposedState, state

return
end
