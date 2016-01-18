pro create_params_file
i_gal_beg = 4

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']

output_file = "~/VIMOS_project/analysis/params.txt"
CLOSE, 1
OPENW, 1, output_file

for i_gal=i_gal_beg, 9 do begin
galaxy = galaxies[i_gal]

tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt'

RDFLOAT, tessellation_File, bin_num, COLUMNS = [3], $
	SKIPLINE = 1, /SILENT 
	
n_bins = max(bin_num)+1


;PRINTF, 1, "idl -e '.compile errors_gas'"

for i = 0, n_bins-1 do begin
PRINTF, 1, "idl -e 'errors, " + strtrim(string(i_gal),2) + ", " + strtrim(string(i),2) + "'"
endfor
endfor
CLOSE,1
print, "Done"
return
end
