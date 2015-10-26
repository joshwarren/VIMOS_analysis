pro create_params_file

galaxy = 'ngc3557'

tessellation_File = '/Data/vimosindi/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt'

RDFLOAT, tessellation_File, x, y, bin_num, COLUMNS = [1,2,3], $
	SKIPLINE = 1, /SILENT 
	
n_bins = max(bin_num)

output_file = "~/VIMOS_project/analysis/params.txt"
CLOSE, 1
OPENW, 1, output_file
PRINTF, 1, "idl -e '.compile errors'"
for i = 0, n_bins-1 do begin
PRINTF, 1, "idl -e 'errors, " + strtrim(string(i),2) + "'"
endfor
CLOSE,1
return
end
