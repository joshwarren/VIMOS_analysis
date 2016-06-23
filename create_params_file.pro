pro create_params_file;, gal
;gal = 5
i_gal_beg = 0
code = "python"
;code = "IDL"

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
gals=[0,1,2,3,4,5,6,7,8]
gals=[6]
;gals = where(galaxies eq gal)


output_file = "~/VIMOS_project/analysis/params.txt"
CLOSE, 1
OPENW, 1, output_file
;i_gal=i_gal_beg
for i_gal=i_gal_beg, n_elements(gals)-1 do begin
;galaxy = galaxies[i_gal]
galaxy = galaxies[gals[i_gal]]


;tessellation_File = '/Data/vimosindi/analysis_sav_2016-02-09/' + galaxy + $
tessellation_File = '/Data/vimos/analysis/' + galaxy + $
	'/voronoi_2d_binning_output.txt'

RDFLOAT, tessellation_File, bin_num, COLUMNS = [3], $
	SKIPLINE = 1, /SILENT 
	
n_bins = max(bin_num)+1



for i = 0, n_bins do begin
if code eq "IDL" then PRINTF, 1, "idl -e 'errors, " + $
    strtrim(string(i_gal),2) + ", " + strtrim(string(i),2) + "'"
if code eq "python" then PRINTF, 1, "python errors2.py " + $
;     strtrim(string(i_gal),2) + " " + strtrim(string(i),2)
    strtrim(string(gals[i_gal]),2) + " " + strtrim(string(i),2)
endfor
endfor
CLOSE,1
print, "Done"
return
end
