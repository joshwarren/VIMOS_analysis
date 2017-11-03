pro showheader

galaxy = "ngc1399"
	dataCubeDirectory = FILE_SEARCH('/Data/vimos/cubes/' + $
		Galaxy + '.cube.combined.fits') 
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		Galaxy + '-1/Q1/*[0-9].fits') 

;	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header


;print, header
print, ""
print, ""
	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 0, header, /SILENT)
print, header

return
end
