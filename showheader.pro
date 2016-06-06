pro showheader

galaxy = "ngc1399"
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*_cube.fits') 
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		Galaxy + '-1/Q1/*[0-9].fits') 

;	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header


;print, header
print, ""
print, ""
	galaxy_noise_temp = MRDFITS(dataCubeDirectory[0], 2, header, /SILENT)
print, header

return
end
