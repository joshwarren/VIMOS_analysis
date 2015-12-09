pro showheader

galaxy = "ngc3557"
	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/reduced/' + $
		Galaxy + $
		'/cube/*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits') 
;	dataCubeDirectory = FILE_SEARCH('/Data/vimosindi/' + $
;		Galaxy + '-1/Q1/*[0-9].fits') 

	FITS_READ, dataCubeDirectory[0], galaxy_data_temp, header


print, header


return
end
