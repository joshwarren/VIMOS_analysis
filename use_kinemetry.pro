;################################################################
; This is addapted from kinemetry_example.pro for my own use with
; VIMOS observations.
;###############################################################

PRO do_work, gal, type
	plot = 0
	; Odd or even moment?
	if type eq 'vel' then even=0 else even=1


	print, 'Have you got the most recent files for '+gal+'?'
	; Select correct file
	if type eq 'flux' then begin
	file = '/Data/vimos/analysis/'+gal+'/results/4200-/flux.dat'
	endif else begin
	file = '/Data/vimos/analysis/'+gal+'/results/4200-/gal_stellar_'+type+'.dat'
	endelse
	; read in field
	rdfloat, file, velbin, er_velbin
	if type eq 'flux' then er_velbin = velbin*0+1 else er_velbin=er_velbin*0.01

	; Read in binning
	file = '/Data/vimos/analysis/'+gal+'/voronoi_2d_binning_output.txt'
	rdfloat, file, _,_,bin_num, xbin,ybin, skipline=1

	file = '/Data/vimos/analysis/galaxies.txt'
	readcol, file, galaxy_gals,_,_,_,x0,y0,_, skipline=1,format='A,D,D,D,I,I,I'
	i_gal = where(galaxy_gals eq gal)

	b = uniq(bin_num,sort(bin_num))
	xbin = xbin[b]
	ybin = ybin[b]

	; Center the origin on the center of the galaxy
	x_cent = max(xbin)/2
	y_cent = max(ybin)/2
	xbin = xbin - x_cent
	ybin = ybin - y_cent

	; kinemetry on maps
	t=systime(1)
	KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, x0=x0[i_gal]-x_cent, $
		y0=y0[i_gal]-y_cent, ntrm=6, scale=0.33, name=gal,er_cf=er_cf, $
		er_pa=er_pa, even=even, ERROR=er_velbin, er_q=er_q, /verbose
	;catch, caught_error
	;if caught_error ne 0 then catch, /cancel
	print, systime(1) -t, 'seconds'

	; kinemetry parameters as defined in Krajnovic et al. (2006)
	k0 = cf[*,0]
	k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
	k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
	k51 = k5/k1
	erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
	erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
	erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1 


	file = '/Data/vimos/analysis/'+gal+'/kinemetry_'+type+'.txt'
	forprint2, rad, pa, er_pa, q, er_q, k1, erk1, k51, erk51, width=200, TEXTOUT = file, $
		/SILENT, comment='  radius(pix)      pa(deg)        err         ellip        err           k1           err          k51         err'


	if keyword_set(plot) then begin
		; plot coeffs.
		r = GET_SCREEN_SIZE()
		window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
		!p.charsize=3
		!y.style=1
		!p.multi=[0,1,4]
		!Y.MARGIN=[0,0] ; show plots with shared X axis
		!Y.OMARGIN=[5,3] ; allow for space for the axis labels
		ploterror, rad, pa, er_pa, PSYM=-5, TITLE=gal, xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[min(pa),max(pa)]
		ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q'
		ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', YTITLE='k1 [km/s]',YRANGE=[0,245]
		ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [arcsec]', YTITLE='k5/k1', YRANGE=[0,0.13]
		!P.MULTI=0
		!Y.MARGIN=[4,2] ; back to default values
		!Y.OMARGIN=[0,0]
		!p.charsize=1
		WAIT, 3000000000
	endif

END


pro use_kinemetry
;	gal = 'eso443-g024'
	gals=['ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34', 'eso443-g024']
	for i=0,9 do begin
		gal=gals[i]
		print, gal

		do_work, gal, 'flux'
		do_work, gal, 'sigma'
		do_work, gal, 'vel' ; Not working for eso443-g024

	endfor

END