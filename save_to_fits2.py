## Routine to transform the glamdring output to fits format in the 
## same setup as SARUON.

# import astropy.io.fits as pyfits 
# import numpy as np 
# import matplotlib.pyplot as plt 

# from cap_display_bins_generators import display_bins_generators 

# filename = "Cappellari2011a_Atlas3D_Paper1_UnBinned_Cubes_v1.0/"+/
# 	"MS_NGC3412_r3_C2D.fits" 
# hdu = pyfits.open(filename) 
# spectrum = hdu[0].data 
# table = hdu[2].data 

# x = table["A"] # Coordinates of the original spaxels in arcsec (nort is up) 
# y = table["D"] 
# flux = np.mean(spectrum, 1) # surface brightness 

# filename = "atlas3d_stellar_kinematics/cappellari2011a/"+/
# 	"PXF_bin_MS_NGC3412_r3_idl.fits.gz"
# hdu = pyfits.open(filename) 
# table = hdu[1].data 

# xgen = table["XS"] # Voronoi generators 
# ygen = table["YS"] 
# velbin = table["VPXF"] # Mean stellar velocity 

# display_bins_generators(xgen, ygen, velbin, x, y) 
# plt.tricontour(x, y, -2.5*np.log10(flux/np.max(flux)), 
# 	levels=np.arange(20)) # 1 mag contours



from checkcomp import checkcomp
cc = checkcomp()
import os
from astropy.io import fits 
import numpy as np 


def save(galaxy, instrument='vimos', debug=False, stellar=True, emission=True,
	absorption=True, absorption_nomask=True, kin_opt='kin', pop_opt='pop'):
	print galaxy, instrument
	if instrument=='vimos':
		cent_index = 4
	elif instrument == 'muse':
		cent_index = 1
	
	if cc.device == 'glamdring':
		if instrument == 'vimos':
			vin_dir = '%s/analysis' % (cc.home_dir)
		elif instrument == 'muse':
			vin_dir = '%s/analysis_muse' % (cc.home_dir)
	else:
		vin_dir = '%s/Data/%s/analysis' % (cc.base_dir, instrument)
	data_file =  "%s/galaxies.txt" % (vin_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_'+kin_opt)[0][0]
	col2 = np.where(file_headings=='SN_'+pop_opt)[0][0]
	x_cent_gals, y_cent_gals, SN_target_kin_gals, SN_target_pop_gals = \
		np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(cent_index,
		cent_index+1, col,col2))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	SN_target_kin=SN_target_kin_gals[i_gal]
	SN_target_pop=SN_target_pop_gals[i_gal]
	cent = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	vin_dir2 = str(vin_dir + '/%s/%s' % (galaxy, pop_opt)) 
	vin_dir += '/%s/%s' % (galaxy, kin_opt) 

	# if not debug:
	# 	if stellar:
	# 		pickle_file = '%s/pickled' % (vin_dir)
	# 		pickleFile = open("%s/dataObj.pkl" % (pickle_file), 'rb')
	# 		D = pickle.load(pickleFile)
	# 		pickleFile.close()
	# 	if emission or absorption or absorption_nomask:
	# 		pickle_file2 = '%s/pickled' % (vin_dir2)
	# 		pickleFile2 = open("%s/dataObj.pkl" % (pickle_file2), 'rb')
	# 		D2 = pickle.load(pickleFile2)
	# 		pickleFile2.close()
	# else:
	# 	from produce_maps import Ds
	# 	D = Ds()
	# 	D2 = Ds()

	# Stellar Kinematics
	
	if stellar:
		xBar, yBar = np.loadtxt("%s/%s/%s/setup/voronoi_2d_binning_output2.txt" %(
			vin_dir, galaxy, kin_opt), unpack=True, skiprows=1)
		bin_num = = np.loadtxt("%s/%s/%s/setup/voronoi_2d_binning_output.txt" %(
			vin_dir, galaxy, kin_opt), unpack=True, skiprows=1, dtype=int)


		hdr = fits.Header()
		hdr['SNR'] = SN_target_kin
		hdr['COMMENT'] = "Stellar kinematics of %s from %s " % (galaxy.upper(), 
			instrument.upper())+"data reduced/analysed by Warren,J. et.al.(2018)"
		primary_hdu = fits.PrimaryHDU(header=hdr)

		hdu = fits.BinTableHDU.from_columns([
			fits.Column(name='NO', format='I', array=np.arange(len(xBar))),
			fits.Column(name='XS', format='D', unit='arcsec', 
				array=-(xBar-cent[0])),
			fits.Column(name='YS', format='D', unit='arcsec', 
				array=-(yBar-cent[1])),
			fits.Column(name='BINSIZE', format='I', array=np.bincount(bin_num)),
			fits.Column(name='FLUX', format='D', unit='1E-15 erg s^-1 cm^-1', 
				array=D.flux), # Check flux units
			fits.Column(name='SNR', format='D', array=D.SNRatio),
			fits.Column(name='VEL', format='D', unit='km s^-1', 
				array=np.array(D.components['stellar'].plot['vel'])),
			fits.Column(name='E_VEL', format='D', unit='km s^-1', 
				array=D.components['stellar'].plot['vel'].uncert),
			fits.Column(name='SIGMA', format='D', unit='km s^-1', 
				array=np.array(D.components['stellar'].plot['sigma'])),
			fits.Column(name='E_SIGMA', format='D', unit='km s^-1', 
				array=D.components['stellar'].plot['sigma'].uncert)
			])


		f_new = fits.HDUList([primary_hdu, hdu])

		save_dir = '%s/Data/%s/analysed_fits/' % (cc.base_dir, instrument)
		if not os.path.exists(save_dir):
			os.makedirs(save_dir)
		if 'kin_opt' != 'kin':
			f_new.writeto('%s/%s_stellar_kine_%s.fits' %(save_dir, galaxy, kin_opt), 
				overwrite=True)
		else:
			f_new.writeto('%s/%s_stellar_kine.fits' %(save_dir, galaxy), 
				overwrite=True)

	if absorption or absorption_nomask or emission:
		xBar, yBar = np.loadtxt("%s/%s/%s/setup/voronoi_2d_binning_output2.txt" %(
			vin_dir2, galaxy, pop_opt), unpack=True, skiprows=1)
		bin_num = = np.loadtxt("%s/%s/%s/setup/voronoi_2d_binning_output.txt" %(
			vin_dir2, galaxy, pop_opt), unpack=True, skiprows=1, dtype=int)

	# Absorption line
	if absorption:
		if instrument=='vimos':
			plots = [
				"D2.absorption_line('G4300')",
				"D2.absorption_line('G4300',uncert=True)[1]",
				"D2.absorption_line('Fe4383')",
				"D2.absorption_line('Fe4383',uncert=True)[1]",
				"D2.absorption_line('Ca4455')",
				"D2.absorption_line('Ca4455',uncert=True)[1]",
				"D2.absorption_line('Fe4531')",
				"D2.absorption_line('Fe4531',uncert=True)[1]",
				"D2.absorption_line('H_beta')",
				"D2.absorption_line('H_beta',uncert=True)[1]",
				"D2.absorption_line('Fe5015')",
				"D2.absorption_line('Fe5015',uncert=True)[1]",
				"D2.absorption_line('Mg_b')",
				"D2.absorption_line('Mg_b',uncert=True)[1]"
				]
			str_plots = [
				"G4300",
				"e_G4300",
				"Fe4383",
				"e_Fe4383",
				"Ca4455",
				"e_Ca4455",
				"Fe4531",
				"e_Fe4531",
				'H_beta',
				'e_H_beta',
				"Fe5015",
				"e_Fe5015",
				'Mgb',
				'e_Mgb'
				]
			units = [
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom'
				]

		elif instrument == 'muse':
			plots = [
				"D2.absorption_line('H_beta')",
				"D2.absorption_line('H_beta',uncert=True)[1]",
				"D2.absorption_line('Fe5015')",
				"D2.absorption_line('Fe5015',uncert=True)[1]",
				"D2.absorption_line('Mg_b')",
				"D2.absorption_line('Mg_b',uncert=True)[1]",
				"D2.absorption_line('Fe5270')",
				"D2.absorption_line('Fe5270',uncert=True)[1]",
				"D2.absorption_line('Fe5335')",
				"D2.absorption_line('Fe5335',uncert=True)[1]",
				"D2.absorption_line('Fe5406')",
				"D2.absorption_line('Fe5406',uncert=True)[1]",
				"D2.absorption_line('Fe5709')",
				"D2.absorption_line('Fe5709',uncert=True)[1]",
				"D2.absorption_line('Fe5782')",
				"D2.absorption_line('Fe5782',uncert=True)[1]",
				"D2.absorption_line('NaD')",
				"D2.absorption_line('NaD',uncert=True)[1]",
				"D2.absorption_line('TiO1',remove_badpix=True)",
				"D2.absorption_line('TiO1',uncert=True,remove_badpix=True)[1]",
				"D2.absorption_line('TiO2',remove_badpix=True)",
				"D2.absorption_line('TiO2',uncert=True,remove_badpix=True)[1]"
				]
			str_plots = [
				'H_beta',
				'e_H_beta',
				"Fe5015",
				"e_Fe5015",
				'Mgb',
				'e_Mgb',
				'Fe5270',
				'e_Fe5270',
				'Fe5335',
				'e_Fe5335',
				'Fe5406',
				'e_Fe5406',
				'Fe5709',
				'e_Fe5709',
				'Fe5782',
				'e_Fe5782',
				'NaD',
				'e_NaD',
				'TiO1',
				'e_TiO1',
				'TiO2',
				'e_TiO2'
				]
			units = [
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'mag',
				'mag',
				'mag',
				'mag'
				]

		hdr = fits.Header()
		hdr['SNR'] = SN_target_pop
		hdr['COMMENT'] = "Absorption line strengths and best-fitting stellar "\
			+"populations of %s from %s " % (galaxy.upper(), instrument.upper()
			) + "data reduced/analysed by Warren,J. et.al.(2018)"
		primary_hdu = fits.PrimaryHDU(header=hdr)

		cols = [
			fits.Column(name='NO', format='I', array=np.arange(len(xBar))),
			fits.Column(name='XS', format='D', unit='arcsec', 
				array=-(xBar-cent[0])),
			fits.Column(name='YS', format='D', unit='arcsec', 
				array=-(yBar-cent[1])),
			fits.Column(name='BINSIZE', format='I', array=np.bincount(bin_num)),
			fits.Column(name='FLUX', format='D', unit='1E-15 erg s^-1 cm^-1', 
				array=D2.flux), # Check flux units
			fits.Column(name='SNR', format='D', array=D2.SNRatio)
			]
		cols.extend([
			fits.Column(name=str_plots[i], format='D', unit=units[i], 
				array=eval(p)) for i, p in enumerate(plots)])

		age = np.zeros(D2.number_of_bins)
		met = np.zeros(D2.number_of_bins)
		alp = np.zeros(D2.number_of_bins)
		unc_age = np.zeros(D2.number_of_bins)
		unc_met = np.zeros(D2.number_of_bins)
		unc_alp = np.zeros(D2.number_of_bins)

		if not debug:
			for j in xrange(D2.number_of_bins):
				ag, me, al = np.loadtxt('%s/pop/distribution/%i.dat' % (
					vin_dir2, j), unpack=True)

				for plot, unc_plot, pop in zip([age,met,alp],
					[unc_age,unc_met,unc_alp], [ag,me,al]):

					hist = np.histogram(pop, bins=40)
					x = (hist[1][0:-1]+hist[1][1:])/2
					hist = hist[0]
					plot[j] = x[np.argmax(hist)]

					gt_fwhm = hist >= np.max(hist)/2
					unc_plot[j] = np.max(x[gt_fwhm]) - np.min(x[gt_fwhm])

		str_plots = ['age', 'e_age', 'metalicity', 'e_metalicity', 
			'elpha', 'e_alpha']
		units = ['Gyr', 'Gyr', 'dex', 'dex', 'dex', 'dex']


		cols.extend([fits.Column(name=str_plots[i], format='D', unit=units[i], 
			array=eval(p)) for i, p in enumerate(['age', 'unc_age', 'met', 
			'unc_met', 'alp', 'unc_alp'])])

		hdu = fits.BinTableHDU.from_columns(cols)
		f_new = fits.HDUList([primary_hdu, hdu])

		save_dir = '%s/Data/%s/analysed_fits/' % (cc.base_dir, instrument)
		if not os.path.exists(save_dir):
			os.makedirs(save_dir)
		if pop_opt != 'pop':
			f_new.writeto('%s/%s_absorption_line_%s.fits' %(save_dir, 
				galaxy, pop_opt), overwrite=True)
		else:
			f_new.writeto('%s/%s_absorption_line.fits' %(save_dir, galaxy), 
				overwrite=True)


	# Absorption line - no emission line masks
	if absorption_nomask:
		if instrument=='vimos':
			plots = [
				"D2.absorption_line('G4300')",
				"D2.absorption_line('G4300',uncert=True)[1]",
				"D2.absorption_line('Fe4383')",
				"D2.absorption_line('Fe4383',uncert=True)[1]",
				"D2.absorption_line('Ca4455')",
				"D2.absorption_line('Ca4455',uncert=True)[1]",
				"D2.absorption_line('Fe4531')",
				"D2.absorption_line('Fe4531',uncert=True)[1]",
				"D2.absorption_line('H_beta')",
				"D2.absorption_line('H_beta',uncert=True)[1]",
				"D2.absorption_line('Fe5015')",
				"D2.absorption_line('Fe5015',uncert=True)[1]",
				"D2.absorption_line('Mg_b')",
				"D2.absorption_line('Mg_b',uncert=True)[1]"
				]
			str_plots = [
				"G4300",
				"e_G4300",
				"Fe4383",
				"e_Fe4383",
				"Ca4455",
				"e_Ca4455",
				"Fe4531",
				"e_Fe4531",
				'H_beta',
				'e_H_beta',
				"Fe5015",
				"e_Fe5015",
				'Mgb',
				'e_Mgb'
				]
			units = [
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom'
				]

		elif instrument == 'muse':
			plots = [
				"D2.absorption_line('H_beta')",
				"D2.absorption_line('H_beta',uncert=True)[1]",
				"D2.absorption_line('Fe5015')",
				"D2.absorption_line('Fe5015',uncert=True)[1]",
				"D2.absorption_line('Mg_b')",
				"D2.absorption_line('Mg_b',uncert=True)[1]",
				"D2.absorption_line('Fe5270')",
				"D2.absorption_line('Fe5270',uncert=True)[1]",
				"D2.absorption_line('Fe5335')",
				"D2.absorption_line('Fe5335',uncert=True)[1]",
				"D2.absorption_line('Fe5406')",
				"D2.absorption_line('Fe5406',uncert=True)[1]",
				"D2.absorption_line('Fe5709')",
				"D2.absorption_line('Fe5709',uncert=True)[1]",
				"D2.absorption_line('Fe5782')",
				"D2.absorption_line('Fe5782',uncert=True)[1]",
				"D2.absorption_line('NaD')",
				"D2.absorption_line('NaD',uncert=True)[1]",
				"D2.absorption_line('TiO1',remove_badpix=True)",
				"D2.absorption_line('TiO1',uncert=True,remove_badpix=True)[1]",
				"D2.absorption_line('TiO2',remove_badpix=True)",
				"D2.absorption_line('TiO2',uncert=True,remove_badpix=True)[1]"
				]
			str_plots = [
				'H_beta',
				'e_H_beta',
				"Fe5015",
				"e_Fe5015",
				'Mgb',
				'e_Mgb',
				'Fe5270',
				'e_Fe5270',
				'Fe5335',
				'e_Fe5335',
				'Fe5406',
				'e_Fe5406',
				'Fe5709',
				'e_Fe5709',
				'Fe5782',
				'e_Fe5782',
				'NaD',
				'e_NaD',
				'TiO1',
				'e_TiO1',
				'TiO2',
				'e_TiO2'
				]
			units = [
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'Angstrom',
				'mag',
				'mag',
				'mag',
				'mag'
				]

		hdr = fits.Header()
		hdr['SNR'] = SN_target_pop
		hdr['COMMENT'] = "Absorption line strengths with no detection limit on "\
			+"subtracted emission lines of %s from %s " % (galaxy.upper(), 
			instrument.upper()) + "data reduced/analysed by Warren,J. et.al.(2018)"
		primary_hdu = fits.PrimaryHDU(header=hdr)

		cols = [
			fits.Column(name='NO', format='I', array=np.arange(len(xBar))),
			fits.Column(name='XS', format='D', unit='arcsec', 
				array=-(xBar-cent[0])),
			fits.Column(name='YS', format='D', unit='arcsec', 
				array=-(yBar-cent[1])),
			fits.Column(name='BINSIZE', format='I', array=np.bincount(bin_num)),
			fits.Column(name='FLUX', format='D', unit='1E-15 erg s^-1 cm^-1', 
				array=D2.flux), # Check flux units
			fits.Column(name='SNR', format='D', array=D2.SNRatio)
			]
		for i, p in enumerate(plots):
			print i, p, str_plots[i], units[i]
		cols.extend([
			fits.Column(name=str_plots[i], format='D', unit=units[i], 
				array=eval(p.replace(')',',nomask=True)'))) for i, p in 
				enumerate(plots)])

		# age = np.zeros(D2.number_of_bins)
		# met = np.zeros(D2.number_of_bins)
		# alp = np.zeros(D2.number_of_bins)
		# unc_age = np.zeros(D2.number_of_bins)
		# unc_met = np.zeros(D2.number_of_bins)
		# unc_alp = np.zeros(D2.number_of_bins)

		# if not debug:
		# 	for j in xrange(D2.number_of_bins):
		# 		ag, me, al = np.loadtxt('%s/pop/distribution/%i.dat' % (
		# 			vin_dir2, j), unpack=True)

		# 		for plot, unc_plot, pop in zip([age,met,alp],
		# 			[unc_age,unc_met,unc_alp], [ag,me,al]):

		# 			hist = np.histogram(pop, bins=40)
		# 			x = (hist[1][0:-1]+hist[1][1:])/2
		# 			hist = hist[0]
		# 			plot[j] = x[np.argmax(hist)]

		# 			gt_fwhm = hist >= np.max(hist)/2
		# 			unc_plot[j] = np.max(x[gt_fwhm]) - np.min(x[gt_fwhm])

		# str_plots = ['age', 'e_age', 'metalicity', 'e_metalicity', 
		# 	'elpha', 'e_alpha']
		# units = ['Gyr', 'Gyr', 'dex', 'dex', 'dex', 'dex']


		# cols.extend([fits.Column(name=str_plots[i], format='D', unit=units[i], 
		# 	array=eval(p)) for i, p in enumerate(['age', 'unc_age', 'met', 
		# 	'unc_met', 'alp', 'unc_alp'])])

		hdu = fits.BinTableHDU.from_columns(cols)
		f_new = fits.HDUList([primary_hdu, hdu])

		save_dir = '%s/Data/%s/analysed_fits/' % (cc.base_dir, instrument)
		if not os.path.exists(save_dir):
			os.makedirs(save_dir)
		if pop_opt != 'pop':
			f_new.writeto('%s/%s_absorption_line_nomask_%s.fits' %(save_dir, 
				galaxy, pop_opt), overwrite=True)
		else:
			f_new.writeto('%s/%s_absorption_line_nomask.fits' %(save_dir, galaxy), 
				overwrite=True)


	# Emission lines
	if emission:
		if '[OIII]5007d' in D2.e_components:
			hdr = fits.Header()
			hdr['SNR'] = SN_target_pop
			hdr['COMMENT'] = "Emission line fluxes and kinematics "\
				+"%s from %s " % (galaxy.upper(), instrument.upper()
				) + "data reduced/analysed by Warren,J. et.al.(2018)"
			primary_hdu = fits.PrimaryHDU(header=hdr)

			cols = [
				fits.Column(name='NO', format='I', array=np.arange(len(xBar))),
				fits.Column(name='XS', format='D', unit='arcsec', 
					array=-(xBar-cent[0])),
				fits.Column(name='YS', format='D', unit='arcsec', 
					array=-(yBar-cent[1])),
				fits.Column(name='BINSIZE', format='I', array=np.bincount(bin_num)),
				fits.Column(name='FLUX', format='D', unit='1E-15 erg s^-1 cm^-1', 
					array=D2.flux), # Check flux units
				fits.Column(name='SNR', format='D', array=D2.SNRatio)
				]
			cols.extend([
				fits.Column(name='VEL', format='D', unit='km s^-1', 
					array=np.array(D2.components['[OIII]5007d'].plot['vel'])),
				fits.Column(name='E_VEL', format='D', unit='km s^-1', 
					array=D2.components['[OIII]5007d'].plot['vel'].uncert),
				fits.Column(name='SIGMA', format='D', unit='km s^-1', 
					array=np.array(D2.components['[OIII]5007d'].plot['sigma'])),
				fits.Column(name='E_SIGMA', format='D', unit='km s^-1', 
					array=D2.components['[OIII]5007d'].plot['sigma'].uncert)
				])

			cols.extend(np.array([[
				fits.Column(name=p.replace('[','').replace(']',''), format='D', 
				unit='1E-15 erg s^-1 cm^-1', array=D2.components[p].flux),
				fits.Column(name='e_'+p.replace('[','').replace(']',''), format='D', 
				unit='1E-15 erg s^-1 cm^-1', array=D2.components[p].flux.uncert),
				fits.Column(name='anr_'+p.replace('[','').replace(']',''), 
					format='D', array=D2.components[p].amp_noise)] 
				for i, p in enumerate(D2.e_components)]).flatten())

			hdu = fits.BinTableHDU.from_columns(cols)
			f_new = fits.HDUList([primary_hdu, hdu])

			save_dir = '%s/Data/%s/analysed_fits/' % (cc.base_dir, instrument)
			if not os.path.exists(save_dir):
				os.makedirs(save_dir)
			if pop_opt != 'pop':
				f_new.writeto('%s/%s_emission_line_%S.fits' %(save_dir, galaxy,
					pop_opt), overwrite=True)
			else:
				f_new.writeto('%s/%s_emission_line.fits' %(save_dir, galaxy), 
					overwrite=True)

# Routine to strink the size of the fits file to something more manageable on my 
# laptop
def save_muse(galaxy):
	if cc.device != 'uni':
		raise ValueError('This routine is to be run on your university computer.'
			+' Chump!')

	from errors2_muse import get_dataCubeDirectory
	f = fits.open(get_dataCubeDirectory(galaxy))

	ex0 = f[0]
	ex1 = fits.ImageHDU(np.array([np.nansum(f[1].data, axis=0)]), f[1].header, 
		name='DATA')
	ex2 = fits.ImageHDU(np.array([np.sqrt(np.nansum(f[2].data**2, axis=0))]), 
		f[2].header, name='ERROR')
	ex3 = fits.ImageHDU(np.nansum(f[3].data).astype(bool).astype(int), f[3].header, 
		name='BADPIX')
	f_new = fits.HDUList([ex0, ex1, ex2, ex3])

	corr_fits_file = '%s/Data/muse/%s/%s.clipped_home.fits' % (cc.base_dir, galaxy, 
		galaxy)

	f_new.writeto(corr_fits_file, overwrite=os.path.isfile(corr_fits_file))

if __name__=='__main__':
	if 'home' in cc.device:
		for galaxy in ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 
			'ngc0612', 'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']:
			save(galaxy, debug=False, stellar=False, absorption=True, 
				absorption_nomask=False, emission=False)
	elif cc.device == 'uni':
		for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
			# save_muse(galaxy)
			save(galaxy, instrument='muse', debug=False, stellar=True, 
				absorption=True, absorption_nomask=True, emission=True)