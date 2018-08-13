# warrenj 20180203 Object to use in place of Data object which will load data 
# from fits files when 
# needed.

import numpy as np
# import ppxf_util as util
# from absorption import absorption
# from glob import glob
from checkcomp import checkcomp
cc = checkcomp()
from tools import moving_weighted_average
# from scipy.ndimage.filters import gaussian_filter1d
# from scipy.interpolate import interp1d
# from itertools import groupby, count
from astropy.io import fits
import warnings, os



c = 299792.458 # speed of light in km/s

def remove_brackets(s):
	return s.replace('[','').replace(']','')

def add_brackets(s):
	out = str(s)
	for l in ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731', 
		'[OIII]5007d', '[NI]d', '[OI]6300d', '[NII]6583d']:
		l_no_bracket = remove_brackets(l)
		if s == l_no_bracket:
			out = l
	return out


class Data(object):
# Attributes:
# number_of_bins: int
# bin: array of bin objects (see below) - holds most of the data - is a 
#	'setter' class
# e_line: dictionary of emission_data objects (see below) - is a 'getter' class
# e_components: string list of emission lines
# unbinned_flux: float array of unbinned flux in 2D array
# flux: float array of total integrated flux (observed - not fitted)
# xBar: float array of central point's x-coord for each bin
# yBar: as xBar
# n_spaxels_in_bin: int array of number of spaxels in each bin
# norm_method: str, normalisation methods for velocity fields:
#		lwv: 	(default) luminosity weighted mean of the whole 
#				field is set to 0.
#		lum: 	velocity of the brightest spaxel is set 
#				to 0.
#		sig: 	Noralised to the mean velocity of bins with the
#				5% highest velocity dispersion.
# vel_norm: float normalisation factor for finding the rest frame of the galaxy.
# center_bin: int, bin number of the brightest bin.
# common_range: array of min and max of common wavelength range
# galaxy_spectrum: array of spectrum of whole galaxy
# galaxy_continuum: array of spectrum of whole galaxy with emission lines 
#	removed.
# SNRatio: (array) Signal to Noise Ratio of each bin. 
# gas_opt: 	0	(default) No emission lines
#			1	All emission lines set with LOSVD fixed to each other
#			2	All Balmers lines set with the same LOSVD, and all others set with 
#				fixed LOSVD ([OIII] and [NI])
#			3	All emission lines set with independant LOSVD
# independent_components: list components with inderpendant LOSVDs.
#
# Methods:
# add_e_line (line name, line wavelength): adds a new emission line object 
#	to e_line.
# set_spaxels_in_bins (spaxel x-coord, spaxel y-coord, bin membership of spaxel): 
#	sets which bin contains which spaxel.
# absorption_line (absorption line): returns absorption line indice level
# 	from Lick like methods.
	def __init__(self, galaxy, instrument='vimos', opt='kin', overide=False):
		if galaxy == 'ngc1316' and (opt == 'kin' or opt == 'pop') and not overide:
			opt = 'pop_no_Na'
			print 'Setting opt to pop_No_Na'
		self.galaxy, self.instrument, self.opt = galaxy, instrument, opt
		self.MCdir = '%s/Data/%s/analysis/%s/%s/MC' % (cc.base_dir, instrument,
			galaxy, opt)
		self._components = {'stellar':stellar_data(self)}
		self.x,self.y,self.bin_num=np.loadtxt('%s/Data/%s/' % (cc.base_dir, 
			instrument) + 'analysis/%s/%s/setup/voronoi_2d_binning_output.txt' % (
			galaxy, opt), usecols=(0,1,2), dtype=int, skiprows=1, unpack=True)
		
		# Fits files
		opt2 = '_'+self.opt if self.opt != 'kin' and self.opt != 'pop' else ''
		if 'kin' not in self.opt:
			self.abs_fits = fits.getdata('%s/Data/%s/' % (cc.base_dir, instrument)
				+'analysed_fits/%s_absorption_line%s.fits' % (galaxy, opt2), 1)
			self.abs_fits_nomask = fits.getdata('%s/Data/%s/' % (cc.base_dir, 
				instrument) +'analysed_fits/%s_absorption_line_nomask%s.fits' % (
				galaxy,opt2), 1)
			self.emi_fits = fits.getdata('%s/Data/%s/' % (cc.base_dir, instrument)
				+'analysed_fits/%s_emission_line%s.fits' % (galaxy, opt2), 1)
			# if os.path.isfile('%s/Data/%s/' % (cc.base_dir, instrument)
			# 	+'analysed_fits/%s_population%s.fits' % (galaxy, opt2)):
			pop_file = '%s/Data/%s/analysed_fits/%s_population%s.fits' % (
				cc.base_dir, instrument, galaxy, opt2)
			if os.path.exists(pop_file):
				self.pop_fits = fits.getdata(pop_file, 1)

			self._components.update({
				add_brackets(l):emission_data(self, add_brackets(l)) for l in 
				self.emi_fits.columns.names if 
				(l not in ['NO', 'E_VEL', 'E_SIGMA', 'BINSIZE', 'FLUX', 'SNR', 
				'XS', 'YS', 'SIGMA', 'VEL']) and ('anr_' not in l) and 
				('e_' not in l) and ('eqw_' not in l)})
			self._gas = 1
		else:
			self._gas = 0

		self.ste_fits = fits.getdata('%s/Data/%s/' % (cc.base_dir, instrument)
			+'analysed_fits/%s_stellar_kine%s.fits' % (galaxy, opt2), 1)

		if 'kin' in self.opt:
			self.base_fits = self.ste_fits
		else:
			self.base_fits = self.emi_fits

		self.set_spaxels_in_bins(self.x, self.y, self.bin_num)
		self.norm_method = 'lwv'
		self.vel_norm = 0.0
		self.common_range = np.array([])
		self.__threshold__ = 3.0

	@property
	def center(self):
		if self.instrument=='vimos':
			i = 4
		elif self.instrument == 'muse':
			i = 1
		file = '%s/Data/%s/analysis/galaxies.txt' % (cc.base_dir, self.instrument)
		x_cent_gals, y_cent_gals = np.loadtxt(file, usecols=(i,i+1), unpack=True,
			dtype=int, skiprows=1)
		galaxy_gals = np.loadtxt(file, usecols=(0,), skiprows=1, dtype=str)
		i_gal = np.where(galaxy_gals==self.galaxy)[0][0]
		return x_cent_gals[i_gal], y_cent_gals[i_gal]

	@property
	def unbinned_flux(self):
		if self.instrument == 'muse':
			ext = 1
			from errors2_muse import get_dataCubeDirectory
			NAXIS3 = 3681
		elif self.instrument == 'vimos':
			ext = 0
			from errors2 import get_dataCubeDirectory
			NAXIS3 = 1916
		cube_fits, header = fits.getdata(get_dataCubeDirectory(self.galaxy), ext,
			header=True)
		return np.nansum(cube_fits, axis=0)# * NAXIS3 * header['CD3_3']


	# def add_e_line(self, line, wav):
	# 	if line not in self._components.keys():
	# 		self._components[line] = emission_data(self, line, wav)

	def set_spaxels_in_bins(self, x, y, bin_num):
		x,y,bin_num = x.astype(int), y.astype(int), bin_num.astype(int)
		self.number_of_bins = int(max(bin_num) + 1) # zero-base
		self.bin = []
		for i in range(self.number_of_bins):
			self.bin.append(Bin(i, self))
		for i, bin in enumerate(bin_num):
			self.bin[bin].xspaxels.extend([x[i]])
			self.bin[bin].yspaxels.extend([y[i]])

	# Calculate the rest frame for velocities.
	def find_restFrame(self):
		if self.norm_method == "lum":
			self.vel_norm = 0.0
			self.vel_norm = self.components['stellar'].\
				plot['vel'][self.center_bin]
		elif self.norm_method == "lwv":
			self.vel_norm = 0.0
			lwv = self.components['stellar'].plot['vel'].unbinned*self.unbinned_flux
			self.vel_norm = np.nanmean(lwv)*np.sum(self.n_spaxels_in_bin)/ \
				np.nansum(self.unbinned_flux) # average lwv / average flux
		elif self.norm_method == "sig":
			self.vel_norm = 0.0
			s_sort = sorted(np.unique(self.components['stellar'].plot['sigma']))
			c = np.where(self.components['stellar'].plot['sigma'] > s_sort[-6])
			self.vel_norm = np.mean(self.components['stellar'].plot['vel'][c])
		elif self.norm_method == 'lws':
			self.vel_norm = 0.0
			d = np.sqrt((self.xBar-self.bin[self.center[0]].xBar)**2 +
				(self.yBar-self.bin[self.center[1]].yBar)**2)
			lws = self.components['stellar'].plot['sigma'] * self.flux * \
				self.n_spaxels_in_bin
			lws[d > 7] = 0
			s_sort = sorted(lws)
			c = np.where(lws > s_sort[int(-np.ceil(self.number_of_bins*0.05))])[0]
			self.vel_norm = np.nanmean(self.components['stellar'].plot['vel'][c])
		elif self.norm_method == 'disk_fit':
			import disk_fit_functions as dfn
			vel = self.unbin(self.components['stellar'].plot['vel'])
			vel_err = self.unbin(self.components['stellar'].plot['vel'].uncert)
			disk,pars=dfn.disk_fit_exp(vel.copy(),vel_err.copy(),leeway=2., 
				verbose=False)
			self.vel_norm = np.nanmean(disk)
		elif self.norm_method is None:
			self.vel_norm = 0.0

	def absorption_line(self, absorption_line, uncert=False, res=None, 
		instrument='vimos', remove_badpix=False, nomask=False):
		if absorption_line == 'Mg_b' or absorption_line == 'mgb':
			absorption_line = 'Mgb'
		elif absorption_line == 'Hbeta' or absorption_line == 'Hb':
			absorption_line = 'H_beta'
		if res is None:
			if nomask:
				if uncert:
					return self.abs_fits_nomask[absorption_line.upper()], \
						self.abs_fits_nomask['e_'+absorption_line.upper()]
				else:
					return self.abs_fits_nomask[absorption_line]
			else:
				if uncert:
					return self.abs_fits[absorption_line.upper()], \
						self.abs_fits['e_'+absorption_line.upper()]
				else:
					return self.abs_fits[absorption_line]
		else:
			ab = np.zeros(self.number_of_bins)
			ab_uncert = np.zeros(self.number_of_bins)
			for i, bin in enumerate(self.bin):
				if uncert:
					ab[i], ab_uncert[i] = bin.absorption_line(absorption_line, 
						uncert=uncert, res=res, instrument=instrument, 
						remove_badpix=remove_badpix)
				else:
					ab[i] = bin.absorption_line(absorption_line, uncert=uncert,
						res=res, instrument=instrument, remove_badpix=remove_badpix)
			if uncert:
				return ab, ab_uncert
			else:
				return ab

	# Take output some atttribute and return 2D 'unbinned' versions i.e. that value 
	# applied to each spaxel within the bin. Norm keyword is for moment-0 quantities
	def unbin(self, attr, imshow=False):
		if isinstance(attr, str):
			attr = self.__getattribute__(self, attr)

		out = np.full(self.unbinned_flux.shape, np.nan)

		for bin in self.bin:
			out[bin.xspaxels, bin.yspaxels] = attr[bin.bin_number]
		# Transform such that imshow displays as plot_velfield_nointerp does.
		if imshow:
			out = np.rot90(out[::-1,:])
		return out

	def rebin(self, field, flux_weighted=True):
		field = np.rot90(field, 3)[::-1,:] # remove transform applied in unbin
		new = np.full(self.number_of_bins, np.nan)
		for i in range(self.number_of_bins):
			if flux_weighted:
				new[i] = np.average(field[self.bin[i].xspaxels, 
					self.bin[i].yspaxels], 
					weights=self.unbinned_flux[self.bin[i].xspaxels, 
					self.bin[i].yspaxels])
			else:
				new[i] = np.average(field[self.bin[i].xspaxels, 
					self.bin[i].yspaxels])
		return new

	@property
	def center_bin(self):
		# return np.nanargmax(self.flux)
		return np.where((self.xBar == self.center[0]) * \
			(self.yBar == self.center[1]))[0][0]
	
	@property
	def e_components(self):
		return list(self.e_line.keys())
	@property
	def list_components(self):
		return list(self.components.keys())
	@property
	def list_components_no_mask(self):
		return list(self.components_no_mask.keys())

	@property
	def e_line(self):
		return {k:v for k,v in self._components.iteritems() if k!='stellar' 
			and not all(v.mask)}

	@property
	def e_line_no_mask(self):
		return {k:v for k,v in self._components.iteritems() if k!='stellar'}

	@property
	def components(self):
		return {k:v for k,v in self._components.iteritems() if k=='stellar' 
			or not all(v.mask)}

	@property
	def components_no_mask(self):
		return self._components

	@property
	def flux(self):
		return self.base_fits['FLUX']

	# @property
	# def gas_flux(self):
	# 	f = np.zeros(self.number_of_bins)
	# 	u = np.zeros(self.number_of_bins)
	# 	for k,v in self.e_line.iteritems():
	# 		l_flux  = v.flux
	# 		f = np.nansum([f, l_flux], axis=0)
	# 		u = np.nansum([u, l_flux.uncert**2], axis=0)
	# 	u = np.sqrt(u)
	# 	return myArray(f, uncert=u)

	"""
	***************************************************************************
	Need to remove the subtraction from center that occured at some stage in the 
	saving of the fits file.
	****************************************************************************
	"""	
	@property
	def xBar(self):
		return self.center[0] - self.base_fits['XS']
		
	@property
	def yBar(self):	
		return self.center[1] - self.base_fits['YS']

	@property
	def n_spaxels_in_bin(self):
		return self.base_fits['BINSIZE']

	@property
	def galaxy_spectrum(self):
		if self.instrument == 'muse':
			ext = 1
			from errors2_muse import get_dataCubeDirectory
		elif self.instrument == 'vimos':
			ext = 0
			from errors2 import get_dataCubeDirectory
		cube_fits = fits.getdata(get_dataCubeDirectory(self.galaxy), ext)
		return np.nansum(cube_fits, axis=(1,2))

	# @property
	# def galaxy_continuum(self):
	# 	s = np.zeros(len(self.bin[0].continuum))
	# 	for bin in self.bin:
	# 		s += bin.continuum
	# 	return s

	#**** This is uncertainty not carience, but I don't think I use this anywhere
	@property
	def galaxy_varience(self):
		warnings.warn('This returns uncertainty (sigma) NOT varience')
		if self.instrument == 'muse':
			ext = 1
			from errors2_muse import get_dataCubeDirectory
		elif self.instrument == 'vimos':
			ext = 0
			from errors2 import get_dataCubeDirectory
		return np.sqrt(np.sum(fits.getdata(get_dataCubeDirectory(self.galaxy), 
			ext+1)**2, axis=(1,2)))


	@property
	def SNRatio(self):
		return self.base_fits['SNR']

	@property
	def gas(self):
		return self._gas
	@gas.setter
	def gas(self, gas):
		if gas not in [0,1,2,3]:
			raise ValueError('Gas parameter must be an integer between 0 and 3'+
				' (inclusive)')
		else:
			self._gas = gas

	@property
	def independent_components(self):
		if self.gas == 0: 
			return []
		elif self.gas == 1: 
			return ['stellar', 'gas']
		elif self.gas == 2: 
			return ['stellar', 'SF', 'Shocks']
		elif self.gas == 3:
			return self.list_components 

	@property
	def gas_dynamics_SN(self):
		# if self.sauron:
		return self.components['[OIII]5007d'].amp_noise
		# elif self.broad_narrow:
		# 	pass
		# else:
		# 	return np.sort([l.amp_noise for k, l in self.e_line.iteritems()],
		# 		axis=0)[-2,:]


# based on example in http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
class myArray(np.ndarray):
	def __new__(cls, input_array, uncert=None):
		# We first cast to be our class type
		obj = np.asarray(input_array).view(cls)
		# add the new attribute to the created instance
		if uncert is not None or input_array is None:
			obj.uncert = uncert
		return obj

	def __array_finalize__(self, obj):
		self.uncert = getattr(obj, 'uncert', np.full(obj.shape, np.nan))


class _data(object):
# To be used with __init__ method within other classes (stellar_data and 
#	emission_data). NB: these are all 'getter' functions.
# Attributes:
# plot: (dictionary) listing attributes below
# vel (_uncert): (array)
# sigma (_uncert): (array)
# h3 (_uncert): (array)
# h4 (_uncert): (array)
	
	def __init__(self):
		self.moment = 2

	@property
	def plot(self):
		if self.moment == 1:
			return {'vel':self.vel}
		elif self.moment == 2:
			return {'vel':self.vel, 'sigma':self.sigma}
		elif self.moment == 3:
			return {'vel':self.vel, 'sigma':self.sigma, 'h3':self.h3}
		elif self.moment == 4:
			return {'vel':self.vel, 'sigma':self.sigma, 'h3':self.h3, 'h4':self.h4}
			
	# If value is unmeasured.
	def unset(self, attr):
		for i, bin in enumerate(self.__parent__.bin):
			try:
				del bin.components[self.name].__dict__[attr]
			except KeyError:
				pass

	def setkin(self, attr, value):
		# Ensures only acts on kinematics
		if attr in ['vel','sigma','h3','h4']:
			self.moment = max(self.moment, 
				np.where(np.array(['vel','sigma','h3','h4'])==attr)[0][0]+1)
			for i, bin in enumerate(self.__parent__.bin):
				try:
					sav = bin.components[self.name].__dict__[attr].uncert
					bin.components[self.name].__dict__[attr] = myFloat(value[i])
					bin.components[self.name].__dict__[attr].uncert = sav
				except KeyError:
					pass

	def setkin_uncert(self, attr, value):
		# Ensures only acts on kinematics
		if attr in ['vel','sigma','h3','h4']:
			for i, bin in enumerate(self.__parent__.bin):
				try:
					bin.components[self.name].__dict__[attr].uncert = value[i]
				except KeyError:
					pass

	# __getattr__ only called when attribute is not found by other means. 
	def __getattr__(self, attr):
		if attr in ['vel','sigma','h3','h4']:
	# 		m = self.mask_dynamics
	# 		k = np.array([])
	# 		for bin in self.__parent__.bin:
	# 			if not m[bin.bin_number]:
	# 				k = np.append(k, getattr(bin.components[self.name], attr))
	# 			else:
	# 				k = np.append(k, np.nan)

	# 		
	# 		kinematics = myArray(k)
	# 		kinematics.uncert = myArray(
	# 			[getattr(bin.components[self.name], attr, myFloat(np.nan)).uncert 
	# 			if not m[i] else np.nan 
	# 			for i, bin in enumerate(self.__parent__.bin)])

			if self.name == 'stellar':
				kinematics = myArray(self.__parent__.ste_fits[attr.upper()], 
					uncert = myArray(self.__parent__.ste_fits['e_'+attr.upper()]))
			else:
				kinematics = myArray(self.__parent__.emi_fits[attr.upper()], 
					uncert = myArray(self.__parent__.emi_fits['e_'+attr.upper()]))

			# Normalise to rest frame of stars
			if attr == 'vel':
				kinematics -= self.__parent__.vel_norm

			kinematics.unbinned = self.__parent__.unbin(kinematics)
			kinematics.uncert.unbinned = self.__parent__.unbin(kinematics.uncert)

			return kinematics
		else:
			return object.__getattribute__(self,attr)


	


class stellar_data(_data):
# Attributes:
# name: (str) 'stellar'
# mask(_dynamics): (boolean array) all set to False
# Inherited attributes from _data object for reading ppxf fitted stellar
# 	kinematics
	def __init__(self, parent):
		self.__parent__ = parent
		self.name = 'stellar'
		_data.__init__(self)
	@property
	def mask(self):
		return np.array([False]*self.__parent__.number_of_bins)
	@property
	def mask_dynamics(self):
		return self.mask # np.array([False]*self.__parent__.number_of_bins)

	@property
	def age(self):
		return myArray(self.__parent__.pop_fits['age'], 
			uncert=self.__parent__.pop_fits['e_age'])

	@property
	def metalicity(self):
		return myArray(self.__parent__.pop_fits['metalicity'], 
			uncert=self.__parent__.pop_fits['e_metalicity'])

	@property
	def alpha(self):
		return myArray(self.__parent__.pop_fits['alpha'], 
			uncert=self.__parent__.pop_fits['e_alpha'])



class emission_data(_data):
# Attributes:
# flux: array of integrated fitted spectrum
# name: emission line name (will also be directory key in Data object)
# wav: emission line quoted wavelength
# equiv_width: array of equivelent width for each bin
# mask: boolean array of if the emission line is masked in bin
# mask_dynamics: boolean array of if all emission lines fixed to the same LOSVD are 
#	masked in bin
# amp_noise: array of Amplitude to Noise ratio for each bin
#
# Inherited attributes from _data object for reading ppxf fitted kinematics
	def __init__(self, parent, name):#, wav):
		self.__parent__ = parent
		self._name = name
		# self._wav = wav
		_data.__init__(self)

	@property	
	def flux(self):
		p = []
		p_uncert = []
		for bin in self.__parent__.bin:
			try:
				f = bin.e_line[self._name].flux
				p.append(f)
				p_uncert.append(f.uncert)
			except KeyError:
				p.append(np.nan)
				p_uncert.append(np.nan)
		p = myArray(p)
		p.uncert = np.array(p_uncert)
		return p

	@property
	def name(self):
		return self._name

	@property
	def wav(self):
		return self._wav

	@property
	def equiv_width(self):
		# cont = np.array([bin.continuum[np.argmin(np.abs(bin.lam/
		# 	(1 + bin.components['stellar'].vel/c + self.__parent__.vel_norm/c) 
		# 	- self.wav))] for bin in self.__parent__.bin])
		# e = myArray(self.flux/cont)
		# cont_uncert = np.array([
		# 	bin.continuum.uncert[np.argmin(np.abs(bin.lam/
		# 	(1 + bin.components['stellar'].vel/c + self.__parent__.vel_norm/c) 
		# 	- self.wav))] for bin in self.__parent__.bin])
		# e.uncert = np.array(e) * np.sqrt(((self.flux.uncert/self.flux)**2 + 
		# 	(cont_uncert/cont)**2))
		# return e
		return myArray(
			self.__parent__.emi_fits['eqw_'+remove_brackets(self.name)], 
			uncert=self.__parent__.emi_fits[
			'e_eqw_'+remove_brackets(self.name)])

	@property
	def mask(self):
		p = np.ones(self.__parent__.number_of_bins).astype(bool) # all masked
		# if self.__parent__.gas == 3:
		for bin in self.__parent__.bin:
			try:
				p[bin.bin_number] = bin.e_line[self._name].mask
			except KeyError:
				p[bin.bin_number] = True
		return np.array(p)
		# elif self.__parent__.gas == 2:
		# 	pass
		# elif self.__parent__.gas == 1:
		# 	# Check all emission lines for a detection (i.e. mask = False)
		# 	for bin in self.__parent__.bin:
		# 		for c in bin.e_line.keys():
		# 			if not bin.e_line[c].mask:
		# 				p[bin.bin_number] = bin.e_line[c].mask
		# return p

	@property
	def mask_dynamics(self):
		# p = np.ones(self.__parent__.number_of_bins).astype(bool) # all masked
		# if self.__parent__.sauron:
		p = self.__parent__.components['[OIII]5007d'].mask
		# elif self.__parent__.gas == 3:
		# 	p = self.mask
		# elif self.__parent__.gas == 2:
		# 	for k, i in self.__parent__.e_line_no_mask.iteritems():
		# 		if ('H' in k) ^ ('H' in self.name):
		# 			p *= i.mask
		# elif self.__parent__.gas == 1:
		# 	# # mask if less than 2 Balmer or forbidden lines
		# 	# n_Balmer = np.zeros(self.__parent__.number_of_bins, dtype=int)
		# 	# n_forbidden = np.zeros(self.__parent__.number_of_bins, dtype=int)
		# 	# for k, i in self.__parent__.e_line_no_mask.iteritems():
		# 	# 	if 'H' in k:
		# 	# 		n_Balmer += (~i.mask).astype(int)
		# 	# 	else:
		# 	# 		n_forbidden += (~i.mask).astype(int)
		# 	# p = (n_Balmer < 2) + (n_forbidden < 2) 

		# 	# Check all emission lines for at least one detection (i.e. mask = False)
		# 	# for k, i in self.__parent__.e_line_no_mask.iteritems():
		# 	# 	p *= i.mask

		# 	# Check at least 2 detections in emission lines. 
		# 	p = np.zeros(self.__parent__.number_of_bins)
		# 	for k, i in self.__parent__.e_line_no_mask.iteritems():
		# 		p += (~i.mask).astype(int)
		# 	p = p < 2
		return p

	@property
	def amp_noise(self):
		return self.__parent__.emi_fits['anr_'+remove_brackets(self.name)]




class Bin(object):
# Attributes:
# bin_number: (int)
# e_line: (dictionary) of emission_line objects
# temp_weight: (dictionary) of template weights
# lam: (array) x axis of spectrum
# loglam: (array) log of lam
# limLimits: (array) limits of lam
# FWHM_gal: (float) Full-width-half-maximum used for calculating emission lines
# bestfit: (array) bestfit from ppxf
# spectrum: (array) observed spectrum (with cuts from errors2.py)
# noise: (array) observed noise (with cuts from errors2.py)
# contiuum: (array) observed spectrum with emission lines removed
# xspaxels: (array) x-coords of spaxels in bin
# yspaxels: as xspaxels
# xBar: (float) x-coord of center of bin
# yBar: as xBar
# n_spaxels_in_bin: (int) number of spaxels in this bin
# stellar: (_bin_data) object for storing ppxf fitted stellar kinematics
# apweight: (array) weighting of *addative* polynomial used in pPXF fit
# mpweight: (array) weighting of *multiplicative* polynomial used in pPXF fit
# unconvolved_spectrum: (array) The unconvolved spectrum created using the template
#	weights and ap and mp weights.
# unconvolved_lam: (array) wavelength array corresponding to unconvolved_spectrum 
#	(see above)
#
# Methods:
# **DEPRECATED** set_emission_lines (FWHM of observations): requires self._lam is  
# 	set. Uses ppxf utilites to determine which lines are within the observed 
# 	spectrum.
# set_templates (template name, ppxf fitted weighting): Sets the template weights
#	particularly in the emission_line objects. Also now loads the unconvolved 
#	spectrum instead of method below.
# **DEPRECATED** Unconvolved spectrum (wavelength of templates, templates used in 
# 	dictionary): (array) Returns the unconvolved templates. NB: this are not in 
# 	the same wavelength bins as spectrum, bestfit, lam etc, but are binned as the 
# 	stellar templates.
# absorption_line (line name): returns absorption line strength (of LICK style 
#	indicies). Use uncert (boolean ) keyword to return uncertainty as well 
#	(defualt False). 

	def __init__(self, bin_number, parent):
		self.__parent__ = parent
		self.bin_number = int(bin_number)
		# e_line set when set_emission_lines is called
		#self._e_line = {}
		self.components = {'stellar':_bin_data(self)}
		if self.__parent__.opt != 'kin':
			self.components.update({l:emission_line(self, l) for l in 
				self.__parent__.components_no_mask.keys()})

		# temp_weight set when set_templates called
		# self.temp_weight = {}
		# self._lam = np.array([])
		# self.bestfit = np.array([])		# ppxf bestfit
		# self.spectrum = np.array([]) 	# VIMOS Data
		# self.noise = np.array([])		# VIMOS noise
		self.xspaxels = []
		self.yspaxels = []
		# self.apweight = np.array([])
		# self.mpweight = np.array([])
		# self.unconvolved_spectrum = np.array([])
		# self.unconvolved_lam = np.array([])

	@property
	def residual_noise(self):
		residuals = self.spectrum - self.bestfit
		_, residuals, _ = moving_weighted_average(self.lam, residuals, 
			step_size=3., interp=True)
		return np.sqrt(residuals**2 + self.noise**2)

	@property
	def e_line(self):
		return {k:v for k,v in self.components.iteritems() if k!='stellar'}

	@property
	def stellar(self):
		return self.components['stellar']
	
	@property
	def lam(self):
		return np.loadtxt('%s/lambda/%s.dat' % (self.__parent__.MCdir, 
			self.bin_number))
	@property
	def spectrum(self):
		return np.loadtxt('%s/input/%s.dat' % (self.__parent__.MCdir, 
			self.bin_number))
	@property
	def noise(self):
		return np.loadtxt('%s/noise_input/%s.dat' % (self.__parent__.MCdir, 
			self.bin_number))
	@property
	def bestfit(self):
		return np.loadtxt('%s/bestfit/%s.dat' % (self.__parent__.MCdir, 
			self.bin_number))
	# @lam.setter
	# def lam(self, lam):
	# 	self._lam = np.array(lam)
	# 	if len(self.__parent__.common_range) != 0:
	# 		self.__parent__.common_range[0] = max(
	# 			self.__parent__.common_range[0],min(lam))
	# 		self.__parent__.common_range[1] = min(
	# 			self.__parent__.common_range[1],max(lam))
	# 	else:
	# 		self.__parent__.common_range = np.array([min(lam), max(lam)])

	# @property
	# def loglam(self):
	# 	return np.log(self._lam)
	# @loglam.setter
	# def loglam(self, loglam):
	# 	self._lam = np.exp(loglam)

	# @property
	# def lamLimits(self):
	# 	return np.array([self._lam[0],self._lam[-1]])

	# @property
	# def continuum(self):
	# 	# NB: Masks not used
	# 	c = myArray(self.spectrum - np.nansum([line.spectrum 
	# 		for key, line in self.e_line.iteritems()],axis=0))
	# 	c.uncert = np.sqrt(self.noise**2 + np.nansum([line.uncert_spectrum**2
	# 		for key, line in self.e_line.iteritems()],axis=0))
	# 	return c

	@property
	def flux(self):
		# NB: only calculated for the common wavelength range.
		return self.__parent__.flux[self.bin_number]

	@property
	def xBar(self):
		return self.__parent__.xBar[self.bin_number]

	@property
	def yBar(self):
		return self.__parent__.yBar[self.bin_number]

	@property
	def n_spaxels_in_bin(self):
		return len(self.xspaxels)

	# def set_templates(self, name, weight):
	# 	weight = weight.astype(float)
	# 	## Solving an error in the Uni system with loading large numbers of files  
	# 	## from the home directory - they have to come from the Data partition 
	# 	## instead.
	# 	if cc.getDevice() == 'uni':
	# 		files = glob('%s/Data/idl_libraries/ppxf/MILES_library/'%(cc.base_dir) +
	# 			'm0[0-9][0-9][0-9]V')
	# 	else:
	# 		files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (
	# 			cc.home_dir))
	# 	wav = np.loadtxt(files[0], usecols=(0,), unpack=True)
	# 	# ***** NB: spec is binned as the stellar templates are binned, NOT as 
	# 	# self.spectrum, bestfit etc i.e. wav != self.lam  *****
	# 	a = [min(np.where(wav>=self.lam[0])[0]), max(np.where(wav<=self.lam[-1])[0])]
	# 	unc_spec = np.zeros(a[1]-a[0])

	# 	for i,n in enumerate(name):
	# 		self.temp_weight[n] = weight[i]
	# 		if not n.isdigit() and weight[i] != 0:
	# 			self.components[n].weight = weight[i]
	# 			self.__parent__.add_e_line(n, self.e_line[n].wav)
	# 		elif n.isdigit() and weight[i] !=0:
	# 			template =  np.loadtxt(files[int(n)], usecols=(1,), unpack=True)
	# 			unc_spec += template[a[0]:a[1]]*weight[i]
	# 	if len(self.mpweight)!=0:
	# 		unc_spec *= np.polynomial.legendre.legval(np.linspace(-1,1,
	# 			len(unc_spec)), np.append(1, self.mpweight))
	# 	if len(self.apweight)!=0:
	# 		unc_spec += np.polynomial.legendre.legval(np.linspace(-1,1,
	# 			len(unc_spec)), self.apweight)
	# 	self.unconvolved_spectrum = unc_spec
	# 	self.unconvolved_lam = wav[a[0]:a[1]]

	# def absorption_line(self, absorption_line, uncert=False, res=None,
	# 	instrument=None, remove_badpix=False):
	# 	convolved = self.bestfit - np.nansum([line.spectrum for key, 
	# 		line in self.e_line.iteritems()], axis=0)
	# 	lam = self.lam/(1+self.components['stellar'].vel/c)
	# 	continuum = np.array(self.continuum)

	# 	if uncert:
	# 		noise = np.array(self.noise)

	# 	if remove_badpix:
	# 		# bad_pix = np.where(continuum < 0)[0]
	# 		bad_pix = np.where(self.noise > np.nanmean(self.noise) 
	# 			+ 3*np.std(self.noise))
	# 		bad_pix = np.unique([[b-1,b,b+1] for b in bad_pix])
	# 		bad_pix = bad_pix[(bad_pix >= 1) * (bad_pix <= len(lam) - 2)]
			
	# 		start = np.array([list(g)[0] for _, g in 
	# 			groupby(bad_pix, lambda n, c=count(): n-next(c))]) - 1
	# 		end = np.array([list(g)[-1] for _, g in 
	# 			groupby(bad_pix, lambda n, c=count(): n-next(c))]) + 1

	# 		for index in zip(start, end):
	# 			interp = interp1d([lam[index[0]], lam[index[1]]], 
	# 				[continuum[index[0]], continuum[index[1]]])
	# 			continuum[index[0]+1:index[1]] = interp(lam[np.arange(
	# 				index[0]+1, index[1])])
	# 			if uncert:
	# 				interp = interp1d([lam[index[0]], lam[index[1]]], 
	# 					[noise[index[0]], noise[index[1]]])
	# 				noise[index[0]+1:index[1]] = interp(lam[np.arange(
	# 				index[0]+1, index[1])])

	# 	temp_res = 2.5 # A (FWHM); Resoluton of Miles templates
	# 	if res is not None:
	# 		if instrument == 'vimos':
	# 			instr_res = 3 # A (FWHM)
	# 		elif instrument == 'muse':
	# 			instr_res = 2.3 # A (FWHM)
	# 		elif instrument == 'sauron':
	# 			instr_res = 4.2 # A (FWHM)
	# 		elif instrument is None:
	# 			raise ValueError('If keyword res is set, so too must the '+
	# 				'instrument keyword')

	# 		if instr_res < res:
	# 			sig_pix = np.sqrt(res**2 - instr_res**2) / 2.355 \
	# 				/ np.median(np.diff(lam))
	# 			continuum = gaussian_filter1d(continuum, sig_pix)
	# 		conv_res = max((instr_res, temp_res))
	# 	else:
	# 		conv_res = temp_res
	# 		continuum = np.array(self.continuum)
			
	# 	if conv_res < temp_res:
	# 		sig_pix = np.sqrt(temp_res**2 - conv_res**2) / 2.355 \
	# 			/ np.median(np.diff(lam))
	# 		convolved = gaussian_filter1d(convolved, sig_pix)
	# 	if uncert:
	# 		return absorption(absorption_line, lam, continuum, 
	# 			noise=self.noise, unc_lam=self.unconvolved_lam, 
	# 			unc_spec=self.unconvolved_spectrum, conv_spec=convolved)
	# 	else:
	# 		return absorption(absorption_line, lam, continuum, 
	# 			unc_lam=self.unconvolved_lam, 
	# 			unc_spec=self.unconvolved_spectrum, conv_spec=convolved)

class myFloat(float):
# Required to add attributes to float object
	def __init__(self, value):
		float.__init__(self)
		# self._uncert = np.nan
		self.uncert = np.nan



class _bin_data(object):
# Attributes:
# vel (.uncert): float
# sigma (.uncert): float
# h3 (.uncert): float
# h4 (.uncert): float
	def __init__(self, parent):
		self.__parent__=parent
		# self.vel = myFloat(np.nan)
		# self.sigma = myFloat(np.nan)
		self.h3 = myFloat(np.nan)
		self.h4 = myFloat(np.nan)
		# self.vel.uncert = np.nan
		# self.sigma.uncert = np.nan
		self.h3.uncert = np.nan
		self.h4.uncert = np.nan

	@property
	def vel(self):
		if hasattr(self, 'name'):
			k = myFloat(self.__parent__.__parent__.components[self.name].plot[
				'vel'][self.__parent__.bin_number])
			k.uncert = self.__parent__.__parent__.components[self.name].plot[
				'vel'].uncert[self.__parent__.bin_number]
		else: # assume stellar object
			k = myFloat(self.__parent__.__parent__.components['stellar'].plot[
				'vel'][self.__parent__.bin_number])
			k.uncert = self.__parent__.__parent__.components['stellar'].plot[
				'vel'].uncert[self.__parent__.bin_number]
		return k
	@property
	def sigma(self):
		if hasattr(self, 'name'):
			k = myFloat(self.__parent__.__parent__.components[self.name].plot[
				'sigma'][self.__parent__.bin_number])
			k.uncert = self.__parent__.__parent__.components[self.name].plot[
				'sigma'].uncert[self.__parent__.bin_number]
		else: # assume stellar object
			k = myFloat(self.__parent__.__parent__.components['stellar'].plot[
				'sigma'][self.__parent__.bin_number])
			k.uncert = self.__parent__.__parent__.components['stellar'].plot[
				'sigma'].uncert[self.__parent__.bin_number]
		return k
	# @property
	# def vel(self):
	# 	return self._vel
	# @property
	# def sigma(self):
	# 	return self._sigma
	# @property
	# def h3(self):
	# 	return self._h3
	# @property
	# def h4(self):
	# 	return self._h4
	# @vel.setter
	# def vel(self, value):
	# 	self._vel = value
	# @sigma.setter
	# def sigma(self, value):
	# 	self._sigma = value
	# @h3.setter
	# def h3(self, value):
	# 	self._h3 = value
	# @h4.setter
	# def h4(self, value):
	# 	self._h4 = value



class emission_line(_bin_data):
# Attributes:
# name: str, name of emission line
# _weight: float, weighting from ppxf fit
# wav: wavelength of emission
# flux: float, integrated fitted spectrum.
# spectrum: array if not masked of fitted spectrum
# spectrum_nomask: as spectrum overides masking
# AmpNoi: float, amplitude of fitted spectrum vs noise from Bin object
# mask: boolean: True if emission line is masked in this bin
# Inherited attributes from _bin_data object for storing ppxf fitted kinematics
	def __init__(self, parent, name):#, wav, _spectrum):
		self.__parent__ = parent
		self.name = str(name)
		# self.weight = 0.0
		# self._spectrum = np.array(_spectrum)
		# self.uncert_spectrum = np.zeros(len(self._spectrum))*np.nan
		# self.wav = wav
		# _bin_data.__init__(self, self.__parent__)

		# self.vel = myFloat(self.__parent__.__parent__.emi_fits['VEL'])
		# self.sigma = myFloat(self.__parent__.__parent__.emi_fits['SIGMA'])
		# self.h3 = myFloat(np.nan)
		# self.h4 = myFloat(np.nan)
		# self.vel.uncert = self.__parent__.__parent__.emi_fits['e_VEL']
		# self.sigma.uncert = self.__parent__.__parent__.emi_fits['e_SIMGA']
		# self.h3.uncert = np.nan
		# self.h4.uncert = np.nan
		self.h3 = myFloat(np.nan)
		self.h4 = myFloat(np.nan)

	
		
	@property
	def flux(self):
		f = myFloat(self.__parent__.__parent__.emi_fits[remove_brackets(self.name)]
			[self.__parent__.bin_number])
		f.uncert = self.__parent__.__parent__.emi_fits['e_'
			+ remove_brackets(self.name)][self.__parent__.bin_number]
		return f

	@property
	def spectrum(self):
		matrix = np.loadtxt("%s/bestfit/matrix/%d.dat" % (
			self.__parent__.__parent__.MCdir, self.__parent__.bin_number), 
		dtype=str)
		i = np.where(matrix[:,0] == self.name)[0][0]
		spec = matrix[i,1:].astype(float)

		temp_name, temp_weight = np.loadtxt("%s/temp_weights/%d.dat" % (
			self.__parent__.__parent__.MCdir, self.__parent__.bin_number), 
			unpack=True, dtype=str)
		i = np.where(temp_name == self.name)[0][0]
		return spec * float(temp_weight[i])

	@property
	def uncert_spectrum(self):
		return np.loadtxt('%s/gas_uncert_spectrum/%s/%i.dat' % (
			self.__parent__.__parent__.MCdir, self.name, 
			self.__parent__.bin_number))



	# @property
	# def spectrum(self):
	# 	if not self.mask:
	# 		return self._spectrum*self.weight
	# 	else:
	# 		return self._spectrum*np.nan

	# @property
	# def spectrum_nomask(self):
	# 	return self._spectrum*self.weight

	
	@property
	def AmpNoi(self):
		# return max(self.spectrum_nomask)/np.median(
		# 	self.__parent__.residual_noise[
		# 	(self.__parent__.loglam > np.log(self.wav) + (self.vel - 300)/c) *
		# 	(self.__parent__.loglam < np.log(self.wav) + (self.vel + 300)/c)])
		return self.__parent__.__parent__.emi_fits[
			'anr_'+remove_brackets(self.name)][self.__parent__.bin_number]

	@property
	def mask(self):
		if 'OIII' in self.name:
			return self.AmpNoi < 4
		elif 'NI' in self.name:
			return (self.AmpNoi < 4) + self.__parent__.e_line['Hbeta'].mask
		else:
			return (self.AmpNoi < 3) + \
				self.__parent__.e_line['[OIII]5007d'].mask
	

# Find the propergated uncertainty from numpy.trapz().
def trapz_uncert(uncert, x=None):
	if x is None:
		x = np.arange(len(uncert))
	return np.sqrt(np.nansum([(x[i+1]-x[i])**2 * (uncert[i+1]**2 + uncert[i]**2)
		 for i in range(len(uncert)-1)]))/2
