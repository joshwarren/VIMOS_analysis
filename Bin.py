## ==================================================================
## 		Emission lines
## ==================================================================
## warrenj 20160725 An object to hold emission lines. 

import numpy as np
import ppxf_util as util

class Data(object):
# Attributes:
# number_of_bins: int
# bin: array of bin objects (see below) - holds most of the data - is a 
#	'setter' class
# e_line: dictionary of emission_data objects (see below) - is a 'getter' class
# e_components: string array of emission lines
# unbinned_flux: float array of unbinned flux in 2D array
# flux: float array of total integrated flux (observed - not fitted)
# xBar: float array of central point's x-coord for each bin
# yBar: as xBar
# n_spaxels_in_bin: int array of number of spaxels in each bin
#
# Methods:
# add_e_line (line name, line wavelength): adds a new emission line object 
#	to e_line.
# set_spaxels_in_bins (spaxel x-coord, spaxel y-coord, bin membership of spaxel): 
#	sets which bin contains which spaxel.

	def __init__(self, xyb_turple):
		x,y,bin_num = xyb_turple
		self.x,self.y,self.bin_num=x,y,bin_num
		self.set_spaxels_in_bins(x, y, bin_num)
		self.unbinned_flux = [[] for i in range(2)] # 2D blank array
		self.components = {'stellar':stellar_data(self)}


	def add_e_line(self, line, wav):
		if line not in self.e_components:
			self.components[line] = emission_data(self, line, wav)

	def set_spaxels_in_bins(self, x, y, bin_num):
		x,y,bin_num = x.astype(int), y.astype(int), bin_num.astype(int)
		self.number_of_bins = int(max(bin_num) + 1) # zero-base
		self.bin = []
		for i in range(self.number_of_bins):
			self.bin.append(Bin(i, self))
		for i, bin in enumerate(bin_num):
			self.bin[bin].xspaxels.extend([x[i]])
			self.bin[bin].yspaxels.extend([y[i]])

	@property
	def e_components(self):
		return list(self.e_line.keys()) 
	@property
	def list_components(self):
		return list(self.components.keys()) 

	@property
	def e_line(self):
		return {k:v for k,v in self.components.iteritems() if k!='stellar'}
	

	@property
	def flux(self):
		return [bin.flux for bin in self.bin]

	@property
	def xBar(self):
		return [bin.xBar for bin in self.bin]	
	@xBar.setter
	def xBar(self, xBar):
		for i, x in enumerate(xBar):
			self.bin[i].xBar = x

	@property
	def yBar(self):	
		return [bin.yBar for bin in self.bin]
	@yBar.setter
	def yBar(self, yBar):
		for i, y in enumerate(yBar):
			self.bin[i].yBar = y

	@property
	def n_spaxels_in_bin(self):
		return [bin.n_spaxels_in_bin for bin in self.bin]



class _data(object):
# Attributes:
# plot: (dictionary) listing attributes below
# vel (_uncert): (array)
# sigma (_uncert): (array)
# h3 (_uncert): (array)
# h4 (_uncert): (array)
	
	def __init__(self, parent, species):
		self.__parent__ = parent
		self.species = species
		self.plot = {'vel':self.vel, 'sigma':self.sigma, 'h3':self.h3, 
			'h4':self.h4, 'vel_uncert':self.vel_uncert, 
			'sigma_uncert':self.sigma_uncert, 'h3_uncert':self.h3_uncert, 
			'h4_uncert':self.h4_uncert}	


	@property
	def vel(self):
		return [bin.components[self.species].vel if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@vel.setter
	def vel(self, vel):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].vel = vel[i]

	@property
	def vel_uncert(self):
		return [bin.components[self.species].vel_uncert if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@vel_uncert.setter
	def vel_uncert(self, vel_uncert):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].vel_uncert = vel_uncert[i]

	@property
	def sigma(self):
		return [bin.components[self.species].sigma if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@sigma.setter
	def sigma(self, sigma):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].sigma = sigma[i]

	@property
	def sigma_uncert(self):
		return [bin.components[self.species].sigma_uncert if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@sigma_uncert.setter
	def sigma_uncert(self, sigma_uncert):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].sigma_uncert = sigma_uncert[i]

	@property
	def h3(self):
		return [bin.components[self.species].h3 if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@h3.setter
	def h3(self, h3):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].h3 = h3[i]

	@property
	def h3_uncert(self):
		return [bin.components[self.species].h3_uncert if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@h3_uncert.setter
	def h3_uncert(self, h3_uncert):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].h3_uncert = h3_uncert[i]

	@property
	def h4(self):
		return [bin.components[self.species].h4 if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@h4.setter
	def h4(self, h4):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].h4 = h4[i]

	@property
	def h4_uncert(self):
		return [bin.components[self.species].h4_uncert if not 
			self.__parent__.mask[i] else np.nan 
			for i, bin in enumerate(self.__parent__.__parent__.bin)]
	@h4_uncert.setter
	def h4_uncert(self, h4_uncert):
		for i, bin in enumerate(self.__parent__.__parent__.bin):
			bin.components[self.species].h4_uncert = h4_uncert[i]

	


class stellar_data(Data, _data):
	def __init__(self, parent):
		self.__parent__ = parent
		self.species = 'stellar'
		_data.__init__(self, self, self.species)
	@property
	def mask(self):
		return [False]*self.__parent__.number_of_bins


	
class emission_data(Data, _data):
# Attributes:
# flux: array of integrated fitted spectrum
# name: emission line name (will also be directory key in Data object)
# wav: emission line quoted wavelength
# equiv_width: array of equivelent width for each bin
# mask: boolean array of if the emission line is masked in bin
# vel: (array)
# sigma: (array)
# h3: (array)
# h4: (array)
	def __init__(self, parent, name, wav):
		self.__parent__ = parent
		self._name = name
		self._wav = wav
		_data.__init__(self, self.__parent__, self._name)

	@property	
	def flux(self):
		p = []
		for bin in self.__parent__.bin:
			try:
				p.append(bin.e_line[self._name].flux)
			except KeyError:
				p.append(np.nan)
		return np.array(p)

	@property
	def name(self):
		return self._name

	@property
	def wav(self):
		return self._wav

	@property
	def equiv_width(self):
		return np.array(self.flux/[bin.continuum[np.argmin(np.abs(bin.lam-self.wav))] for bin in self.__parent__.bin])

	@property
	def mask(self):
		p = []
		for bin in self.__parent__.bin:
			try:
				p.append(bin.e_line[self._name].mask)
			except KeyError:
				p.append(True)
		return np.array(p)

	# @property
	# def vel(self):
	# 	return [bin.e_line[self._name].vel for bin in self.__parent__.bin]
	# @vel.setter
	# def vel(self, vel):
	# 	for i, bin in enumerate(self.__parent__.bin):
	# 		bin.e_line[self._name].vel = vel[i]
	# @property
	# def sigma(self):
	# 	return [bin.e_line[self._name].sigma for bin in self.__parent__.bin]
	# @sigma.setter
	# def sigma(self, sigma):
	# 	for i, bin in enumerate(self.__parent__.bin):
	# 		bin.e_line[self._name].sigma = sigma[i]
	# @property
	# def h3(self):
	# 	return [bin.e_line[self._name].h3 for bin in self.__parent__.bin]
	# @h3.setter
	# def h3(self, h3):
	# 	for i, bin in enumerate(self.__parent__.bin):
	# 		bin.e_line[self._name].h3 = h3[i]
	# @property
	# def h4(self):
	# 	return [bin.e_line[self._name].h4 for bin in self.__parent__.bin]
	# @h4.setter
	# def h4(self, h4):
	# 	for i, bin in enumerate(self.__parent__.bin):
	# 		bin.e_line[self._name].h4 = h4[i]






class Bin(Data):
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
# xspaxels: (array) x-coords of spaxels in bin
# yspaxels: as xspaxels
# xBar: (float) x-coord of center of bin
# yBar: as xBar
# n_spaxels_in_bin: (int) number of spaxels in this bin
# stellar: (stellar) object
#
# Methods:
# set_emission_lines (FWHM of observations): requires self._lam is set. Uses ppxf
#	utilites to determine which lines are within the observed spectrum.
# set_templates (template name, ppxf fitted weighting): Sets the template weights
#	particularly in the emission_line objects. 

	def __init__(self, bin_number, parent):
		self.__parent__ = parent
		self.bin_number = int(bin_number)
		# e_line set when set_emission_lines is called
		#self._e_line = {}
		self._components = {'stellar':_bin_data(self)}
		# temp_weight set when set_templates called
		self.temp_weight = {}
		self._lam = np.array([])
		self.bestfit = np.array([])		# ppxf bestfit
		self.spectrum = np.array([]) 	# VIMOS Data
		self.noise = np.array([])		# VIMOS noise
		self.xspaxels = []
		self.yspaxels = []
		self._xBar = np.nan
		self._yBar = np.nan

	@property
	def components(self):
		return self._components
	

	@property
	def e_line(self):
		return {k:v for k,v in self.components.iteritems() if k!='stellar'}

	@property
	def stellar(self):
		return self.components['stellar']
	
	@property
	def lam(self):
		return self._lam
	@lam.setter
	def lam(self, lam):
		self._lam = np.array(lam)

	@property
	def loglam(self):
		return np.log(self._lam)
	@loglam.setter
	def loglam(self, loglam):
		self._lam = np.exp(loglam)

	@property
	def lamLimits(self):
		return [self._lam[0],self._lam[-1]]

	@property
	def continuum(self):
		# NB: Masks not used as emission lines still used in bestfit.
		c = np.array(self.bestfit)
		for key,line in self.e_line.iteritems():
			c -= line.spectrum
		return c

	@property
	def flux(self):
		return np.trapz(self.spectrum, x=self.lam)

	@property
	def xBar(self):
		return self._xBar
	@xBar.setter
	def xBar(self, xBar):
		self._xBar = xBar

	@property
	def yBar(self):
		return self._yBar
	@yBar.setter
	def yBar(self, yBar):
		self._yBar = yBar
	

	@property
	def n_spaxels_in_bin(self):
		return len(self.xspaxels)
	

	def set_emission_lines(self, FWHM_gal):
	# Sets emission lines
		self.FWHM_gal = float(FWHM_gal)

		line_spectrums, line_names, line_wavs = util.emission_lines(
			self.loglam, self.lamLimits, FWHM_gal, quiet=True)

		for i in range(len(line_wavs)):
			line = emission_line(self, line_names[i], line_wavs[i], 
				line_spectrums[:,i].flatten())
			self.components[line_names[i]] = line
			if line_names[i] not in self.__parent__.e_components:
				self.__parent__.add_e_line(line_names[i], line_wavs[i])

	def set_templates(self, name, weight):
		weight = weight.astype(float)
		for i in range(len(name)):
			self.temp_weight[name[i]] = weight[i]
			if not name[i].isdigit():
				self.components[name[i]].weight = weight[i]
				if name[i] not in self.__parent__.e_components:
					self.__parent__.add_e_line(name[i], self.e_line[name[i]].wav)


class _bin_data(object):
# Attributes:
# vel (_uncert): float
# sigma (_uncert): float
# h3 (_uncert): float
# h4 (_uncert): float
	def __init__(self, parent):
		self.__parent__=parent
		self.vel = np.nan
		self.sigma = np.nan
		self.h3 = np.nan
		self.h4 = np.nan
		self.vel_uncert = np.nan
		self.sigma_uncert = np.nan
		self.h3_uncert = np.nan
		self.h4_uncert = np.nan


class emission_line(Bin, _bin_data):
# Attributes:
# name: str, name of emission line
# weight: float, weighting from ppxf fit
# wav: wavelength of emission
# flux: float, integrated fitted spectrum
# spectrum: array if not masked of fitted spectrum
# AmpNoi: float, amplitude of fitted spectrum vs noise from Bin object
# mask: boolean: True if emission line is masked in this bin
# vel: float
# sigma: float
# h3: float
# h4: float
	def __init__(self, parent, name, wav, spectrum):
		self.__parent__ = parent
		self.name = str(name)
		self.weight = 0.0
		self._spectrum = np.array(spectrum)
		self.wav = wav
		self.__threshold__ = 4.0
		_bin_data.__init__(self, self.__parent__)

		
	@property
	def flux(self):
		if not self.mask:
			return np.trapz(self.spectrum, x=self.__parent__.lam)
		else:
			return np.nan

	@property
	def spectrum(self):
		if not self.mask:
			return self._spectrum*self.weight
		else:
			return self._spectrum*np.nan

	@spectrum.setter
	def spectrum(self, spectrum):
		self._spectrum = np.array(spectrum)
		
	@property
	def AmpNoi(self):
		return max(self._spectrum)/self.__parent__.noise[np.argmin(np.abs(
			self.__parent__.lam-self.wav))]

	@property
	def mask(self):
		return self.AmpNoi < self.__threshold__

