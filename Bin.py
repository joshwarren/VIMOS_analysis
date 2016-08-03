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
# components: array of emission lines

	def __init__(self,number_of_bins):
		self.number_of_bins = number_of_bins
		self.e_line = {}

		bins = []
		for i in range(number_of_bins):
			bins.append(Bin(i, self))
		self.bin = bins


	def add_e_line(self, line, wav):
		if line not in self.components:
			self.e_line[line] = emission_data(self, line, wav)

	@property
	def components(self):
		return list(self.e_line.keys()) 

	@property
	def flux(self):
		return [bin.flux for bin in self.bin]

	
class emission_data(Data):
# Attributes:
# flux: array of integrated fitted spectrum
# name: emission line name (will also be directory key in Data object)
# wav: emission line quoted wavelength
# equiv_width: array of equivelent width for each bin
# mask: boolean array of if the emission line is masked in bin
	def __init__(self, parent, name, wav):
		self.__parent__ = parent
		self._name = name
		self._wav = wav

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
		#if self._name == 'Hdelta':
			#print [bin.continuum[np.argmin(np.abs(bin.lam-self.wav))] for bin in self.__parent__.bin]
			#print ''
			#print np.where(self.flux < 0)[0]
			#print ''
			#print np.where (self.flux/[bin.continuum[np.argmin(np.abs(bin.lam-self.wav))] for bin in self.__parent__.bin] <0)[0]
			#print ''
			#print ''
			#print ''
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


class Bin(Data):
# Attributes:
# bin_number: int
# e_line: dictionary of emission_line objects
# temp_weight: dictionary of template weights
# lam: (array) x axis of spectrum
# loglam: (array) log of lam
# limLimits: (array) limits of lam
# FWHM_gal: (float) Full-width-half-maximum used for calculating emission lines
# bestfit: (array) bestfit from ppxf

	def __init__(self, bin_number, parent):
		self.__parent__ = parent
		self.bin_number = int(bin_number)
		# e_line set when set_emission_lines is called
		self.e_line = {}
		# temp_weight set when set_templates called
		self.temp_weight = {}
		self._lam = np.array([])
		self.bestfit = np.array([])		# ppxf bestfit
		self.spectrum = np.array([]) 	# VIMOS Data
		self.noise = np.array([])		# VIMOS noise

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
	
	def set_emission_lines(self, FWHM_gal):
	# Sets emission lines
		self.FWHM_gal = float(FWHM_gal)

		line_spectrums, line_names, line_wavs = util.emission_lines(
			self.loglam, self.lamLimits, FWHM_gal, quiet=True)

		for i in range(len(line_wavs)):
			line = emission_line(self, line_names[i], line_wavs[i], 
				line_spectrums[:,i].flatten())
			self.e_line[line_names[i]] = line
			if line_names[i] not in self.__parent__.components:
				self.__parent__.add_e_line(line_names[i], line_wavs[i])

	def set_templates(self, name, weight):
		weight = weight.astype(float)
		for i in range(len(name)):
			self.temp_weight[name[i]] = weight[i]
			if not name[i].isdigit():
				self.e_line[name[i]].weight = weight[i]
				if name[i] not in self.__parent__.components:
					self.__parent__.add_e_line(name[i], self.e_line[name[i]].wav)


class emission_line(Bin):
# Attributes:
# flux: float, integrated fitted spectrum
# spectrum: array if not masked of fitted spectrum
# AmpNoi: float, amplitude of fitted spectrum vs noise from Bin object
# maks: boolean: True if emission line is masked in this bin
	def __init__(self, parent, name, wav, spectrum):
		self.__parent__ = parent
		self.name = name
		self.weight = 0.0
		self._spectrum = np.array(spectrum)
		self.wav = wav
		self.__threshold__ = 4.0
		
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

