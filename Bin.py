## ==================================================================
## 		Emission lines
## ==================================================================
## warrenj 20160725 An object to hold emission lines. 

import numpy as np
import ppxf_util as util

class Bins_list(object):
	def __init__(self,number_of_bins):
		self.number_of_bins = number_of_bins
		self.components = []

		bins = []
		for i in range(number_of_bins):
			bins.append(Bin(i))
		self.bin = bins

	def components1(self, line):
		if line not in self.components:
			self.components.append[line]

	@property
	def e_lines(self, line):

		def flux(self):
			return [self.bin[i].e_line[line].flux for i in range(self.number_of_bins)]

		

	

class emission_line(object):
	def __init__(self, name, wav, spectrum, lam):
		self.name = name
		self.spectrum = spectrum
		self.wav = wav
		self.weight = 0.0
		self.lam = lam
		
	@property
	def flux(self):
		return np.trapz(self.spectrum, x=self.lam)
	



		
class Bin(Bins_list):
# Attributes:
# bin_number: int
# e_line: dictionary of emission_line objects
# temp_weight: dictionary of template weights
# lam: (array) x axis of spectrum
# loglam: (array) log of lam
# limLimits: (array) limits of lam
# FWHM_gal: (float) Full-width-half-maximum used for calculating emission lines
# bestfit: (array) bestfit from ppxf

	def __init__(self, bin_number):
		self.bin_number = int(bin_number)
		# e_line set when set_emission_lines is called
		self.e_line = {}
		# temp_weight set when set_templates called
		self.temp_weight = {}
		self._lam = []
		self.bestfit = []

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
	
	


	
	def set_emission_lines(self, FWHM_gal):
	# Sets emission lines
		self.FWHM_gal = float(FWHM_gal)

		line_spectrums, line_names, line_wavs = util.emission_lines(
			self.loglam, self.lamLimits, FWHM_gal, quiet=True)

		for l in range(len(line_wavs)):
			line = emission_line(line_names[l], line_wavs[l], 
				line_spectrums[:,l].flatten(), self.lam)
			self.e_line[line_names[l]] = line

	def set_templates(self, name, weight):
		weight = weight.astype(float)
		for i in range(len(name)):
			self.temp_weight[name[i]] = weight[i]
			if not name[i].isdigit():
				self.e_line[name[i]].weight = weight[i]
				super(Bin, self).components1(name[i])



