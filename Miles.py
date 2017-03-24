# Routine to create class for comparison to Miles models.
from astropy.io import fits
import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

class miles(object):
	def __init__(self, galaxy):
		self.f = fits.open('%s/Data/vimos/cubes/%s.cube.combined.corr.fits' % (
			cc.base_dir, galaxy))
		s = self.f[0].data.shape

		self.gal_spec = np.nansum(self.f[0].data, axis=(1,2))
		self.gal_noise = np.sqrt(np.nansum(self.f[1].data**2, axis=(1,2)))

		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
		# different data types need to be read separetly
		z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
			unpack=True, skiprows=1, usecols=(1,2,3,4,5))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		i_gal = np.where(galaxy_gals==galaxy)[0][0]
		self.vel = vel_gals[i_gal]
		self.sig = sig_gals[i_gal]
		self.z = z_gals[i_gal]

		self.lam = np.arange(s[0])*self.f[0].header['CDELT3'] + self.f[0].header['CRVAL3']

	def get_spec(self):
		return self.gal_spec, self.gal_noise