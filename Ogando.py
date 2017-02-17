## warrenj 20170209
## Routine to move Ogando methods into.
from checkcomp import checkcomp
cc = checkcomp()
from astropy.io import fits
import numpy as np
from tools import funccontains, slit, get_slit


class ogando(object):
	def __init__(self, galaxy, slit_h=4.1, slit_w=2.5, slit_pa=30, debug=False):
		self.debug = debug
## ----------============== Load galaxy info ================---------
		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
		# different data types need to be read separetly
		z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
			unpack=True, skiprows=1, usecols=(1,2,3,4,5))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		i_gal = np.where(galaxy_gals==galaxy)[0][0]
		self.vel = vel_gals[i_gal]
		self.sig = sig_gals[i_gal]
		self.z = z_gals[i_gal]

		self.f = fits.open('%s/Data/vimos/cubes/%s.cube.combined.corr.fits' % (cc.base_dir, 
			galaxy))
		self.lam = np.arange(self.f[0].header['NAXIS3'])*self.f[0].header['CDELT3'] + self.f[0].header['CRVAL3']
## ----------============ Find data within slit =============---------

		x = (np.arange(self.f[0].header['NAXIS1']) * self.f[0].header['CDELT1']).repeat(
			self.f[0].header['NAXIS2'])
		y = np.tile(np.arange(self.f[0].header['NAXIS2']) * self.f[0].header['CDELT2'],
			self.f[0].header['NAXIS1'])

		slit_x, slit_y = get_slit(galaxy, slit_h, slit_w, slit_pa)

		frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
		if self.debug:
			frac = funccontains(slit, (slit_x,slit_y), x=x, y=y).contains.astype(int)
		else:
			frac = funccontains(slit, (slit_x,slit_y), x=x, y=y).fraction
		frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(self.f[0].header['NAXIS2']),
			np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].header['NAXIS1'])] = frac

		cube = self.f[0].data
		cube[self.f[3].data==1] = 0
		noise_cube = self.f[1].data
		noise_cube[self.f[3].data==1] = 0.000000001
		cube[~np.isfinite(noise_cube)] = 0
		noise_cube[~np.isfinite(noise_cube)] = 0.000000001
	
		self.gal_spec = np.einsum('ijk,jk->i', cube, frac_image)
		self.gal_noise = np.sqrt(np.einsum('ijk,jk->i', noise_cube**2, frac_image**2))

	def get_spec(self):
		return self.gal_spec, self.gal_noise