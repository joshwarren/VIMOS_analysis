## warrenj 20170209 Routines needed for rebinning into a slit-like results. 
## Also include testing of binning by modelling the galaxy as 2D gaussian.
import numpy as np
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
from tools import get_slit, get_slit2, slit, funccontains

class rampazzo_binning(object):
	def __init__(self, galaxy, slit_h=4.5*60, slit_w=2, slit_pa=30, method='aperture',
		debug=False):
		print galaxy
		self.galaxy = galaxy
		self.slit_w = slit_w
		self.slit_h = slit_h
		self.slit_pa = slit_pa
		self.debug = debug
## ----------============== Load galaxy info ================---------
		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
		# different data types need to be read separetly
		z_gals, vel_gals, sig_gals, self.x_gals, self.y_gals = np.loadtxt(data_file, 
			unpack=True, skiprows=1, usecols=(1,2,3,4,5))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		self.i_gal = np.where(galaxy_gals==galaxy)[0][0]
		vel = vel_gals[self.i_gal]
		sig = sig_gals[self.i_gal]
		self.z = z_gals[self.i_gal]


		data_file = "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)
		# different data types need to be read separetly
		ellip = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(2,))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		i_gal2 = np.where(galaxy_gals==galaxy)[0][0]
		self.ellip = ellip[i_gal2]






		self.f = fits.open('%s/Data/vimos/cubes/%s.cube.combined.corr.fits' % (
			cc.base_dir, galaxy))
		self.lam = np.arange(self.f[0].header['NAXIS3'])*self.f[0].header['CDELT3'] + \
			self.f[0].header['CRVAL3']

		self.cube = self.f[0].data
		self.cube[self.f[3].data==1] = 0
		self.noise_cube = self.f[1].data
		self.noise_cube[self.f[3].data==1] = 0.000000001
		self.cube[~np.isfinite(self.noise_cube)] = 0
		self.noise_cube[~np.isfinite(self.noise_cube)] = 0.000000001
## ----------============ Find data within slit =============---------
		self.x = (np.arange(self.f[0].header['NAXIS1']) * self.f[0].header['CDELT1']).repeat(
			self.f[0].header['NAXIS2'])
		self.y = np.tile(np.arange(self.f[0].header['NAXIS2']) * self.f[0].header['CDELT2'],
			self.f[0].header['NAXIS1'])

		self.slit_res = 0.82 # "/px


		if method == 'aperture':
			self.aperture()
		elif method == 'gradient1':
			self.gradient1()









	def gradient1(self):
		slit_corners = get_slit(self.galaxy, self.slit_h, 2, self.slit_pa)
		frac = funccontains(slit, (slit_corners), x=self.x, y=self.y, 
			fraction=not self.debug).astype(float)

		frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
		frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(self.f[0].header['NAXIS2']),
			np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].header['NAXIS1'])] = frac
			
		self.cube = np.einsum('ijk,jk->ijk', self.cube, frac_image)
		self.noise = np.sqrt(np.einsum('ijk,jk->ijk', self.noise_cube**2, 
				frac_image**2))

		spec = np.trapz(self.cube, dx=self.f[0].header['CDELT2'],axis=2)
		spec = np.trapz(spec, dx=self.f[0].header['CDELT1'],axis=1)




	def aperture(self):
		slit_r = np.arange(-self.slit_h/2 + self.slit_res/2, self.slit_h/2 - self.slit_res/2, self.slit_res)
		n_spaxels = len(slit_r)
		slit_x_cent = slit_r * np.sin(np.radians(self.slit_pa)) + (self.x_gals[self.i_gal] * 
			f[0].header['CDELT1'])
		slit_y_cent = slit_r * np.cos(np.radians(self.slit_pa)) + (self.y_gals[self.i_gal] * 
			f[0].header['CDELT2'])

		
		slit_pixels = get_slit2(self.slit_res, 2, self.slit_pa, slit_x_cent, slit_y_cent)

		self.spec = np.zeros((n_spaxels, self.f[0].header['NAXIS3']))
		self.noise = np.zeros((n_spaxels, self.f[0].header['NAXIS3']))

		for i in xrange(n_spaxels):
			frac = funccontains(slit, (slit_pixels[0][:,i], slit_pixels[1][:,i]), x=self.x, y=self.y, 
				fraction=not self.debug).astype(float)
			
			frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
			frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(self.f[0].header['NAXIS2']),
				np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].header['NAXIS1'])] = frac
			
			scaling_factor = np.pi * np.sqrt(1 - self.ellip) * ((slit_r[i] + self.slit_res/2)**2 - 
				(slit_r[i] - self.slit_res/2)**2)/(self.slit_res * self.slit_w) 

			self.spec[i,:] = np.einsum('ijk,jk->i', self.cube, frac_image) * scaling_factor
			self.noise[i,:] = np.sqrt(np.einsum('ijk,jk->i', self.noise_cube**2, 
				frac_image**2)) * scaling_factor













		


if __name__=='__main__':
	rampazzo_binning('ic1459', method='gradient1',debug=True)

