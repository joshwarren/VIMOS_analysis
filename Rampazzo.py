## warrenj 20170209 Routines needed for rebinning into a slit-like results. 
## Also include testing of binning by modelling the galaxy as 2D gaussian.
import numpy as np
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
from tools import get_slit, get_slit2, slit, funccontains

class rampazzo(object):
	def __init__(self, galaxy, slit_h=4.5*60, slit_w=2, slit_pa=30, method='aperture',
		r1=0.0, r2=1.0, debug=False):
		self.galaxy = galaxy
		self.slit_w = slit_w
		self.slit_h = slit_h
		self.slit_pa = slit_pa
		self.debug = debug
		self.method = method
		self.r1 = r1
		self.r2 = r2
## ----------============== Load galaxy info ================---------
		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
		# different data types need to be read separetly
		z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
			unpack=True, skiprows=1, usecols=(1,2,3,4,5))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		self.i_gal = np.where(galaxy_gals==galaxy)[0][0]
		self.vel = vel_gals[self.i_gal]
		self.sig = sig_gals[self.i_gal]
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

		self.x_cent = x_gals[self.i_gal] * self.f[0].header['CDELT1'] 
		self.y_cent = y_gals[self.i_gal] * self.f[0].header['CDELT2']

		self.cube = self.f[0].data
		self.cube[self.f[3].data==1] = 0
		self.noise_cube = self.f[1].data
		self.noise_cube[self.f[3].data==1] = 0.000000001
		self.cube[~np.isfinite(self.noise_cube)] = 0
		self.noise_cube[~np.isfinite(self.noise_cube)] = 0.000000001
## ----------============ Find data within slit =============---------
		self.x = (np.arange(self.f[0].header['NAXIS1']) * self.f[0].header['CDELT1']
			).repeat(self.f[0].header['NAXIS2'])
		self.y = np.tile(np.arange(self.f[0].header['NAXIS2']) * 
			self.f[0].header['CDELT2'], self.f[0].header['NAXIS1'])

		self.slit_res = 0.82 # "/px


	def get_spec(self):
		exec('self.'+self.method+'()')
		return self.spec, self.noise


	# Sets the slit to be the box self.r1 < r < self.r2 and finds the contribution to that 
	# slit
	def gradient(self):
		slit_x_cent = ((self.r2 - self.r1)/2 + self.r1) * np.sin(np.radians(self.slit_pa)
			) + self.x_cent
		slit_y_cent = ((self.r2 - self.r1)/2 + self.r1) * np.cos(np.radians(self.slit_pa)
			)+ self.y_cent
		slit_corners = get_slit2(self.slit_h, 2, self.slit_pa, slit_x_cent, slit_y_cent)
		if self.debug:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).contains.astype(
				float)
		else:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).fraction

		frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
		frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(
			self.f[0].header['NAXIS2']),
			np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].
			header['NAXIS1'])] = frac

		# Multiply by r for integral (eq (3))
		x = np.outer((np.arange(self.f[0].header['NAXIS1'])* f[0].header['CDELT1'] - 
			self.x_cent), np.ones(self.f[0].header['NAXIS2']))
		y = np.outer(np.ones(self.f[0].header['NAXIS1']), 
			(np.arange(self.f[0].header['NAXIS2'])* f[0].header['CDELT2'] - self.y_cent) )
		frac_image *= np.sqrt(x**2 + y**2)

		cube = np.einsum('ijk,jk->ijk', self.cube, frac_image)
		noise = np.sqrt(np.einsum('ijk,jk->ijk', self.noise_cube**2, 
				frac_image**2))

		half_int = np.trapz(cube, dx=self.f[0].header['CDELT2'],axis=2)
		self.spec = np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1)

		half_int = np.trapz(noise_cube**2, dx=self.f[0].header['CDELT2'],axis=2)
		# sqrt occurs at end of routine
		self.noise = np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1)


		# Other half of the slit
		slit_x_cent = -((self.r2 - self.r1)/2 + self.r1) * np.sin(np.radians(self.slit_pa)
			) + self.x_cent
		slit_y_cent = -((self.r2 - self.r1)/2 + self.r1) * np.cos(np.radians(self.slit_pa)
			) + self.y_cent
		slit_corners = get_slit2(self.slit_h, 2, self.slit_pa, slit_x_cent, slit_y_cent)
		if self.debug:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).contains.astype(
				float)
		else:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).fraction			
		frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
		frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(
			self.f[0].header['NAXIS2']),
			np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].
			header['NAXIS1'])] = frac

		# Multiply by r for integral (eq (3))
		x = np.outer((np.arange(self.f[0].header['NAXIS1'])* f[0].header['CDELT1'] - 
			self.x_cent), np.ones(self.f[0].header['NAXIS2']))
		y = np.outer(np.ones(self.f[0].header['NAXIS1']), 
			(np.arange(self.f[0].header['NAXIS2'])* f[0].header['CDELT2'] - self.y_cent))
		frac_image *= np.sqrt(x**2 + y**2)
			
		cube = np.einsum('ijk,jk->ijk', self.cube, frac_image)
		noise = np.sqrt(np.einsum('ijk,jk->ijk', self.noise_cube**2, 
				frac_image**2))

		half_int = np.trapz(cube, dx=self.f[0].header['CDELT2'],axis=2)
		self.spec += np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1)

		half_int = np.trapz(noise_cube**2, dx=self.f[0].header['CDELT2'],axis=2)
		self.noise += np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1)
		self.noise = np.sqrt(self.noise)


	# Finds the contribution to the entire slit and then sets any pixel outside 
	# of self.r1 < r < self.r2 to zero.
	def gradient2(self):
		slit_corners = get_slit(self.galaxy, self.slit_h, 2, self.slit_pa)
		if self.debug:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).contains.astype(
				float)
		else:
			frac = funccontains(slit, (slit_corners), x=self.x, y=self.y).fraction

		frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
		frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(
			self.f[0].header['NAXIS2']),
			np.tile(np.arange(self.f[0].header['NAXIS2']),self.f[0].
			header['NAXIS1'])] = frac

		eff_r = frac_image * self.f[0].header['CDELT1'] * self.f[0].header['CDELT2'] / \
			self.slit_w

		# Multiply by r for integral (eq (3))
		x = np.outer((np.arange(self.f[0].header['NAXIS1'])* self.f[0].header['CDELT1'] - 
			self.x_cent), np.ones(self.f[0].header['NAXIS2']))
		y = np.outer(np.ones(self.f[0].header['NAXIS1']), 
			(np.arange(self.f[0].header['NAXIS2'])* self.f[0].header['CDELT2'] - 
			self.y_cent))
		r = np.sqrt(x**2 + y**2)
		r[(r > self.r2) + (r < self.r1)] = 0
		set_r = np.array(r)
		set_r[set_r != 0] = 1
		frac_image *= set_r

		cube = np.einsum('ijk,jk->ijk', self.cube, frac_image*eff_r)
		noise = np.sqrt(np.einsum('ijk,jk->ijk', self.noise_cube**2, 
				frac_image**2*eff_r))

		self.spec = np.sum(cube, axis=(1,2))
		self.noise = np.sqrt(np.sum(noise**2, axis=(1,2)))

		# half_int = np.trapz(cube, dx=self.f[0].header['CDELT2'],axis=2)
		# self.spec = np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1)

		# half_int = np.trapz(noise_cube**2, dx=self.f[0].header['CDELT2'],axis=2)
		# self.noise = np.sqrt(np.trapz(half_int, dx=self.f[0].header['CDELT1'],axis=1))

		cube = np.einsum('ijk,jk->ijk', self.cube, frac_image*eff_r*r)
		self.r = np.mean(np.sum(cube, axis=(1,2))/self.spec)

		self.spec /= abs(self.r2-self.r1)
		self.noise /= abs(self.r2-self.r1)




	def aperture(self):
		slit_r = np.arange(-self.slit_h/2 + self.slit_res/2, self.slit_h/2 - 
			self.slit_res/2, self.slit_res)
		n_spaxels = len(slit_r)
		slit_x_cent = slit_r * np.sin(np.radians(self.slit_pa)) + self.x_cent
		slit_y_cent = slit_r * np.cos(np.radians(self.slit_pa)) + self.y_cent

		
		slit_pixels = get_slit2(self.slit_res, 2, self.slit_pa, slit_x_cent, slit_y_cent)

		self.spec = np.zeros((n_spaxels, self.f[0].header['NAXIS3']))
		self.noise = np.zeros((n_spaxels, self.f[0].header['NAXIS3']))
		self.r = np.zeros((n_spaxels, self.f[0].header['NAXIS3']))

		for i in xrange(n_spaxels):
			if self.debug:
				frac = funccontains(slit, (slit_pixels[0][:,i], slit_pixels[1][:,i]), 
					x=self.x, y=self.y).contains.astype(float)
			else:
				frac = funccontains(slit, (slit_pixels[0][:,i], slit_pixels[1][:,i]), 
					x=self.x, y=self.y).fraction
			
			frac_image = np.zeros((self.f[0].header['NAXIS1'], self.f[0].header['NAXIS2']))
			frac_image[np.arange(self.f[0].header['NAXIS1']).repeat(
				self.f[0].header['NAXIS2']),
				np.tile(np.arange(self.f[0].header['NAXIS2']),
				self.f[0].header['NAXIS1'])] = frac

			x = np.linspace(self.x_cent - slit_r[i] - self.slit_res, 
				self.x_cent + slit_r[i] + self.slit_res, 1000).repeat(1000)
			y = np.tile(np.linspace(self.y_cent - slit_r[i] - self.slit_res, 
				self.y_cent + slit_r[i] + self.slit_res, 1000), 1000)

			in_annulus = in_ellipse(x, y, self.x_cent, 
				self.y_cent, slit_r[i] + self.slit_res/2, 
				self.ellip, self.slit_pa) ^ in_ellipse(x, y, self.x_cent, self.y_cent, 
				slit_r[i] - self.slit_res/2, self.ellip, self.slit_pa)
			
			# Circular aperture (a=b=self.r2 <=> e=0, pa=0)
			in_aperture = in_ellipse(x, y, self.x_cent, self.y_cent, self.r2, 0, 0)

			# import matplotlib.pyplot as plt 
			# f,ax = plt.subplots()
			# ax.scatter(x[in_aperture],y[in_aperture])
			# ax.scatter(x[~in_aperture],y[~in_aperture], color='r')
			# ax.scatter(x[in_annulus], y[in_annulus], color='g')
			# ax.scatter(x[in_annulus*in_aperture], y[in_annulus*in_aperture], color='y')			
			# ax.set_aspect('equal')
			# plt.show()

			# Fraction of area of elliptical annulus that is within the aperture
			if np.sum(in_annulus) != 0:
				area_frac = float(np.sum(in_annulus*in_aperture))/np.sum(in_annulus)
			else: area_frac = 0

			scaling_factor = area_frac * np.pi * np.sqrt(1 - self.ellip) * abs(
				(slit_r[i] + self.slit_res/2)**2 - (slit_r[i] - self.slit_res/2)**2)/(
				self.slit_res * self.slit_w) / 2

			self.spec[i,:] = np.einsum('ijk,jk->i', self.cube, frac_image) * \
				scaling_factor
			self.noise[i,:] = np.einsum('ijk,jk->i', (self.noise_cube * 
				scaling_factor)**2, frac_image**2)

			self.r[i,:] = np.einsum('ijk,jk->i', self.cube, frac_image) * \
				scaling_factor * abs(slit_r[i])



		# Taking the mean is not in the paper, but otherwise <r> is wavelength dependent...
		self.r = np.mean(np.sum(self.r, axis=0) / np.sum(self.spec, axis=0))
		self.spec = np.sum(self.spec, axis=0) / np.pi * self.r2**2
		self.noise = np.sqrt(np.sum(self.noise, axis=0)) / (np.pi * self.r2**2)






# Much faster and simpler than using funccontains for aperture method
def in_ellipse(xp, yp, x0, y0, a, e, pa):

	cosa=np.cos(pa)
	sina=np.sin(pa)
	b = a * np.sqrt(1 - e**2)

	x_prime = cosa * (xp - x0) + sina * (yp - y0)
	y_prime = sina * (xp - x0) - cosa * (yp - y0)
	ellipse = (x_prime/a)**2 + (y_prime/b)**2

	return ellipse <= 1


def fake_galaxy():
	#cube = np.zeros((5,1000,1000))
	x = np.exp(-(np.arange(1000) - 500.0)**2 / (2 * 200**2))
	y = np.exp(-(np.arange(1000) - 500.0)**2 / (2 * 150**2))
	cube = np.outer(x,y)
	import matplotlib.pyplot as plt 
	# plt.imshow(cube)
	# plt.show()

	cube = np.array([cube,cube,cube,cube,cube])
	fits_cube = np.zeros((5,40,40))
	res = 1000/40
	for i in xrange(40):
		for j in xrange(40):
			fits_cube[:, i, j] = np.sum(cube[:, i*res:(i+1)*res, j*res:(j+1)*res], 
				axis=(1,2))

	data=rampazzo('ngc3557', slit_pa=0)
	data.f[0].data=fits_cube
	data.f[0].header['NAXIS3'] = 5

	fits_s, _ = data.get_spec()



	slit_r = np.arange(-(20*0.67) + data.slit_res/2, (20*0.67) - data.slit_res/2, 
		data.slit_res)
	n_spaxels = len(slit_r)
	slit_x_cent = slit_r * np.sin(np.radians(0)) + (20*0.67) # half the size of IFU FoV
	slit_y_cent = slit_r * np.cos(np.radians(0)) + (20*0.67)








if __name__=='__main__':
	#rampazzo('ic1459', method='gradient1',debug=True)
	fake_galaxy()


