## ==================================================================
## Unsharp mask HST images
## ==================================================================
## warrenj 20151020 Routine to apply an unsharp mask on various 
## scales to HST images.


from prefig import Prefig 
Prefig()
import matplotlib.pyplot as plt # used for plotting
import glob # for searching for files
from astropy.io import fits
import numpy as np # for array handling
from scipy import ndimage

#-----------------------------------------------------------------------------

def unsharpmask(galaxy):


	sigma = 10
	weight = 0.6
	clip = 3

	fig, ax = plt.subplots(1,2)
   
	img_file = glob.glob("/Data/hst/%s/*.fits" %(galaxy))

	for file in img_file:
		f = fits.open(file)
		img = f[1].data 
		hst_header = f[1].header
  
  		try:
			pa = hst_header['OORIENTA']
			img = ndimage.rotate(img, pa, reshape=False)
		except KeyError:
			pa = 0
				
		fig, ax = plt.subplots(1,2)
		m = np.nanmean(img)
		s = np.nanstd(img)
		ax[0].imshow(np.rot90(img.clip(m - clip*s,m + clip*s)), cmap='gist_yarg')
		ax[0].set_title('Original image for %s' % (galaxy.upper()))
	
		# fig.text(0.02,0.9, "Galaxy: " + galaxy.upper())

		blurred = ndimage.filters.gaussian_filter(img,sigma)

		unsharp = img - blurred*weight

		ax[1].imshow(np.rot90(unsharp.clip(m - clip*s,m + clip*s)), cmap='gist_yarg')
		ax[1].set_title('Unsharped image for %s' % (galaxy.upper()))

		fig.savefig(file.replace('fits', 'png'))


if __name__ == '__main__':
	for g in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399', 'ngc3557']:
		unsharpmask(g)


