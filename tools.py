## warrenj 20161214 A selection of python tools 
import numpy as np

##############################################################################
## warrenj 20161214 A routine to find if a point or array of points in contained 
##	within (including on the boundary of) a given function.
## Parameters/keywords:
## func:				The function to evaluate. This should return an upper and a 
##						lower value for its output if within the x range, otherwise 
##						an nan, None or Error will be treated as outside of the 
##						range.
##						NB: if func returns a single value this will also be counted 
##						as outside of the boundary. If func returns more than 2 
##						results, this rountine will not work
## points 		None	An array of points: points[:,0] = x values, points[:,1] = y
##		or
## x			None	A list or array of the x values. NB y must also be supplied
##						and points must not be supplied.
## y 			None	As x above
## fraction		False	If True routine returns approx fraction of bin within func.
##############################################################################
class funccontains(object):
	def __init__(self, func, *args, **kwargs):
		# obvious *args syntax is a only avaliable for Python 3. In Python 2, arguments 
		# and keywords can both be set and use *args in this way. Therefore, default 
		# values of keywords are set below.
		self.points = kwargs.pop('points', None) 
		self.x = kwargs.pop('x', None) 
		self.y = kwargs.pop('y', None)
		self.func = func
		self.args = args
		self.kwargs = kwargs

		if self.points is None and self.x is None and self.y is None:
			raise ValueError('Points not supplied')
		elif (self.x is not None) ^ (self.y is not None):
			raise ValueError('Both x and y values must be supplied')
		elif self.x is not None and self.points is not None:
			raise ValueError('Please supply only x and y or points')
		
		if self.x is not None and self.y is not None:
			self.points = np.vstack([self.x,self.y]).T
		# else:
		# 	x = points[:,0]
		# 	y = points[:,1]

	@property
	def contains(self):
		result = []
		yf = self.func(self.points[:,0], *self.args)


		for p in self.points:
			yf = self.func(p[0], *self.args)
			yf = np.array(yf)[~np.isnan(yf)]
			if length(yf)%2 == 0:
				r = False
				yf = np.sort(yf)
				for i in range(len(yf)/2):
					r += p[1] <= np.nanmax(yf[i*2:i*2+2]) and \
						p[1] >= np.nanmin(yf[i*2:i*2+2])
				result.append(r)
			else: result.append(False)
		return np.array(result)

	@property
	def fraction(self):
		## sampling with n * 100 points, where n is the number of points supplied
		x_sample = np.linspace(min(self.x), max(self.x), 
			np.ceil(np.sqrt(len(self.x))*10)).repeat(np.ceil(np.sqrt(len(self.y))*10))
		y_sample = np.tile(np.linspace(min(self.y), max(self.y), 
			np.ceil(np.sqrt(len(self.y))*10)), int(np.ceil(np.sqrt(len(self.x))*10)))

		sample_ownership = np.argmin(np.subtract.outer(x_sample,self.x)**2 + 
			np.subtract.outer(y_sample,self.y)**2, axis=1)
		
		contained = funccontains(self.func, *self.args, x=x_sample, y=y_sample, 
			**self.kwargs).contains

		# Find fraction
		inside, counts_in = np.unique(sample_ownership[contained], return_counts=True)
		total, counts_tot = np.unique(sample_ownership, return_counts=True)

		frac = np.zeros(len(self.x))
		frac[inside] = counts_in
		frac /= counts_tot
		
		return frac
##############################################################################

##############################################################################
## warrenj 20161215 A routine to adapt len to return 1 if the argument is a float 
##	or int, similar to IDL's n_elements command.
def length(f):
	try:
		r = len(f)
	except:
		if f is None:
			r = 0
		else:
			r = 1
	return r
##############################################################################

##############################################################################
## Ryan Houghton nearest.py
def roundToNearest(x, dx=None, ceil=False, floor=False):
	"""
	Purpose: to round a number to the nearest multiple of the dx you specify
			 Particularly useful for plots 
			 e.g. nearest(234,50) = 250 but nearest(234,50,/floor) = 200
			 if dx==None, round to nearest Order of Mag
	"""

	# sanity check
	if (ceil==True & floor==True): raise ValueError("Can't CEIL and FLOOR, just one or"+
		" t'other")

	
	if dx==None:
		dx = 10.0**np.floor(np.log10(np.fabs(x)))
		if (~np.isfinite(np.log10(dx))): dx=10.0 # was rounding zero value

	near = float(x) / float(dx)
	if ceil:
		result = np.ceil(near)*dx
	elif floor:
		result = np.floor(near)*dx
	else:
		result = round(near)*dx

	return result
##############################################################################

##############################################################################
## warrenj 20170126 
## Purpose:	to return the index of an array which is closest to a given number.
def nearest(x, a, ceil=False, floor=False):
	if (ceil==True & floor==True): raise ValueError("Can't CEIL and FLOOR, just one"+
		" or t'other")
	if ~np.isfinite(a): raise ValueError("Comparison value must be finite")

	x = np.array(x)
	if ceil:
		x[x - a < 0] = np.inf
	elif floor:
		x[x - a > 0] = np.inf

	result = np.argmin(np.abs(x - a))

	return result
##############################################################################
## warrenj 20170201
"""
Purpose: 	to generate a list of the neighbouring cells in a given 2d array, m.
Example:
$ m = np.zeros((3,3))
$ for out in all_neighbours(m):
$     print out[0], out[1], out[2]
This will list the x,y coordinates of each point, followed by a list of the x 
coords of the neighbouring points and then the y coords.
"""
def all_neighbours(m):
	NEIGHBOURS = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
	s = m.shape
	if len(s) != 2:
		raise ValueError('m must be 2D for the routine to work')

	for i in xrange(len(m.flatten())):
		nx = []
		ny = []
		y, x = divmod(i, s[0])
		for u, v in NEIGHBOURS:
			ux = u + x
			vy = v + y
			if 0 <= ux < s[0] and 0 <= vy < s[1]:
				nx.append(ux)
				ny.append(vy)
		yield [x,y], nx,ny
##############################################################################

##############################################################################
## warrenj 20170205 
"""
Purpose:	Return f(x) where f represents a quadrilateral defined by its corners
			args = (x_points, y_points)
			NB: not sure how this would perform with a concave shape
"""
def slit(x, args):
	from tools import length as len
	x_points,y_points = args
	result = np.zeros((len(x),2))

	x_points = np.asarray(x_points).astype(float)
	y_points = np.asarray(y_points).astype(float)
	x = np.asarray(x).astype(float)

	# Bin is inside slit
	above = x_points < np.outer(x, np.ones(4))
	inside = (0 < np.sum(above, axis=1)) * (4 > np.sum(above, axis=1))

	# A Reference Index array
	ref = np.outer(np.ones(len(x)),np.arange(4))
	ref[above] = np.nan
	ref[~inside,:] = np.nan # If outside of the quadrilateral.


	# Lines which are interected by vertical lines at x
	# NB: use_lines is only transfered to int type later in order to retain the nan's for
	#	as long as possible
	use_lines = np.dstack([np.vstack([np.nanmax(ref, axis=1), np.nanmin(ref, axis=1)]).T,
		np.vstack([np.nanmax(ref, axis=1)+1, np.nanmin(ref,axis=1)-1]).T])
	use_lines[use_lines==4] = 0 # Cyclic boundary
	use_lines[use_lines==-1] = 3

	# Calc
	x0 = np.full((len(x),2),np.nan)
	x0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(
		x_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape(
		(np.sum(inside),2))
	x1 = np.full((len(x),2),np.nan)
	x1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(
		x_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape(
		(np.sum(inside),2))

	y0 = np.full((len(x),2),np.nan)
	y0[~np.isnan(use_lines[:,0,0]),:] = np.asarray(
		y_points[use_lines[:,:,0][~np.isnan(use_lines[:,0,:])].astype(int)]).reshape(
		(np.sum(inside),2))
	y1 = np.full((len(x),2),np.nan)
	y1[~np.isnan(use_lines[:,0,1]),:] = np.asarray(
		y_points[use_lines[:,:,1][~np.isnan(use_lines[:,1,:])].astype(int)]).reshape(
		(np.sum(inside),2))

	result = y0 + (y1 - y0)/(x1 - x0) * (np.outer(x, np.ones(2)) - x0)

	return result[:,0], result[:,1]
##############################################################################

##############################################################################
# warrenj 20170206 Routine to get the coordinates of a slit from hight (h), width(w) 
#	and postion angle (pa). NB: pa must be in degrees
def get_slit(galaxy, h, w, pa):
	from checkcomp import checkcomp
	cc = checkcomp()
	pa = np.radians(pa)
	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	res=0.67
	x_cent, y_cent = x_gals[i_gal]*res, y_gals[i_gal]*res
	return get_slit2(h, w, np.degrees(pa), x_cent, y_cent)

	

def get_slit2(h, w, pa, x_cent, y_cent):
	pa = np.radians(pa)
	# NB: RA (x) is increasing right to left
	x = np.array([x_cent + w/2, x_cent - w/2, x_cent - w/2, x_cent + w/2])
	y = np.array([y_cent + h/2, y_cent + h/2, y_cent - h/2, y_cent - h/2])
	# Apply pa (angle N to E or top to left)
	rotx = x_cent + np.cos(pa) * (x - x_cent) - np.sin(pa) * (y - y_cent)
	roty = y_cent + np.sin(pa) * (x - x_cent) + np.cos(pa) * (y - y_cent)

	return rotx,roty
##############################################################################


def fwhm(x, y, k=10, debug=False):
	from scipy.interpolate import splrep, sproot
	"""
	Determine full-with-half-maximum of a peaked set of points, x and y.

	Assumes that there is only one peak present in the datasset.  The function
	uses a spline interpolation of order k.

	Taken from https://stackoverflow.com/questions/10582795/finding-the-full-
	width-half-maximum-of-a-peak, written by user: jdg
	"""

	class MultiplePeaks(Exception): pass
	class NoPeaksFound(Exception): pass

	half_max = np.max(y)/2.0
	s = splrep(x, y - half_max)

	if debug:
		x2 = np.linspace(np.min(x), np.max(x), 200)
		from scipy.interpolation import splev
		y2 = splev(x2, s)
		import matplotlib.pyplot as plt
		plt.plot(x, y-half_max, 'o')
		plt.plot(x2, y2)
		plt.show()
	roots = sproot(s)

	if len(roots) > 2:
		raise MultiplePeaks("The dataset appears to have multiple peaks, and "
				"thus the FWHM can't be determined.")
	elif len(roots) < 2:
		raise NoPeaksFound("No proper peaks were found in the data set; likely "
				"the dataset is flat (e.g. all zeros).")
	else:
		return abs(roots[1] - roots[0])
##############################################################################
# Taken from answer by EOL from
# https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def weighted_avg_and_std(values, weights, axis=None):
	"""
	Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.
	"""
	average = np.average(values, axis=axis, weights=weights)
	# variance = np.average(values2**2, axis=axis, weights=weights)
	variance = np.average(np.subtract(values.T, average).T**2, 
		axis=axis, weights=weights)
	return (average, np.sqrt(variance))

def gaussian(x,amp=1,mean=0,sigma=1):
	return amp*np.exp(-(x-mean)**2/(2*sigma**2))

# Taken from answer by Jamie from 
# https://stackoverflow.com/questions/18517722/weighted-moving-average-in-python
def moving_weighted_average(x, y, step_size=.1, weights=None, interp=True):
	from numpy.lib.stride_tricks import as_strided

	# Default is gaussian
	if weights is None:
		x1 = np.arange(-5,6)
		weights = gaussian(x1,sigma=7.)

	# This ensures that all samples are within a bin
	number_of_bins = int(np.ceil(np.ptp(x) / step_size))
	bins = np.linspace(np.min(x), np.min(x) + step_size*number_of_bins,
					   num=number_of_bins+1)
	bins -= (bins[-1] - np.max(x)) / 2 # centering range
	bin_centers = bins[:-int(len(weights))] + step_size*len(weights)/2.
	# bin_centers = (bins[:-1] + bins[1:])/2.

	counts, _ = np.histogram(x, bins=bins)
	vals, _ = np.histogram(x, bins=bins, weights=y)
	bin_avgs = vals / counts
	n = len(bin_avgs)
	# windowed_bin_avgs = as_strided(bin_avgs, 
	# 	shape=(n, len(weights)), strides=bin_avgs.strides*2)
	windowed_bin_avgs = as_strided(bin_avgs,
		(n-len(weights)+1, len(weights)), bin_avgs.strides*2)
	

	weighted_average, weighted_std = weighted_avg_and_std(
		windowed_bin_avgs, axis=1, weights=weights)

	if interp:
		from scipy.interpolate import interp1d
		cut = 1
		inter = interp1d(bin_centers, weighted_average, bounds_error=False, 
			fill_value='extrapolate')
		weighted_average = inter(x)

		inter = interp1d(bin_centers, weighted_std, bounds_error=False, 
			fill_value='extrapolate')
		weighted_std = inter(x)

		bin_centers = np.array(x)

	return bin_centers, weighted_average, weighted_std



def myerrorbar(ax, x, y, xerr=None, yerr=None, zorder=0, color=None, 
	marker=None, size=None, colorbar=False, cmap=None):
	import matplotlib.pyplot as plt
	from matplotlib import cm
	
	sc = ax.scatter(x, y, marker=marker, s=size, c=color, 
		zorder=zorder, cmap=cmap)

	if colorbar:
		#create colorbar according to the scatter plot
		clb = plt.colorbar(sc)

	#create errorbar plot and return the outputs to a,b,c

	#convert time to a color tuple using the colormap used for scatter
	if color is not None:
		if not isinstance(color, str):
			a,b,c = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.',
				zorder=zorder-1)

			color = np.array(color).astype(float)
			if colorbar:
				time_color = clb.to_rgba(color)
			else:
				if isinstance(cmap, str):
					cmap = cm.get_cmap(cmap)
				elif cmap is None:
					cmap = cm.get_cmap()
				color = cmap((color - np.nanmin(color))/nanptp(color))

			# adjust the color of c[0], which is a LineCollection, 
			#	to the colormap
			if xerr is not None:
				m = np.isfinite(x) * np.isfinite(y) * np.isfinite(xerr)
				c[0].set_color(color[m])

				if yerr is not None:
					m = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr)
					c[1].set_color(color[m])

			elif yerr is not None:
				m = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr)
				c[0].set_color(color[m])


			
		else:
			a,b,c = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='',
				zorder=zorder-1, color=color)


def nanptp(v):
	m = np.isfinite(v)
	return np.nanmax(v[m]) - np.nanmin(v[m])


def e_binomial(p, n, str=False):
	p = float(p)/n
	e = np.sqrt((p * (1 - p))/n)
	if str:
		print '%.4f+/-%.4f' % (p, e)
	return e


##############################################################################
# Python implimentation of IDL code from 
# https://idlastro.gsfc.nasa.gov/ftp/pro/astro/calz_unred.pro
# requested by a Jane Rigby on Python users in Astronomy facebook page 

import numpy as np
import warnings

def calz_unred(wave, flux, ebv, R_V=4.05):

    w1 = np.where((wave >= 6300) * (wave <= 22000))[0]
    w2 = np.where((wave >= 912) * (wave < 6300))[0]

    x  = 10000.0/wave                      # Wavelength in inverse microns

    if len(w1) + len(w2) != len(wave):
        warnings.warn('Warning - some elements of wavelength vector outside valid'
            + ' domain')
        flux[(wave < 912) + (wave > 22000)] = 0

    klam = flux*0.0

    if len(w1) > 0:
        klam[w1] = 2.659*(-1.857 + 1.040*x[w1]) + R_V

    if len(w2) > 0:
        klam[w2] = 2.659 * np.poly1d([-2.156, 1.509, -0.198, 0.011][::-1])(x[w2])\
            + R_V

    funred = flux * 10.0**(0.4 * klam * ebv)
    return funred