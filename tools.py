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
		## sampling with n**2 * 100 points, where n is the number of points supplied
		x_sample = np.linspace(min(self.x), max(self.x), np.sqrt(len(self.x))*10).repeat(
			np.sqrt(len(self.y))*10)
		y_sample = np.tile(np.linspace(min(self.y), max(self.y), np.sqrt(len(self.y))*10), 
			np.sqrt(len(self.x))*10)

		xdelt = np.subtract.outer(x_sample,self.x)
		ydelt = np.subtract.outer(y_sample,self.y)
		sample_ownership = np.argmin(xdelt**2+ydelt**2, axis=1)
		
		contained = funccontains(func, *self.args, x=x_sample, y=y_sample, **self.kwargs
			).contains

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
    if (ceil==True & floor==True): raise "Can't CEIL and FLOOR, just one or t'other"

    
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
	if (ceil==True & floor==True): raise "Can't CEIL and FLOOR, just one or t'other"
	if ~np.isfinite(a): raise "Comparison value must be finite"

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
		raise 'm must be 2D for the routine to work'

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
