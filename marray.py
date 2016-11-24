import numpy as np
# Class to allow easy manipulation of array and mask
# taken from http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
class marray(np.ndarray):
	def __new__(cls, input_array):
		# Input array is an already formed ndarray instance
		# We first cast to be our class type
		obj = np.asarray(input_array).view(cls)
		# add the new attribute to the created instance
		obj.mask = np.zeros(len(obj))
		# Finally, we must return the newly created object:
		return obj

	def __new__(cls, input_array, mask):
		# Input array is an already formed ndarray instance
		# We first cast to be our class type
		obj = np.asarray(input_array).view(cls)
		# add the new attribute to the created instance
		obj.mask = mask
		# Finally, we must return the newly created object:
		return obj		

#	def __array_finalize__(self, obj):
#		# see InfoArray.__array_finalize__ for comments
#		if obj is None: return
#		self.uncert = getattr(obj, 'uncert', None)
	#pass
	def __array_wrap__(self,obj):
		a=np.where(obj.mask)
		b=[i for i in a-1 if i not in a]
		c=[i for i in a+1 if i not in a]
		for i in range(len(b)):
			if b[i] < 0:
				obj[:c[i]]=obj[c[i]]
			elif c[i] > len(obj):
				obj[b[i]+1:]=obj[b[i]]
			else:
				obj[b[i]+1:c[i]]=np.interp(range(b[i]+1,c[i]), [b[i],c[i]], 
					[obj[b[i]],obj[c[i]]])
		return np.ndarray.__array_wrap__(self,obj)
