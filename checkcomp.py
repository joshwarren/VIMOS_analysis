## ==================================================================
## 		Find which computer
## ==================================================================
## warrenj 20160725 A routine to check which computer a given program 
## is currently being run on.
## warrenj 20160916 Added override since glamdring cores do not know 
## they are part of glamdring.
##
## Attributes:
## uname			device 			base_dir
## Hal				home			/cygdrive/x
## glamdring 		glamdring 		~
## asosx146 		uni 			/
## asosx134 		uni 			/

import platform
import warnings
import os

class checkcomp(object):

	def __init__(self, override=None):
		if override is None:
			self.device = self.getDevice()
		else:
			self.device = override
		self.overriden = override is not None
		#self.base_dir = self.getBase_dir()
		self.home_dir = os.environ['HOME'] # self.getHome_dir() 
		self.remote = os.environ.get('DISPLAY') is None
		
	def getDevice(self):
		uname = platform.uname()
		if 'Hal' in uname and 'CYGWIN_NT-10.0' in uname: 
			d = 'cygwin_home'
		elif 'HAL' in uname and 'Linux' in uname:
			d = 'bash_home'
		elif 'glamdring' in uname:
			d = 'glamdring'
		elif 'asosx146' in uname or 'asosx134' in uname:
			d = 'uni'
		else:
			try:
				with open('/etc/clustername', 'r') as f:
					if 'glamdring' in f.read():
						d = 'glamdring'
					else:
						d=-1
			except IOError:
				d = -1 
		return d

	@property
	def base_dir(self):
		comp = self.device
		if comp == 'cygwin_home':
			#b_dir = '/cygdrive/x' # Moved from external to laptop
			b_dir = ''
		elif comp == 'bash_home':
			b_dir = '/mnt/x'
		elif comp == 'glamdring':
			if not self.overriden:
				warnings.warn('This may not setup to run on Glamdring')
			b_dir = self.home_dir
		elif comp == 'uni':
			b_dir = ''
		else:
			b_dir = ''

		return b_dir

	# def getHome_dir(self):
	# 	comp = self.device
	# 	if comp == 'cygwin_home':
	# 		h_dir = '/home/Home'
	# 	elif comp == 'bash_home':
	# 		h_dir = '/root'
	# 	elif comp == 'glamdring':
	# 		h_dir = '/users/warrenj'
	# 	elif comp == 'uni':
	# 		h_dir = '/home/warrenj'
	# 	else:
	# 		h_dir = ''

	# 	return h_dir
