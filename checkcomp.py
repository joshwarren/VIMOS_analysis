## ==================================================================
## 		Find which computer
## ==================================================================
## warrenj 20160725 A routine to check which computer a given program 
## is currently being run on.
##
## Attributes:
## uname			device 			base_dir
## Hal				home			/cygdrive/x
## glamdring 		glamdring 		~
## asosx146 		uni 			

import platform

class checkcomp(object):

	def __init__(self):
		self.device = self.getDevice()
		self.base_dir = self.getBase_dir()

		
	def getDevice(self):
		uname = platform.uname()
		if 'Hal' in uname and 'CYGWIN_NT-10.0' in uname: 
			d = 'cygwin_home'
		elif 'HAL' in uname and 'Linux' in uname:
			d = 'bash_home'
		elif 'glamring' in uname:
			d = 'glamdring'
		elif 'asosx146' in uname:
			d = 'uni'
		else:
			d = -1 
		
		return d


	def getBase_dir(self):
		comp = self.device
		if comp == 'cygwin_home':
			b_dir = '/cygdrive/x'
		elif comp == 'bash_home':
			b_dir = '/mnt/x'
		elif comp == 'glamdring':
			raise ValueError('This may not setup to run on Glamdring')
			b_dir = '~'
		elif comp == 'uni':
			b_dir = ''
		else:
			b_dir = ''

		return b_dir