## ==================================================================
## 		Find which computer
## ==================================================================
## warrenj 20160725 A routine to check which computer a given program 
## is currently being run on.
import platform

class checkcomp(oject):
	def __init__(self):
		self.device = getDevice()
		self.base_dir = getBase(self)

	def getDevice():
	    uname = platform.uname()
    	if uname contains 'Hal': 
	    	device = 'home'
    	elif uname contain 'glamring':
	    	device = 'glamdring'
	    else:
	    	device = -1

        return device

    def getBase(self):
    	if comp == 'home':
    		base_dir = '/cygdrive/x'
	    elif comp == 'glamdring':
    		raise ValueError('This may not setup to run on Glamdring')
    		base_dir = '~'
    	elif comp == 'uni':
	    	base_dir = ''
	    return base_dir





##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
	checkcomp()