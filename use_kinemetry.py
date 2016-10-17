import numpy as np
import matplotlib.pyplot as plt
from checkcomp import checkcomp
import os
cc=checkcomp()

def use_kinemetry(gal):
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	f = '%s/%s/kinemetry.txt' % (out_dir, gal)
	if os.path.exists(f):
		rad, pa, er_pa, q, er_q, k1, erk1, k51, erk51 = np.loadtxt(f, unpack=True, skiprows=1)

		f, ax = plt.subplots()
		a=ax.plot(rad,pa,label='PA')

		ax2=ax.twinx()
		b=ax2.plot(rad,k1,'r', label='k1')
		ax2.spines['right'].set_color('red')
		ax2.yaxis.label.set_color('red')
		ax2.tick_params(axis='y', colors='red')


		ax2.spines['left'].set_color('blue')
		ax.yaxis.label.set_color('blue')
		ax.tick_params(axis='y', colors='blue')

		ax.ylabel('PA','b')
		ax2.ylabel('k1','r')

		# lns = a+b
		# labs = [l.get_label() for l in lns]
		# ax.legend(lns, labs, loc=0)

		ax.set_title('KINEMETRY output')



		f.savefig('%s/%s/kinemetry.png'%(out_dir,gal))
	else:
		print 'There is no KINEMETRY file for %s.' %(gal)

##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	use_kinemetry('ic1459')
