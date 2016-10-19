import numpy as np
import matplotlib.pyplot as plt
from checkcomp import checkcomp
import os
cc=checkcomp()

def use_kinemetry(gal):
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	colors=['b','y','g']
	fig, ax = plt.subplots()
	pl_ob=[] # array for the plot objects
	missing=0 # number of missing plots (out of 3)

	# Plot all avaliable types
	for i, type in enumerate(['flux','vel','sigma']):
		f = '%s/%s/kinemetry_%s.txt' % (out_dir, gal, type)
		if os.path.exists(f):
			rad, pa, er_pa, q, er_q, k1, erk1, k51, erk51 = np.loadtxt(f, unpack=True, 
				skiprows=1)
			rad*=0.33 # Change to arcsec

			# Align pa[0] as closely as possible with flux pa[0] by +/- 360 deg
			if type == 'flux': pa0 = pa[0]
			else:
				a = np.argmin(np.abs(pa[0]-pa0+[-360,0,360]))
				if a == 0:
					pa -= 360
					print '-%s'%(type)
				elif a == 2:
					print '+%s'%(type)
					pa += 360

			# Optimizing PA
			# cut = np.arange(360,0,-10)
			# r = np.array([list(pa)]*len(cut))
			# for j, c in enumerate(cut):
			# 	r[j, np.where(r[j,:] > c)[0]] -=360
			# l = np.argmin(np.ptp(r,axis=1))
			# pa = r[l]

			# Finding smoothest pa by add or subtracting 360 deg	
			for j in range(1,len(pa)):
			    test = np.array([pa[j]-360, pa[j], pa[j]+360])
			    a = np.argmin(np.abs(test-pa[j-1]))
			    if a==0: 
			        pa[j:]-=360
			    elif a==2: 
			        pa[j:]+=360
			


			pl = ax.plot(rad,pa,colors[i],label='%s PA' % (type))
			pl_ob.append(pl)

			# Plot k1 from velocity fit.
			if type == 'vel':
				ax2=ax.twinx()
				b=ax2.plot(rad,k1,'r', label='k1')
				ax2.spines['right'].set_color('red')
				ax2.yaxis.label.set_color('red')
				ax2.tick_params(axis='y', colors='red')
				ax2.set_ylabel('k1',color='r', rotation=270, labelpad=10)

				ax2.spines['left'].set_color('blue')

		else:
			print 'There is no %s KINEMETRY file for %s.' %(type, gal)
			missing +=1

	# If not all missing
	if missing != 3:
		# Get redshift of galaxy from data_file
		data_file =  "%s/galaxies.txt" % (out_dir)
		# different data types need to be read separetly
		z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,))
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		i_gal = np.where(galaxy_gals==gal)[0][0]
		z = z_gals[i_gal]

		# Set secondary x axis in kpc.
		ax3=ax.twiny()
		c = 299792 #km/s
		#H = 67.8 #(km/s)/Mpc # From Planck
		H = 70.0 # value used by Bolonga group.
		r_kpc = np.radians(rad/(60.0*60.0)) * z*c/H *1000
		ax3.set_xlim(r_kpc[0],r_kpc[-1])
		ax3.set_xlabel('Radius (kpc)')

		#ax.yaxis.label.set_color('blue')
		#ax.tick_params(axis='y', colors='blue')
		ax.set_ylabel('PA')#,color='b')
		ax.set_xlabel('radius (arcsec)')

		# Legend
		try:
			lns = b
			for l in pl_ob: lns+=l
		except UnboundLocalError:
			if missing==1:lns=pl_ob[0]+pl_ob[1]
			else:lns=pl_ob[0]

		labs = [l.get_label() for l in lns]
		ax.legend(lns, labs, loc=0)

		# Moves title clear of upper x axis
		plt.subplots_adjust(top=0.85)
		ax.set_title('KINEMETRY output', y=1.12)

		fig.savefig('%s/%s/kinemetry.png'%(out_dir,gal))
	
##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	use_kinemetry('ic1459')
