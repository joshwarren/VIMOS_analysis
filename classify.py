## ==================================================================
## 		Calculate classifications 
## ==================================================================
## warrenj 20161114 Pull together several calculations from the 
## kinematics analysis to classify each galaxy

from pickler import pickler
import numpy as np # for array handling
from checkcomp import checkcomp
cc = checkcomp()
import re # for regex expressions



def classify(galaxy):
	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)
	classify_file = "%s/galaxies_classify.txt" % (analysis_dir)
	galaxiesFile_Re =  "%s/galaxies_R_e.txt" % (analysis_dir)

	galaxy_gals, RR, NF, NR, KT, M2, KDC = np.loadtxt(classify_file, skiprows=1, 
		unpack=True, usecols=(0,1,2,3,4,5,6), dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0]

	log_R_e_RC3_gals, R_e_2MASS_gals = np.loadtxt(galaxiesFile_Re, unpack=True, 
		skiprows=1, usecols=(1,2))
	R_e_RC3 = 6*10**log_R_e_RC3_gals[i_gal]/2 # convert to arcsec
	R_e_2MASS = R_e_2MASS_gals[i_gal]

	R_e = np.nanmean([R_e_RC3,R_e_2MASS])
# ------------================= RR/NRR ===================----------
	file = '%s/%s/kinemetry_vel.txt' % (analysis_dir,galaxy)
	rad, pa, k1, k51 = np.loadtxt(file, usecols=(0,1,5,7), skiprows=1, unpack=True)
	rad *= 0.67 # Pix to arcsec


	# Finding smoothest pa by add or subtracting 360 deg	
	for j in range(1,len(pa)):
	    test = np.array([pa[j]-360, pa[j], pa[j]+360])
	    a = np.argmin(np.abs(test-pa[j-1]))
	    if a==0: 
	        pa[j:]-=360
	    elif a==2: 
	        pa[j:]+=360

	###### NEEDS TO BE CHANGED TO LUMINOSITY WEIGHTED MEAN
	k5k1 = np.average(k51/k1, weights=np.ones(len(k1)))
	if k5k1 < 0.04:
		RR[i_gal] = 'RR'
	else: RR[i_gal] = 'NRR'
# ------------=========== kinematic features =============----------
	i_rad = np.where(rad < R_e)[0][-1] # Largest radius below R_e
	feature = False

	# Non-rotators
	if k1[i_rad] < 5: 
		NR[i_gal] = 'NR'
		feature = True
	else:
		NR[i_gal] = '-'

	# kinematcs twist
	# smooth: are 80% of pionts increaded or decreased from their neighbour?
	difference = np.append(pa[:-1]-pa[1:], 0) 
	change = difference < 0
	smooth = 0.7 	# percentage of points to be moving the same direction to be 
					# considered smooth
	if np.sum(change) > smooth*len(pa):
		twist = abs(np.arcsin(np.sin(max(pa[change])-min(pa[change]))))
		if twist > 10: 
			KT[i_gal] = 'KT'
			feature= True
		else: KT[i_gal] = '-'
	elif np.sum(change) < (1-smooth)*len(pa):
		twist = abs(np.arcsin(np.sin(max(pa[~change])-min(pa[~change]))))
		if twist > 10: 
			KT[i_gal] = 'KT'
			feature = True
		else: KT[i_gal] = '-'
	else: KT[i_gal] = '-'

	# Kinematically distinct/decoupled core
	# CRC not included yet
	sharp = np.abs(difference) > 30
	kdc_location = np.logical_and(sharp, k1 < 0.15*max(k1))
	if any(kdc_location) and any(difference[0:np.median(np.where(kdc_location)[0])]<3):
		KDC[i_gal] = 'KDC'
		feature = True
	else: KDC[i_gal] = '-'

	# Double maximum (M2)
	k1_change = np.append(pa[1:-1]-pa[0:-2], 0) > 0
	k1_change = str(k1_change.astype(int))
	# at least 3 True's followed by at least 2 False' followed by at least 3 True's
	if re.search('1 1( 1)+ 0( 0)+ 1 1( 1)+', k1_change) is not None:
		M2[i_gal] = '2M'
		feature = True
	else: M2[i_gal] = '-'

	# No features
	if not feature: NF[i_gal] = 'NF' 
	else: NF[i_gal] = '-'


# ------------============== Save outputs ================----------
	template3 = "{0:13}{1:7}{2:5}{3:5}{4:5}{5:5}{6:5}\n"

	f3 = open(classify_file, 'wb')
	f3.write(template3.format('Galaxy', 'RR/NRR', 'NF', 'NR', 'KT', '2M', 'KDC'))

	for i in range(len(galaxy_gals)):
		f3.write(template3.format(galaxy_gals[i], RR[i], NF[i], NR[i], KT[i], 
			M2[i], KDC[i]))

	f3.close()


##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc3557'

	classify(galaxy)