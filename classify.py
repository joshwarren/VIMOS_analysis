## ==================================================================
## 		Calculate classifications 
## ==================================================================
## warrenj 20161114 Pull together several calculations from the 
## kinematics analysis to classify each galaxy

import numpy as np # for array handling
from checkcomp import checkcomp
cc = checkcomp()
import re # for regex expressions
from rolling_stats import rollmed

def get_R_e(galaxy):
	if cc.device != 'glamdring':
		galaxiesFile_Re =  "%s/Data/galaxies_R_e.txt" % (cc.base_dir)
	else:
		galaxiesFile_Re =  "%s/galaxies_R_e.txt" % (cc.base_dir)
	galaxy_gals = np.loadtxt(galaxiesFile_Re, dtype=str, usecols=(0,), skiprows=1,
		unpack=True)
	i_gal = np.where(galaxy_gals==galaxy)[0]

	log_R_e_RC3_gals, R_e_2MASS_gals = np.loadtxt(galaxiesFile_Re, unpack=True, 
		skiprows=1, usecols=(1,2))
	R_e_RC3 = 6*10**log_R_e_RC3_gals[i_gal]/2 # convert to arcsec
	R_e_2MASS = R_e_2MASS_gals[i_gal]

	return np.nanmean([R_e_RC3,R_e_2MASS])


def classify(galaxy, opt='kin'):
	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)
	classify_file = "%s/galaxies_classify.txt" % (analysis_dir)

	try:
		galaxy_gals, RR, NF, M2, KT, KDC, S2, LV, group = my_loadtxt(classify_file, 
			skiprows=1, dtype=str)

		i_gal = np.where(galaxy_gals==galaxy)[0]
		if len(i_gal) == 0:
			galaxy_gals = np.append(galaxy_gals, [galaxy])
			print galaxy_gals
			i_gal = np.where(np.array(galaxy_gals)==galaxy)[0]
			RR = list(RR) + ['-']
			NF = list(NF) + ['-']
			M2 = list(M2) + ['-']
			KT = list(KT) + ['-']
			KDC = list(KDC) + ['-']
			S2 = list(S2) + ['-']
			LV = list(LV) + ['-']
			group = list(group) + ['-']
		i_gal = i_gal[0]

	except:
		galaxy_gals = np.array([galaxy])
		RR = np.array(['-'], dtype='S3')
		NF = np.array(['-'], dtype='S3')
		M2 = np.array(['-'], dtype='S3')
		KT = np.array(['-'], dtype='S3')
		KDC = np.array(['-'], dtype='S5')
		S2 = np.array(['-'], dtype='S3')
		LV = np.array(['-'], dtype='S3')
		group = np.array(['-'])
		i_gal = np.array([0])

	R_e = get_R_e(galaxy)
# ------------================= RR/NRR ===================----------
	file = '%s/%s/%s/kinemetry/kinemetry_stellar_vel.txt' % (analysis_dir,galaxy,opt)
	rad, pa, pa_err, k1, k51 = np.loadtxt(file, usecols=(0,1,2,5,7), skiprows=1, 
		unpack=True)
	
	pa = rollmed(pa, 5)
	k1 = rollmed(k1, 5)


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
	possible_no_features = True
# Low Velocity 
	if np.mean(k1) < 5: 
		LV[i_gal] = 'LV'
	else:
		LV[i_gal] = '-'

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
				possible_no_features = False
			else: KT[i_gal] = '-'
		elif np.sum(change) < (1-smooth)*len(pa):
			twist = abs(np.arcsin(np.sin(max(pa[~change])-min(pa[~change]))))
			if twist > 10: 
				KT[i_gal] = 'KT'
				possible_no_features = False
			else: KT[i_gal] = '-'
		else: KT[i_gal] = '-'


# Double maximum (M2)
		k1_change = np.append(k1[rad>1.0][1:-1]-k1[rad>1.0][0:-2], 0) > 0
		k1_change = str(k1_change.astype(int))
		# at least 3 True's followed by at least 3 False' followed by at least 3 True's
		if re.search('1 1 1 ( 1)+ 0 0( 0)+ 1 1( 1)+', k1_change) is not None:
			M2[i_gal] = '2M'
			possible_no_features = False
		else: M2[i_gal] = '-'


# Kinematically distinct/decoupled core
	# CRC not included yet
		sharp = np.abs(difference) > 30
		kdc_location = np.logical_and(sharp, k1 < 0.15*max(k1))
		if any(kdc_location) and any(difference[0:int(np.median(
			np.where(kdc_location)[0]))]<3):
			KDC[i_gal] = str(round(np.median(rad[kdc_location]), 3))
			possible_no_features = False
		else: KDC[i_gal] = '-'

	
# 2-sigma galaxies - none seen in our sample
		S2[i_gal] = '-'

# No Features
		if possible_no_features:
			contant_fit = np.poly1d(np.polyfit(rad[20:], pa[20:], 0, 
				w=1/pa_err[20:]))(rad[20:])
			if np.sqrt(np.sum((contant_fit-pa[20:])**2)) < 3 * \
				np.sqrt(np.sum(pa_err[20:]**2)):
				NF[i_gal] = 'NF'
			else:
				NF[i_gal] = '-'
		else:
			NF[i_gal] = '-'

	if RR[i_gal] == 'NRR' and LV[i_gal] == 'LV':
		group[i_gal] = 'a'
	elif RR[i_gal] == 'NRR' and NF[i_gal] == 'NF':
		group[i_gal] = 'b'
	elif KDC[i_gal] != '-':
		group[i_gal] = 'c'
	elif S2[i_gal] != '-':
		group[i_gal] = 'd'
	elif (RR[i_gal] == 'RR' and NF[i_gal] == 'NF') or (
		RR[i_gal] == 'RR' and M2[i_gal] == 'M2') or (
		RR[i_gal] == 'RR' and KT[i_gal] == 'KT'):
		group[i_gal] = 'e'
	else:
		group[i_gal] = 'f'	


# ------------============== Save outputs ================----------
	template3 = "{0:13}{1:7}{2:5}{3:5}{4:5}{5:13}{6:5}{7:5}{8:7}\n"	
	f3 = open(classify_file, 'wb')
	f3.write(template3.format('Galaxy', 'RR/NRR', 'NF', '2M', 'KT', 'KDC (arcsec)', 
		'2S', 'LV', 'group'))

	for i in range(len(galaxy_gals)):
		f3.write(template3.format(galaxy_gals[i], RR[i], NF[i], M2[i], KT[i], KDC[i],
			S2[i], LV[i], group[i]))

	f3.close()
##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc1399'

	classify(galaxy)