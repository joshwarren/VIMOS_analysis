from find_galaxy import find_galaxy #part of mge package: fits photometry
import numpy as np


def find_centre(galaxy, discard=0, opt="kin", D=None, 
	instrument='vimos'):
	analysis_dir = "%s/Data/%s/analysis" % (cc.base_dir, instrument)
	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)
	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()

	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)

	d = np.loadtxt(galaxiesFile, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	if instrument == 'vimos':
		z_gals, vel_gals, sig_gals = d[1][1:].astype(float), \
			d[2][1:].astype(float), d[3][1:].astype(float)
		x_gals, y_gals = d[4][1:].astype(int), d[5][1:].astype(int)
		SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(6,
			len(d))}
		i_gal = np.where(galaxy_gals==galaxy)[0][0]
	elif instrument == 'muse':
		x_gals, y_gals = d[1][1:].astype(int), d[2][1:].astype(int)
		SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(3,
			len(d))}
		i_gal = np.where(galaxy_gals==galaxy)[0][0]

	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)

	x_gals[i_gal] = int(round(f.xmed, 0))
	y_gals[i_gal] = int(round(f.ymed, 0))

	if instrument == 'vimos':
		template = "{0:12}{1:11}{2:10}{3:15}{4:4}{5:4}" + ''.join(
			['{%i:%i}'%(i+6,len(t)+1) for i, t in 
			enumerate(SN_gals.keys())]) + "\n"

		SN_titles = list(SN_gals.keys())
		with open(galaxiesFile, 'w') as file:
			file.write(template.format("Galaxy", "z", "velocity", 
				"sigma", "x", "y", *(s for s in SN_titles)))
			for i in range(len(galaxy_gals)):
				file.write(template.format(galaxy_gals[i], 
					str(round(z_gals[i],7)), str(round(vel_gals[i],4)), 
					str(round(sig_gals[i],4)), str(int(x_gals[i])), 
					str(int(y_gals[i])), 
					*(str(round(SN_gals[s][i],2)) for s in SN_titles)))

	elif instrument == 'muse':
		template = "{0:12}{1:4}{2:4}" + ''.join(
			['{%i:%i}'%(i+3,len(t)+1) for i, t in 
			enumerate(SN_gals.keys())]) + "\n"

		SN_titles = list(SN_gals.keys())
		with open(galaxiesFile, 'w') as file:
			file.write(template.format("Galaxy",  "x", "y", 
				*(s for s in SN_titles)))
			for i in range(len(galaxy_gals)):
				file.write(template.format(galaxy_gals[i], 
					str(int(x_gals[i])), str(int(y_gals[i])), 
					*(str(round(SN_gals[s][i],2)) for s in SN_titles)))

	return D

