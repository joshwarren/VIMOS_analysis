import numpy as np
import cPickle as pickle
from checkcomp import checkcomp
cc = checkcomp()



galaxies = np.array(['ngc3557',
		'ic1459', 
		'ic1531', 
		'ic4296', 
		'ngc0612',
		'ngc1399',
		'ngc3100',
		'ngc7075', 
		'pks0718-34', 
		'eso443-g024'])

for galaxy in galaxies:
	bin=25
	wav_range='4200-'

	dir = '%s/Data/vimos' % (cc.base_dir)
	out_pickle = '%s/analysis/%s/results/%s/pickled' % (dir,galaxy,wav_range)
	pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	x = []
	y = []
	bin_num = []
	xBin = []
	yBin = []
	xBar = []
	yBar = []


	for i, bin in enumerate(D.bin):
		x.extend(bin.xspaxels)
		y.extend(bin.yspaxels)
		for a in range(bin.n_spaxels_in_bin):
			bin_num.append(i)
			xBin.append(bin.xBar)
			yBin.append(bin.yBar)

		xBar.append(bin.xBar)
		yBar.append(bin.yBar)

	x = np.array(x)
	y = np.array(y)
	xBin = np.array(xBin)
	yBin = np.array(yBin)
	bin_num = np.array(bin_num)

	xBar = np.array(xBar)
	yBar = np.array(yBar)

	o = np.argsort(x*36+y)

	x = x[o]
	y = y[o]
	bin_num = bin_num[o]
	xBin = xBin[o]
	yBin = yBin[o]

	temp = "{0:3}{1:3}{2:8}{3:9}{4:9}\n"
	temp2 = "{0:9}{1:9}\n"

	with open("%s/analysis/%s/voronoi_2d_binning_output_%s.txt" % (dir,galaxy,'pop'), 
		'w+') as f:
		f.write(temp.format('X"', 'Y"', 'BIN_NUM', 'XBIN', 'YBIN'))
		for i in range(len(xBin)):
			int(x[i])
			f.write(temp.format(str(int(x[i])), str(int(y[i])), str(int(bin_num[i])), 
				str(round(xBin[i],5)), str(round(yBin[i],5))))


	with open("%s/analysis/%s/voronoi_2d_binning_output2_%s.txt" % (dir,galaxy,'pop'), 
		'w+') as f:
		f.write(temp2.format('XBAR','YBAR'))
		for i in range(len(xBar)):
			f.write(temp2.format(str(round(xBar[i],5)), str(round(yBar[i],5)))) 