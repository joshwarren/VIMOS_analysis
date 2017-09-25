import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

code = "python" 	# python or IDL
opt = 'pop' 		# kin, pop

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
	'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
# gals=[2, 4, 8]
gals=range(10)
# gals = [6]

output_file = "%s/VIMOS_project/analysis/params.txt" % (cc.home_dir)
f = open(output_file, 'w')
for gal in gals:
	galaxy = galaxies[gal]

	tessellation_File = '%s/Data/vimos/analysis/%s/%s/' % (cc.base_dir, galaxy, opt) + \
		'setup/voronoi_2d_binning_output.txt' 

	bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
	n_bins = int(max(bin_num)+1)

	for i in range(n_bins):
		f.write("python errors2.py " + str(gal) + " " + opt + " " + str(i) + 
			"\n")
		
if 'pop' in opt:
	for gal in gals:
		galaxy = galaxies[gal]

		tessellation_File = '%s/Data/vimos/analysis/%s/%s/' % (cc.base_dir, galaxy, opt) + \
			'setup/voronoi_2d_binning_output.txt' 

		bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
		n_bins = int(max(bin_num)+1)

		for i in range(n_bins):
			# if 'kin' in opt:
			f.write("python pop.py " + str(gal) + " " + opt + " " + str(i) + 
				"\n")
f.write("push.sh 'Glamdring VIMOS run finished'")

f.close()
print "Done"
