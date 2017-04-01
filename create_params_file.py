import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

code = "python" 	# python or IDL
opt = 'pop' 		# kin, abs or pop

galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 
	'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
gals=[2, 4, 8]

output_file = "params.txt"
f = open(output_file, 'w')
# for gal in gals:
for gal in range(10):
	galaxy = galaxies[gal]

	tessellation_File = '%s/Data/vimos/analysis/%s/voronoi_2d_binning_output_%s.txt' % (
		cc.base_dir, galaxy, opt)

	bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
	n_bins = int(max(bin_num)+1)


	for i in range(n_bins):
		if code == "IDL":
			f.write("idl -e 'errors, " + str(gal) + ", " + str(i) + "' \n")
		if code == "python":
			if opt == "kin":
				f.write("python errors2.py " + str(gal) + " " + str(i) + "\n")
			elif opt == "abs":
				f.write("python errors3.py " + str(gal) + " " + str(i) + "\n")
			elif opt == "pop":
				f.write("python pop.py " + str(gal) + " " + str(i) + "\n")
f.write("push 'Glamdring run finished'")

f.close()
print "Done"
