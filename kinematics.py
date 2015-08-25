## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa #fit kinemetry
#---------------------------------------------------------------------------
wav_range=None


galaxy = "ngc3557"
discard = 2 # rows of pixels to discard- must have been the same 
            #    for all routines 
wav_range="4200-"




if wav_range:
    wav_range_dir = wav_range + "/"
else:
    wav_range_dir = ""

dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
    "*crcl_oextr1_fluxcal_vmcmb_darc_cexp_cube.fits" % (galaxy)) 

tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output.txt"
tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
"voronoi_2d_binning_output2.txt"

output_v = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
"%sgal_vel.dat" % (wav_range_dir)

# ------------=============== Photometry =================----------
# ------------========== Reading the data cube ===========----------

galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)


s = galaxy_data.shape
rows_to_remove = range(discard)
rows_to_remove.extend([s[1]-1-i for i in range(discard)])
cols_to_remove = range(discard)
cols_to_remove.extend([s[2]-1-i for i in range(discard)])

galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)


galaxy_data = np.sum(galaxy_data, axis=0)
galaxy_data /= np.median(galaxy_data)


# ------------============= Fit photometry ===============----------
f = find_galaxy(galaxy_data, quiet=True)
print f.eps
print f.theta





# ------------================ Kinemetry =================----------
# ------------======== Reading the velocity field ========----------
# Read tessellation file
x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
    skiprows = 1) 


xBar, yBar = np.loadtxt(tessellation_File2, unpack=True, 
    skiprows = 1) 


v_field = np.loadtxt(output_v)
v_field -= np.median(v_field)


# ------------============== Fit kinemetry ===============----------
k = fit_kinematic_pa(xBar, yBar, v_field)#, quiet=True, plot=True)
print k[0]




