## ==================================================================
## Propergate uncertainty
## ==================================================================
## warrenj 20150216 Process to progerate the uncertainty using Monty
## Carlo methods to get uncertainty in velocity space.
## warrenj 20160118 Python version of errors.pro

import numpy as np # for array handling
import glob # for searching for files
import pyfits # reads fits files (is from astropy)

def errors(i_gal, bin):
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
# 	galaxy = galaxies[1]
    galaxy = galaxies[i_gal]
    reps = 5000 ## number of monte carlo reps per bin.
    discard = 2
    range = [4200,10000]
    c = 299792.458
#    z = 0.01 # redshift to move galaxy spectrum to its rest frame 
#    vel = 114.0d # Initial estimate of the galaxy velocity and
#    sig = 269.0d # velocity dispersion in km/s in the rest frame
    FWHM_gal = 4*0.571 # The fibre FWHM on VIMOS is
                           # about 4px with a dispersion of
                           # 0.571A/px. (From: http://www.eso.org
                           # /sci/facilities/paranal/instruments
                           # /vimos/inst/ifu.html)
 
    moments = 4 # number of componants to calc with ppxf (see 
                # keyword moments in ppxf.pro for more details)
    degree = 4  # order of addative Legendre polynomial used to 
                #; correct the template continuum shape during the fit 
## File for output: an array containing the calculated dynamics of the
## galaxy. 



#    dir = '~/'
    dir = '/Data/vimosindi/'
#    dir2 = '~/'
    dir2 = '/Data/idl_libraries/'


    data_file =  "/Data/vimosindi/analysis/galaxies.txt"
    # different data types need to be read separetly
    z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,2,3,4,5))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    vel = vel_gals[i_gal]
    sig = sig_gals[i_gal]
    z = z_gals[i_gal]



    tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output2.txt"



    FWHM_gal = FWHM_gal/(1+z) # Adjust resolution in Angstrom

## ----------===============================================---------
## ----------=============== Run analysis  =================---------
## ----------===============================================---------


## ----------=============== Miles library ================---------
# Finding the template files
    templatesDirectory = dir2 + 'ppxf/MILES_library/'
    templateFiles = glob.glob(templatesDirectory + \
		'm0[0-9][0-9][0-9]V') #****** no longer have nfiles.


# v1 is wavelength, v2 is spectrum
    v1, v2 = np.loadtxt(templateFiles[0], unpack='True')

# Using same keywords as fits headers
    CRVAL_temp = v1[0]		# starting wavelength
    NAXIS_temp = np.shape(v2)[0]   # Number of entries
# wavelength increments (resolution?)
    CDELT_temp = (v1[NAXIS_temp-1]-v1[0])/(NAXIS_temp-1)





























































##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    errors(0,0)
