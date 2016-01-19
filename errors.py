## ==================================================================
## Propergate uncertainty
## ==================================================================
## warrenj 20150216 Process to progerate the uncertainty using Monty
## Carlo methods to get uncertainty in velocity space.
## warrenj 20160118 Python version of errors.pro

import numpy as np # for array handling
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
import glob # for searching for files

from ppxf import ppxf
import ppxf_util as util

def use_templates(galaxy):
    template_weighting = '/Data/vimosindi/analysis/' + galaxy + \
	'/templates.txt' 

    templatesToUse = np.loadtxt(template_weighting, usecols=(0,), dtype='i')
    return templatesToUse







def errors(i_gal, bin):
## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
# 	galaxy = galaxies[1]
    galaxy = galaxies[i_gal]
    reps = 5000 ## number of monte carlo reps per bin.
    discard = 2
    set_range = [4200,10000]
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

    lamRange_template = CRVAL_temp + [0, CDELT_temp*(NAXIS_temp-1)]

    log_temp_template, logLam_template, velscale = \
        util.log_rebin(lamRange_template, v1)


## ****************************************************************
## NB: shouldn't this be 0.9A as this is resolution?
    FWHM_tem = 2.5     # Miles spectra have a resolution
                           # FWHM of 2.5A.


## Which templates to use are given in use_templates.pro. This is
## transfered to the array templatesToUse.
    templatesToUse = use_templates(galaxy)
    nfiles = len(templatesToUse)
    templates = np.zeros((len(log_temp_template), nfiles))


## Reading the contents of the files into the array templates. 
## Including rebinning them.
    for i in range(nfiles):
        v1, v2 = np.loadtxt(templateFiles[templatesToUse[i]], unpack='True')
## Rebinning templates logarthmically
        log_temp_template, logLam_template, velscale = \
            util.log_rebin(lamRange_template, v2)
## ****************************************************************
## ^^^ this has changed from the last time we called this: we called v1 before...
        templates[:,i] = log_temp_template

    templates /= np.median(log_temp_template)

## ----------========= Reading Tessellation  =============---------

## Reads the txt file containing the output of the binning_spaxels
## routine. 
    x,y,bin_num = np.loadtxt(tessellation_File, usecols=(0,1,2), \
        unpack=True, skiprows=1)

    n_bins = max(bin_num) + 1
## Contains the order of the bin numbers in terms of index number.
    order = np.sort(bin_num)

## ----------========= Reading the spectrum  =============---------

    dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
        "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % (galaxy)) 
        
    galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 0, header=True)
    galaxy_noise = pyfits.getdata(dataCubeDirectory[0], 1)

## write key parameters from header - can then be altered in future	
    CRVAL_spec = header['CRVAL3']
    CDELT_spec = header['CD3_3']
    s = galaxy_data.shape

    rows_to_remove = range(discard)
    rows_to_remove.extend([s[1]-1-i for i in range(discard)])
    cols_to_remove = range(discard)
    cols_to_remove.extend([s[2]-1-i for i in range(discard)])

    galaxy_data = np.delete(galaxy_data, rows_to_remove, axis=1)
    galaxy_data = np.delete(galaxy_data, cols_to_remove, axis=2)
    galaxy_noise = np.delete(galaxy_noise, rows_to_remove, axis=1)
    galaxy_noise = np.delete(galaxy_noise, cols_to_remove, axis=2)


    n_spaxels = len(galaxy_data[0,0,:])*len(galaxy_data[0,:,0])

## ----------========== Spatially Binning =============---------
    spaxels_in_bin = np.where(bin_num == bin)[0]
    n_spaxels_in_bin = len(spaxels_in_bin)

    bin_lin_temp = np.zeros(s[0])
    bin_lin_noise_temp = np.zeros(s[0])

    for i in range(n_spaxels_in_bin):
        x_i = x[spaxels_in_bin[i]]
        y_i = y[spaxels_in_bin[i]]
        for k in range(s[0]):
            bin_lin_temp[k] += galaxy_data[k,y_i,x_i]
            bin_lin_noise_temp[k] += galaxy_noise[k,y_i,x_i]

    bin_lin_noise_temp = np.sqrt(bin_lin_noise_temp)

## --------======== Finding limits of the spectrum ========--------
## limits are the cuts in pixel units, while lamRange is the cuts in
## wavelength unis.
    gap=12
    ignore = int((5581 - CRVAL_spec)/CDELT_spec) + np.arange(-gap,gap)  
    ignore2 = int((5199 - CRVAL_spec)/CDELT_spec) + np.arange(-gap,gap) 

## h is the spectrum with the peak enclosed by 'ignore' removed.
    h = np.delete(bin_lin_temp, ignore)
    h = np.delete(h,ignore2)
    






















## Change to pixel units
#    if keyword_set(range) THEN range = FIX((range - CRVAL_spec)/CDELT_spec)
















































##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    errors(0,0)
