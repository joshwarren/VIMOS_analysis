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
from scipy import ndimage # for gaussian blur
import math
import os
import sys

from ppxf import ppxf
import ppxf_util as util


#-----------------------------------------------------------------------------
def set_lines (lines, logLam_temp, FWHM_gal):
# In this routine all lines are free to have independent intensities.
# One can fix the intensity ratio of different lines (e.g. the [OIII] doublet)
# by placing them in the same emission template
    lam = np.exp(logLam_temp)
#    lines = lines[where((lines gt min(lam)) and (lines lt max(lam)))]
    sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom
    emission_lines = np.zeros((len(logLam_temp),len(lines)))
    for j in range(len(lines)):
        emission_lines[:,j] = np.exp(-0.5*np.power((lam - lines[j])/sigma,2))
    return emission_lines
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
def use_templates(galaxy, glamdring=False):
    if glamdring:
        template_weighting = '/users/warrenj/analysis/' + galaxy + \
	    '/templates.txt' 
    else:
        template_weighting = '/Data/vimosindi/analysis/' + galaxy + \
	    '/templates.txt' 

    templatesToUse = np.loadtxt(template_weighting, usecols=(0,), dtype='i')
    return templatesToUse
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
def determine_goodpixels(logLam, lamRangeTemp, vel, z, gas=False):
# warrenj 20150905 Copied from ppxf_determine_goodPixels.pro
#
# PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of
#	goodPixels to be used as input keyword for the routine
#	PPXF. This is useful to mask gas emission lines or atmospheric
#	absorptions. It can be trivially adapted to mask different
#	lines. 
# 
# INPUT PARAMETERS:
# - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
#     of each pixel of the log rebinned *galaxy* spectrum.
# - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the
#     minimum and maximum wavelength in Angstrom in the stellar
#     *template* used in PPXF. 
# - VEL: Estimate of the galaxy velocity in km/s.
# 
# V1.0: Michele Cappellari, Leiden, 9 September 2005
# V1.01: Made a separate routine and included additional common
#   emission lines. MC, Oxford 12 January 2012
# V1.02: Included more lines. MC, Oxford, 7 Januray 2014

    c = 299792.458 # speed of light in km/s

## 20150617 warrenj Added Telluric lines (tell) at 5199 (is a blended sky
## line)

 
#dv = lines*0+800d # width/2 of masked gas emission region in km/s
    dv = 800 # width/2 of masked gas emission region in km/s
#    flag = bytearray([0]*len(logLam))
    flag = logLam < 0

# Marks telluric line
    tell = 5199
    flag |= (logLam > np.log(tell) - z - dv/c) \
        & (logLam < np.log(tell) - z + dv/c) 

    if not gas:
#                         -----[OII]-----   Hdelta    Hgamma   Hbeta;
        lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, \
#            -----[OIII]----- ----??----   [OI]   
            4958.92, 5006.84, 5528, 5535, 6300.30, \
#           -----[NII]-----  Halpha   -----[SII]-----   
            6548.03, 6583.41,  6562.80, 6716.47, 6730.85])

        for j in range(len(lines)):
            flag |= (logLam > np.log(lines[j]) + (vel- dv)/c) \
                & (logLam < np.log(lines[j]) + (vel+ dv)/c)


    flag |= logLam < np.log(lamRangeTemp[0]) + (vel + 900)/c # Mask edges of
    flag |= logLam > np.log(lamRangeTemp[1]) + (vel - 900)/c # stellar library


    flag[0:3] = 1 # Mask edge of data
    flag[-4:]= 1 # to remove edge effects
    return np.where(flag == 0)[0]
#-----------------------------------------------------------------------------






#-----------------------------------------------------------------------------
def errors(i_gal=None, bin=None):
    if i_gal is None: i_gal=int(sys.argv[1])
    if bin is None: bin=int(sys.argv[2])

## ----------===============================================---------
## ----------============= Input parameters  ===============---------
## ----------===============================================---------
    glamdring = False
    gas = False
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
# 	galaxy = galaxies[1]
    galaxy = galaxies[i_gal]
    reps = 0 ## number of monte carlo reps per bin.
    discard = 2
#    set_range = None
    set_range = np.array([4200,10000])
    c = 299792.458
#    z = 0.01 # redshift to move galaxy spectrum to its rest frame 
#    vel = 114.0d # Initial estimate of the galaxy velocity and
#    sig = 269.0d # velocity dispersion in km/s in the rest frame
    FWHM_gal = 4*0.571 # The fibre FWHM on VIMOS is
                           # about 4px with a dispersion of
                           # 0.571A/px. (From: http://www.eso.org
                           # /sci/facilities/paranal/instruments
                           # /vimos/inst/ifu.html)
 
    stellar_moments = 4 # number of componants to calc with ppxf (see 
                        # keyword moments in ppxf.pro for more details)
    gas_moments = 4
    degree = 4  # order of addative Legendre polynomial used to 
                #; correct the template continuum shape during the fit 
## File for output: an array containing the calculated dynamics of the
## galaxy. 


    if glamdring:
        dir = '/users/warrenj/'
        dir2 = '/users/warrenj/'
    else:
        dir = '/Data/vimosindi/'
        dir2 = '/Data/idl_libraries/'


    data_file = dir + "analysis/galaxies.txt"
    # different data types need to be read separetly
    z_gals, vel_gals, sig_gals, x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,2,3,4,5))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    vel = vel_gals[i_gal]
    sig = sig_gals[i_gal]
    z = z_gals[i_gal]



    tessellation_File = dir + "analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = dir + "analysis/%s/" %(galaxy) +\
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
    FWHM_tem = 2.5 # Miles spectra have a resolution
                               # FWHM of 2.5A.


## Which templates to use are given in use_templates.pro. This is
## transfered to the array templatesToUse.
    templatesToUse = use_templates(galaxy, glamdring)
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

    dataCubeDirectory = glob.glob(dir + "reduced/%s/cube/" \
        "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % (galaxy)) 
        
    galaxy_data, header = pyfits.getdata(dataCubeDirectory[0], 1, header=True)
    galaxy_noise = pyfits.getdata(dataCubeDirectory[0], 2)
    
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
            bin_lin_noise_temp[k] += galaxy_noise[k,y_i,x_i]**2

    bin_lin_noise_temp = np.sqrt(bin_lin_noise_temp)
## --------======== Finding limits of the spectrum ========--------
## limits are the cuts in pixel units, while lamRange is the cuts in
## wavelength unis.
    gap=12
    ignore = int((5581 - CRVAL_spec)/CDELT_spec) + np.arange(-gap+1,gap)  
    ignore2 = int((5199 - CRVAL_spec)/CDELT_spec) + np.arange(-gap+1,gap) 

## h is the spectrum with the peak enclosed by 'ignore' removed.
    h = np.delete(bin_lin_temp, ignore)
    h = np.delete(h,ignore2)

    half = s[0]/2
    a = h/np.median(h) - np.append(h[4:],[0,0,0,0])/np.median(h)
    a = np.where(np.isfinite(a), a, 0)

#    if np.where(np.abs(a[:0.5*half]) > 0.2)[0]:
#        lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
#    else: 
#        lower_limit = 0
#        print str(i_gal) + ', ' + str(bin)

    lower_limit = max(np.where(np.abs(a[:0.5*half]) > 0.2)[0])
    upper_limit = min(np.where(np.abs(a[1.5*half:]) > 0.2)[0])+int(1.5*half)

    if upper_limit > ignore2[0]: upper_limit+=gap 
    if upper_limit > ignore[0]: upper_limit+=gap

    if lower_limit < 0:
        lower_limit = min(np.where(a[:half] != 0)[0]) + 5
        if lower_limit < 0: lower_limit = 0
    else:
        lower_limit +=5

    if upper_limit > s[0]-1 or upper_limit < half:
        upper_limit = s[0]-6 
    else:
        upper_limit += -5

## --------========= Using set_range variable =========--------  
    if set_range is not None:
## Change to pixel units
        set_range = ((set_range - CRVAL_spec)/CDELT_spec).astype(int)
        if set_range[0] > lower_limit: lower_limit = set_range[0] 
        if set_range[1] < upper_limit: upper_limit = set_range[1]

    lamRange = np.array([lower_limit, upper_limit])*CDELT_spec + CRVAL_spec

## ----------========= Writing the spectrum  =============---------
    bin_lin = bin_lin_temp[lower_limit:upper_limit]
    bin_lin_noise = bin_lin_noise_temp[lower_limit:upper_limit]
## ----------======== Calibrating the spectrum  ===========---------
## For calibrating the resolutions between templates and observations
## using the gauss_smooth command
    FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_gal**2)
    sigma = FWHM_dif/2.355/CDELT_temp # Sigma difference in pixels

## smooth spectrum to fit with templates resolution
    bin_lin = ndimage.gaussian_filter1d(bin_lin, sigma)
    bin_lin_noise = ndimage.gaussian_filter1d(bin_lin_noise, sigma)
    
    lamRange = lamRange/(1+z)
## rebin spectrum logarthmically
    bin_log, logLam_bin, velscale = util.log_rebin(lamRange, bin_lin)
    bin_log_noise, logLam_bin, velscale = util.log_rebin(lamRange, np.power(bin_lin_noise,2))
    bin_log_noise = np.sqrt(bin_log_noise)

## Normalis the spectrum
    med_bin = np.median(bin_log)
    bin_log /= med_bin
    bin_log_noise /= med_bin
    noise = bin_log_noise+0.0000000000001


    dv = (logLam_template[0]-logLam_bin[0])*c # km/s
# Find the pixels to ignore to avoid being distracted by gas emission
#; lines or atmospheric absorbsion line.  
    goodPixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z) 
    lambdaq = np.exp(logLam_bin)
    start = [vel, sig] # starting guess
    component = [0]*len(templates[0,:])

## ----------============= Emission lines =================---------
    moments = stellar_moments
    if gas:
        emission_file = dir + 'analysis/emission_line.dat'
        line_name = np.loadtxt(emission_file, unpack=True, skiprows=1, 
            usecols=(1,), dtype=str)
        line_wav = np.loadtxt(emission_file, unpack=True, skiprows=1, 
            usecols=(2,))

        outOfRange = np.where((line_wav < lamRange[0]) | \
                              (line_wav > lamRange[1]))[0]

        line_name = np.delete(line_name, outOfRange)
        line_wav = np.delete(line_wav, outOfRange)

        emission_lines = set_lines(line_wav, logLam_template, FWHM_gal)


        component = component + [1]*len(line_name)
        templates = np.column_stack((templates, emission_lines))
       
        start = [start,start]
        moments = [stellar_moments, gas_moments]
        goodPixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z, 
            gas=True) 




## ----------=========== The bestfit part =================---------


    bin_log_sav = bin_log

    pp = ppxf(templates, bin_log, noise, velscale, start, 
              goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, 
              component=component, lam=lambdaq, plot=True, quiet=False)

## ----------================= The MC part ==================---------
    stellar_output = np.zeros((reps, stellar_moments))
    stellar_errors = np.zeros((reps, stellar_moments))
    if gas:
        gas_output = np.zeros((reps, gas_moments))
        gas_errors = np.zeros((reps, gas_moments))

    for rep in range(reps):
        print rep
        random = np.random.randn(len(noise))
#        gaussian = 1/(np.sqrt(2*math.pi)*noise)*np.exp(-0.5*((random)/noise)**2)
#        add_noise = (random/abs(random))* \
#            np.sqrt((-2*np.power(noise,2))*np.log(gaussian*noise))
        add_noise = random*np.abs(noise)
        bin_log = pp.bestfit + add_noise
    
        ppMC = ppxf(templates, bin_log, noise, velscale, start, goodpixels=goodPixels, moments=moments, degree=degree, vsyst=dv, lam=lambdaq, plot=True, quiet=False, bias=0.1, component=component)

        stellar_output[rep,:] = ppMC.sol[0:stellar_moments][0]
        stellar_errors[rep,:] = ppMC.error[0:stellar_moments][0]
        if gas:
            gas_output[rep,:] = ppMC.sol[0:gas_moments][1]
            gas_errors[rep,:] = ppMC.error[0:gas_moments][1]
## ----------============ Write ouputs to file ==============---------

    if not os.path.exists("%sanalysis/%s/gas_MC/stellar/errors" % (dir,galaxy)):
        os.makedirs("%sanalysis/%s/gas_MC/stellar/errors" % (dir, galaxy))
    if not os.path.exists("%sanalysis/%s/gas_MC/gas/errors" % (dir, galaxy)):
        os.makedirs("%sanalysis/%s/gas_MC/gas/errors" % (dir, galaxy))
    bin_file = "%sanalysis/%s/gas_MC/stellar/%s.dat" % (dir, galaxy, str(bin))
    errors_file = "%sanalysis/%s/gas_MC/stellar/errors/%s.dat" % (dir, galaxy, 
        str(bin))
    gas_file = "%sanalysis/%s/gas_MC/gas/%s.dat" % (dir, galaxy, str(bin))
    gas_errors_file = "%sanalysis/%s/gas_MC/gas/errors/%s.dat" % (dir, galaxy, 
        str(bin))

    f = open(bin_file, 'w')
    e = open(errors_file, 'w')
    g = open(gas_file, 'w')
    ger = open(gas_errors_file, 'w')
    for i in range(reps):
        f.write(str(stellar_output[i,0]) + "   " + \
            str(stellar_output[i,1]) + "   " + str(stellar_output[i,2]) + \
            "   " + str(stellar_output[i,3]) + '\n')
        e.write(str(stellar_errors[i,0]) + "   " + str(stellar_errors[i,1]) + \
            "   " + str(stellar_errors[i,2]) + "   " + \
            str(stellar_errors[i,3]) + '\n')
        if gas:
            g.write(str(gas_output[i,0]) + "   " + str(gas_output[i,1]) + \
                "   " + str(gas_output[i,2]) + "   " + \
                str(gas_output[i,3]) + '\n')
            ger.write(str(gas_errors[i,0]) + "   " + str(gas_errors[i,1]) + \
                "   " + str(gas_errors[i,2]) + "   " + \
                str(gas_errors[i,3]) + '\n')

## save bestfit spectrum
    if not os.path.exists("%sanalysis/%s/gas_MC/bestfit" % (dir, galaxy)):
        os.makedirs("%sanalysis/%s/gas_MC/bestfit" % (dir, galaxy)) 
    bestfit_file = "%sanalysis/%s/gas_MC/bestfit/%s.dat" % (dir, galaxy, 
        str(bin))
   
    s = open(bestfit_file, 'w')
    for i in range(len(pp.bestfit)):
        s.write(str(pp.bestfit[i]) + '\n')

## save input
    if not os.path.exists("%sanalysis/%s/gas_MC/input" % (dir, galaxy)):
        os.makedirs("%sanalysis/%s/gas_MC/input" % (dir, galaxy)) 
    input_file = "%sanalysis/%s/gas_MC/input/%s.dat" % (dir, galaxy, str(bin))
   
    inp = open(input_file, 'w')
    for i in range(len(bin_log_sav)):
        inp.write(str(bin_log_sav[i]) + '\n')


## save bestfit output
    bestfit_file = "%sanalysis/%s/gas_MC/%s.dat" % (dir, galaxy, str(bin))
   
    b = open(bestfit_file, 'w')
    if gas: b.write(str(pp.sol[0][0]) + "   " + str(pp.sol[0][1]) + "   " + \
        str(pp.sol[0][2]) + "   " + str(pp.sol[0][3]) + '\n' + \
        str(pp.sol[1][0]) + "   " + str(pp.sol[1][1]) + "   " + \
        str(pp.sol[1][2]) + "   " + str(pp.sol[1][3]) + '\n')
    else: b.write(str(pp.sol[0]) + "   " + str(pp.sol[1]) + "   " + \
        str(pp.sol[2]) + "   " + str(pp.sol[3]) + '\n')





















##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    errors(0,10) if len(sys.argv)<3 else errors()
