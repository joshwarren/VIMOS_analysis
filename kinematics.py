## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import glob # for searching for files
import pyfits # reads fits files (is from astropy)
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import math # for sine functions
import matplotlib.pyplot as plt # used for plotting
import matplotlib.axes as ax # for adding text onto images
from scipy.optimize import curve_fit # for fitting a gaussian
#---------------------------------------------------------------------------
def spxToKpc(x, z):
    H0 = 70.4 # Mpc / kms^-1
    val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
    return val

def kinematics(galaxy, discard=0, wav_range="", 
    corrections=0, plots=False):


    data_file =  "/Data/vimosindi/analysis/galaxies.txt"
    data_file2 =  "/Data/vimosindi/analysis/galaxies2.txt"
    # different data types need to be read separetly
    z_gals, vel_gals, sig_gals, x_gals, y_gals, SN_gals = np.loadtxt(data_file, 
        unpack=True, skiprows=1, usecols=(1,2,3,4,5,6))
    galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
    i_gal = np.where(galaxy_gals==galaxy)[0][0]
    z = z_gals[i_gal]




    ellip_gals, star_psi_gals, gas_psi_gals = np.loadtxt(data_file2,
        unpack=True, skiprows=1, usecols=(1,2,3))



    #discard = 2 # rows of pixels to discard- must have been the same 
            #    for all routines 
    #wav_range = "4200-"
    #corrections = [[13,9],[35,1]]
    find_Re = False # fits radial profile with a gaussian





    if wav_range:
        wav_range_dir = wav_range + "/"
    else:
        wav_range_dir = ""
        wav_range = ""

    dataCubeDirectory = glob.glob("/Data/vimosindi/reduced/%s/cube/" \
        "*crcl_oextr1*vmcmb_darc_cexp_cube.fits" % (galaxy)) 

    tessellation_File = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output.txt"
    tessellation_File2 = "/Data/vimosindi/analysis/%s/" %(galaxy) +\
        "voronoi_2d_binning_output2.txt"

    output_v = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_vel.dat" % (wav_range_dir)
    output_sigma = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_sigma.dat" % (wav_range_dir)
    output_OIII = "/Data/vimosindi/analysis/%s/results/" % (galaxy) +\
        "%sgal_OIII.dat" % (wav_range_dir)

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


    galaxy_data_error = pyfits.getdata(dataCubeDirectory[0], 1)
    galaxy_data_error = np.delete(galaxy_data_error, rows_to_remove, axis=1)
    galaxy_data_error = np.delete(galaxy_data_error, cols_to_remove, axis=2)

    galaxy_data_error = np.sum(galaxy_data_error, axis=0)
    ##galaxy_data_error += galaxy_data

    galaxy_data_error /= np.median(galaxy_data)
    galaxy_data /= np.median(galaxy_data)

#print galaxy_data[12:15,8:11]
#print galaxy_data[34:,0:3]
#galaxy_data[13,9]=np.median(galaxy_data)
#galaxy_data[35,1]=np.median(galaxy_data)
#galaxy_data_error[13,9]=np.median(galaxy_data_error)
#galaxy_data_error[35,1]=np.median(galaxy_data_error)
#print np.min(galaxy_data)
#print np.max(galaxy_data)
    if corrections:
        for i in range(len(corrections)):
            galaxy_data[corrections[i][0],corrections[i][1]] = \
                np.median(galaxy_data)
            galaxy_data_error[corrections[i][0],corrections[i][1]] = \
                np.median(galaxy_data_error)

# ------------============= Fit photometry ===============----------
    save_to = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%splots/photometry_%s.png" % (wav_range_dir, wav_range)
    f = find_galaxy(galaxy_data, quiet=True, plot=plots, 
        galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
    #f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
    x_gals[i_gal] = f.xpeak
    y_gals[i_gal] = f.ypeak
    ellip_gals[i_gal] = f.eps
    print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
    print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))






# ------------================ Kinemetry =================----------
# ------------======== Reading the velocity field ========----------
# Read tessellation file
    x, y, bin_num, xBin, yBin = np.loadtxt(tessellation_File, unpack=True, 
        skiprows = 1) 

    xBar, yBar = np.loadtxt(tessellation_File2, unpack=True, 
        skiprows = 1) 

    xBar += -(40.0-2*discard)/2.0
    yBar += -(40.0-2*discard)/2.0

    v_field = np.loadtxt(output_v, unpack=True)
    v_field -= np.median(v_field)



# ------------============== Fit kinemetry ===============----------
    save_to = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%s/plots/stellar_kinematics_%s.png" % (wav_range_dir, wav_range)

    k = fit_kinematic_pa(xBar, yBar, v_field, quiet=True, plot=plots, \
        sav_fig=save_to)

    print "PA_kin: " + str(k[0]) + "+/-" + str(k[1]/3)


# ------------============== Misalignment ================----------
    phot = math.radians(90-f.theta)
    kine = math.radians(k[0])


    mis = math.asin(abs(math.sin(phot-kine)))
    mis = math.degrees(mis)
    star_psi_gals[i_gal] = mis
    print "Psi: " + str(mis)


# ------------================= Hot gas ================----------
    OIII_vel = np.loadtxt(output_OIII, unpack=True)
    OIII_vel -= np.median(OIII_vel)

    save_to = "/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%s/plots/OIII_kinematics_%s.png" % (wav_range_dir, wav_range)
    gas_k = fit_kinematic_pa(xBar, yBar, OIII_vel, quiet=True, 
        plot=plots, sav_fig=save_to)

    print "Gas_PA_kin: " + str(gas_k[0]) + "+/-" + str(gas_k[1]/3)

# ------------============== Misalignment ================----------
    gas_kine = math.radians(gas_k[0])

    gas_mis = math.asin(abs(math.sin(gas_kine-kine)))
    gas_mis = math.degrees(gas_mis)
    gas_psi_gals[i_gal] = gas_mis
    print "Psi: " + str(gas_mis)



# ------------================= Lambda_R ================----------
# Reload v_field and undo changes above
    xBar += (40.0-2*discard)/2.0
    yBar += (40.0-2*discard)/2.0
    v_field = np.loadtxt(output_v, unpack=True)

    r = np.sqrt(np.square(xBar-f.xmed)+np.square(yBar-f.ymed))
    sigma = np.loadtxt(output_sigma, unpack=True)


    order = np.argsort(r)

    x=np.int_(x)
    y=np.int_(y)

    flux_b = np.zeros(len(v_field))
    for i in range(len(v_field)):
        spaxels_in_bin = (bin_num == i).nonzero()
        flux_b[i] = np.sum(galaxy_data[x[spaxels_in_bin],y[spaxels_in_bin]])


# NB: lam is ordered in terms of increasing R.
    lam_num = flux_b[order]*r[order]*abs(v_field[order]) # numerator
    lam_den = flux_b[order]*r[order]*np.sqrt(np.square(v_field[order]) + 
        np.square(sigma[order])) # denominator

# cumulative summation
    lam_num = np.cumsum(lam_num)
    lam_den = np.cumsum(lam_den)

    lam = lam_num/lam_den
    plt.figure()
    plt.title(r"Radial $\lambda_R$ profile")
    plt.xlabel("Radius (kpc)")
    plt.ylabel(r"$\lambda_R$")
    plt.plot(spxToKpc(r[order],z), lam)
    ax =plt.gca()
    plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), verticalalignment='top',
        transform=ax.transAxes)
    plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
        "%s/plots/lambda_R_%s.png" % (wav_range_dir, wav_range), \
        bbox_inches="tight")
    if plots: 
        plt.show()


# ------------============= Effective Radius =============----------
## Fit galaxy radial flux profile with Gaussian and find half-light. 
## How close does this correspond to half of flux in field of view?
    if find_Re:
        r = ((x-f.xmed)/abs(x-f.xmed)) * np.sqrt(np.square(x-f.xmed) + 
            np.square(y-f.ymed))
        order = np.argsort(r)
   
        x=np.int_(x)
        y=np.int_(y)
    
        xaxis = r[order]
        yaxis = galaxy_data[x[order],y[order]]
        error = galaxy_data_error[x[order],y[order]]


        def gaussian(x, a, c):
            val = a * np.exp(-(x)**2 / c**2)
            return val

        yaxis += -min(yaxis)
        #yaxis[612] = 15

        popt,pcov = curve_fit(gaussian, xaxis, yaxis, 
            sigma = np.sqrt(error))

        p = [1.0/10.0, 0.0, -1.0/3.0, 0.0, 1.0, -math.sqrt(math.pi)/4.0]
        roots = np.roots(p)
        Re = popt[1] * np.real(roots[-1])
        Re_error = np.sqrt(pcov[1, 1]) * np.real(roots[-1])
    
        Re_kpc = spxToKpc(Re, z)
        Re_error_kpc = spxToKpc(Re_error, z)
    

        #print("Scale =  %.3f +/- %.3f" % (popt[0], np.sqrt(pcov[0, 0])))
        #print("Sigma =  %.3f +/- %.3f" % (popt[1], np.sqrt(pcov[1, 1])))
        #print "Re = %.3f +/- %.3f" % (Re, Re_error)
        print "Re = (%.3f +/- %.3f)kpc" % (Re_kpc, Re_error_kpc)

        plt.figure()
        plt.title("Radial flux profile")
        plt.xlabel("Radius (kpc)")
        plt.ylabel('Flux (normalised)')
        plt.plot(spxToKpc(xaxis,z),yaxis, 'y')
        xm = np.linspace(-30., 30., 100)  # 100 evenly spaced points
        plt.plot(spxToKpc(xm,z), gaussian(xm, popt[0], popt[1]), 'r')
        plt.axvline(Re_kpc)
        plt.axvline(-Re_kpc)
        plt.savefig("/Data/vimosindi/analysis/%s/results/" % (galaxy) + \
            "%s/plots/radial_light_profile_%s" % (wav_range_dir, wav_range) + \
            ".png", bbox_inches="tight")
        if plots:
            plt.show()

# ------------============= Save outputs =============----------
    f = open(data_file, 'w')
    f.write('Galaxy      z     velocity     velocity dispersion      x    y     Target SN \n')
    
    for i in range(len(galaxy_gals)):
        f.write(galaxy_gals[i] + '    ' + str(z_gals[i]) + '    ' + \
            str(vel_gals[i]) + '    ' + str(sig_gals[i]) + '    ' + \
            str(int(x_gals[i])) + '    ' + str(int(y_gals[i])) + '    ' + \
            str(SN_gals[i]) + '\n')

    f2 = open(data_file2, 'w')
    f2.write('Galaxy      Ellipticity     Stellar Misalignment      Gas Misalignment \n')
    for i in range(len(galaxy_gals)):
        f2.write(galaxy_gals[i] + '    ' + str(ellip_gals[i]) + '    ' + \
            str(star_psi_gals[i]) + '    ' + str(gas_psi_gals[i]) + '\n')

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
    galaxy = 'ngc3557'
    discard = 2
    wav_range = '4200-'
    corrections = [[13,9],[35,1]]

    kinematics(galaxy, discard=discard, wav_range=wav_range, 
        corrections=corrections, plots=True)
