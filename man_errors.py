## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import array


def man_errors():
    galaxy = 'ic1459'
    glamdring_file = "/Data/vimosindi/analysis/%s/glamdring_result.txt" % (galaxy)

#    bins, rep, vel, sig = np.loadtxt(glamdring_file, usecols=(), unpack=True)
    bins, reps, vel, sig  =np.loadtxt(glamdring_file, unpack=True, usecols=(1,2,3,4))




    for bin in range(int(max(bins))):
        x = np.where(bins == bin)
        if len(x[0]) != 5000:
            print bin, len(x[0])


    



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    man_errors()
