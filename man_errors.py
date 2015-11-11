## ==================================================================
## Manage the uncertain results from glamdring
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 

import numpy as np # for reading files
import array

def man_errors():
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
        'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    galaxy = galaxies[8]
    wav_range = "4200-/"


# ------------========== File 1: vel and sig ==========----------
    glamdring_file = "/Data/vimosindi/analysis/%s/glamdring_result.txt" % (
        galaxy)

    bins, reps, vel, sig  =np.loadtxt(glamdring_file, unpack=True, 
        usecols=(1,2,3,4))
    i_line = np.arange(len(bins),dtype='int') 
    n_bins = np.max(bins)+1

    repeated_bins = np.zeros((1,2), dtype=np.int)
    for bin in range(int(max(bins))):
        x = np.where(bins == bin)
        if len(x[0]) != 5000:
            repeated_bins = np.append(repeated_bins, 
                [[bin, len(x[0])-5000]], axis=0)

# remove empty first line
    repeated_bins = np.delete(repeated_bins,0,0)
    for bin in repeated_bins[:,0]:
        a = np.where(bins == bin)[0]
        d = i_line[a]
        b = np.sort(reps[a])
        c = b[0:-2] - b[1:-1]
        e = b[np.where(c != -1)[0]].astype(int)
#        e = np.append(e,e+1)
        f = d[e]

        bins = np.delete(bins, f)
        reps = np.delete(reps, f)
        vel = np.delete(vel, f)
        sig = np.delete(sig, f)


    v = np.zeros(n_bins)
    v_s = np.zeros(n_bins)
    s = np.zeros(n_bins)
    s_s = np.zeros(n_bins)



    for bin in range(len(v)):
        i_bin = np.where(bins == bin)[0]
        v[bin] = np.mean(vel[i_bin])
        s[bin] = np.mean(sig[i_bin])
        v_s[bin] = np.std(vel[i_bin])
        s_s[bin] = np.std(sig[i_bin])

    v_file = "/Data/vimosindi/analysis/%s/results/%sgal_vel_uncert.dat" % (
        galaxy, wav_range)
    s_file = "/Data/vimosindi/analysis/%s/results/%sgal_sigma_uncert.dat" % (
        galaxy, wav_range)


# ------------========== File 2: h3 and h4 ===========----------
    glamdring_file = "/Data/vimosindi/analysis/%s/glamdring_result2.txt" % (
        galaxy)

    bins, reps, h3s, h4s  =np.loadtxt(glamdring_file, unpack=True, 
        usecols=(1,2,3,4))
    i_line = np.arange(len(bins),dtype='int') 


    repeated_bins = np.zeros((1,2), dtype=np.int)
    for bin in range(int(max(bins))):
        x = np.where(bins == bin)
        if len(x[0]) != 5000:
            repeated_bins = np.append(repeated_bins, 
                [[bin, len(x[0])-5000]], axis=0)

# remove empty first line
    repeated_bins = np.delete(repeated_bins,0,0)
    for bin in repeated_bins[:,0]:
        a = np.where(bins == bin)[0]
        d = i_line[a]
        b = np.sort(reps[a])
        c = b[0:-2] - b[1:-1]
        e = b[np.where(c != -1)[0]].astype(int)
#        e = np.append(e,e+1)
        f = d[e]

        bins = np.delete(bins, f)
        reps = np.delete(reps, f)
        h3s = np.delete(h3s, f)
        h4s = np.delete(h4s, f)


    h3 = np.zeros(n_bins)
    h3_s = np.zeros(n_bins)
    h4 = np.zeros(n_bins)
    h4_s = np.zeros(n_bins)



    for bin in range(len(h3)):
        i_bin = np.where(bins == bin)[0]
        h3[bin] = np.mean(h3s[i_bin])
        h4[bin] = np.mean(h4s[i_bin])
        h3_s[bin] = np.std(h3s[i_bin])
        h4_s[bin] = np.std(h4s[i_bin])

    h3_file = "/Data/vimosindi/analysis/%s/results/%sgal_h3_uncert.dat" % (
        galaxy, wav_range)
    h4_file = "/Data/vimosindi/analysis/%s/results/%sgal_h4_uncert.dat" % (
        galaxy, wav_range)


# ------------============= Save outputs =============----------
    f_v = open(v_file, 'w')
    f_s = open(s_file, 'w')
    f_h3 = open(h3_file, 'w')
    f_h4 = open(h4_file, 'w')
    for i in range(len(v)):
        f_v.write(str(v[i]) + '    ' + str(v_s[i]) + '\n')
        f_s.write(str(s[i]) + '    ' + str(s_s[i]) + '\n')
        f_h3.write(str(h3[i]) + '    ' + str(h3_s[i]) + '\n')
        f_h4.write(str(h4[i]) + '    ' + str(h4_s[i]) + '\n')
    



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':
    man_errors()
