import numpy as np
import matplotlib.pyplot as plt # used for plotting
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)


def removerows(gal):
    d = "/Data/vimos/cubes/"
    f = glob.glob(d + gal + ".cube.combined.fits")
    p = pyfits.open(f[0])
    data = p[0].data
    header = p[0].header
    data_new = data[:,0:40,0:40]

    pyfits.writeto(d+gal+".cube.combined.fits", data_new, header=header, clobber=True)

    for ex in range(1, len(p)):
        data_ex = p[ex].data
        header_ex = p[ex].header
        data_ex_new = data_ex[:,0:40,0:40]
        pyfits.append(d+gal+".cube.combined.fits",data_ex_new,header=header_ex)

    p.close()


    
if __name__ == '__main__':
    galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
    for gal in galaxies:
        print gal
        removerows(gal)
