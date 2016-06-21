from astropy.io import fits
import glob


def radio():
    dir = "/Users/warrenj/vla/"

    fs = glob.glob(dir+"*.FITS")

    f, ax = plt.subplots(2,2, sharex='col', sharey='row')

    
    for i in range(len(fs)):
        r = fits.open(f[i])
        r_img = r[0].data
        r_head = r[0].header


        ax[i] = plt.contour(r_img)

    f.show()


        

        

    
    








if __name__ == '__main__':
    radio()
