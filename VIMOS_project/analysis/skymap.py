import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt # used for plotting
import numpy as np # for array handling
from adjustText import adjust_text

plt.subplot(111, projection='aitoff')



# Pass the coordinates in radians, since otherwise the
# aitoff projection messes up:
galaxy = ['eso443-g024',
          'ic1459',
          'ic1531',
          'ic4296',
          'ngc0612',
          'ngc1399',
          'ngc3100',
          'ngc3557',
          'ngc7075',
          'pks0718-34']

RA = np.array([195.253313,
      344.294200,
      002.398017,
      204.162721,
      023.4905708,
      054.621179,
      150.1701533,
      167.490221,
      322.887483,
      110.198221])

RA[np.where(RA > 180.0)[0]] = RA[np.where(RA>180)[0]]-360

dec = np.array([-32.441381,
        -36.462221,
        -32.276817,
        -33.965917,
        -36.4933075,
        -35.450742,
        -31.664542,
        -37.539172,
        -38.617958,
        -34.118289])





deg2rad = np.pi/180.0

plt.plot(RA*deg2rad, dec*deg2rad, 'g*')

plt.grid(True)

plt.text(0,-120*deg2rad,'Declination (degrees)', ha='center', va='center')
plt.ylabel('Right Ascension (degrees)')
#plt.set_xlimit(0,360)

#    plt.subplots_adjust(top=0.95,bottom=0.0)


texts =[]
for i in range(len(galaxy)):
#        texts.append(ax.annotate(galaxy[i], (RA[i]*deg2rad, dec[i]*deg2rad)))
    texts.append(plt.text(RA[i]*deg2rad, dec[i]*deg2rad, galaxy[i], bbox={'pad':0, 'alpha':0}, size=7))



#adjust_text(RA*deg2rad, dec*deg2rad, texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5), bbox={'pad':0, 'alpha':0}, size=7)
#plot_map()
plt.show()
