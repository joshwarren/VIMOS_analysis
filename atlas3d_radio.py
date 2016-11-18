import numpy as np 
from checkcomp import checkcomp
cc=checkcomp()
import matplotlib.pyplot as plt
from collections import Counter



file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
gal_struc, kinestruc, cat= np.loadtxt(file, usecols=(0, 13,14), dtype=str, unpack =True)

file = '%s/Data/atlas3d/XXXI_tableA1.dat' % (cc.base_dir)
gal_radio, FS, radio = np.loadtxt(file, usecols=(0,3,10), dtype=str, unpack=True, skiprows=2)
limit = np.array(['<' in val for val in radio])
radio = np.array([val.replace('<','') for val in radio]).astype(float)

i_struc = np.array([np.where(gal_struc == gal)[0][0] for gal in gal_radio])
kinestruc = kinestruc[i_struc]

FR = np.array(FS=='F')
RR = np.array(['RR' in val for val in kinestruc])
kinestruc = np.array([val.replace('NRR/','') for val in kinestruc])
kinestruc = np.array([val.replace('RR/','') for val in kinestruc])


detected = Counter(kinestruc[~limit])
nondetected = Counter(kinestruc[limit])

f,ax_array = plt.subplots(1,2)
ax_array[0].pie(nondetected.values(), labels=nondetected.keys(), autopct='%1.1f%%',explode=np.ones(len(nondetected.values()))*0.2,startangle=90)
ax_array[1].pie(detected.values(), labels=detected.keys(), autopct='%1.1f%%',explode=np.ones(len(nondetected.values()))*0.2,startangle=90)
[ax.set_aspect('equal') for ax in ax_array]

plt.show()

