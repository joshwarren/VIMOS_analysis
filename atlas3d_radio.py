import numpy as np 
from checkcomp import checkcomp
cc=checkcomp()
import matplotlib.pyplot as plt



file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
gal_struc, kinestruc, cat= np.loadtxt(file, usecols=(0, 13,14), dtype=str, unpack =True)
#kinestruc=cat

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


labels = np.unique(kinestruc)
detected = np.zeros(len(labels))
undetected = np.zeros(len(labels))
explode = np.zeros(len(labels))
for i, l in enumerate(labels):
	detected[i] = list(kinestruc[~limit]).count(l)
	undetected[i] = list(kinestruc[limit]).count(l)
	if l=='KDC' or l=='KT' or l=='CRC':
		explode[i]=0.2
	if l=='c':
		explode[i]=0.2



f,ax_array = plt.subplots(1,2)
p0 = ax_array[0].pie(undetected, #labels=np.round(undetected/len(kinestruc[limit]),3), 
	autopct='%1.1f%%', explode=explode,startangle=90)
t0 = ax_array[0].set_title('No Radio Detected')
t0.set_position([.5, 1.05])
ax_array[0].text(-1,-1.3,'Sample size: %i' % (int(np.sum(limit))))
ax_array[0].text(-1,-1.4,'Disterbed objects: %i' % (int(np.sum(undetected[explode!=0]))))

p1 = ax_array[1].pie(detected, #labels=np.round(detected/len(kinestruc[~limit]),3), 
	autopct='%1.1f%%', explode=explode,startangle=90)
t1 = ax_array[1].set_title('Radio Detected')
t1.set_position([.5, 1.05])
ax_array[1].text(-1,-1.3,'Sample size: %i' % (int(np.sum(~limit))))
ax_array[1].text(-1,-1.4,'Disterbed objects: %i' % (int(np.sum(detected[explode!=0]))))


[ax.set_aspect('equal') for ax in ax_array]
leg = plt.legend(labels)
f.tight_layout()
f.subplots_adjust(left=0.075, bottom=0, right=1, top=1, wspace=0, hspace=0)


plt.show()

