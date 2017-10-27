import numpy as np 
import matplotlib.pyplot as plt
from prefig import Prefig
Prefig(size=(30,30))
from checkcomp import checkcomp
cc = checkcomp()

fig, axs = plt.subplots(3, 3, sharex=True, sharey='row')

file = '%s/Data/vimos/analysis/pop_sigma.txt' % (cc.base_dir)
value = np.loadtxt(file, unpack=True, skiprows=2, usecols=range(1,25))
galaxy_gals = np.loadtxt(file, unpack=True, skiprows=2, usecols=range(0,), 
	dtype=str)

Atlas_file = '%s/Data/atlas3d/XV_table1.dat' % (cc.base_dir)
atlas3d_sigma = np.loadtxt(Atlas_file, usecols=(1,), unpack=True)
atlas3d_gals = np.loadtxt(Atlas_file, usecols=(0,), unpack=True, dtype=str)
mask = np.array([g not in ['NGC4268, PGC170172', 'NGC1222', 'NGC4684', 
	'NGC5273', 'NGC3998', 'PGC029321'] for g in atlas3d_gals])


for i, str_ap in enumerate(['R_e/8', 'R_e/2', 'R_e']):
	sigma = np.log10(value[2*8+0])
	e_sigma = np.abs(value[2*8+1]/value[2*8+1]/np.log(10))

	Atlas_file = '%s/Data/atlas3d/XXX_table%i.dat' % (cc.base_dir,i+1)

	for j, str_ssp in  enumerate(['age', 'metallicity', 'alpha']):
		ssp = value[i*8+j*2+2]
		e_ssp = value[i*8+j*2+3]

		atlas3d_ssp, atlas3d_e_ssp = np.loadtxt(Atlas_file, unpack=True, 
			usecols=(2*j+9, 2*j+10))

		axs[j, i].errorbar(sigma, ssp, xerr=e_sigma, yerr=e_ssp, c='r', 
			fmt='.')
		axs[j, i].text(0.1, 0.9, str_ap, va='center', 
			transform=axs[j,i].transAxes)

		axs[j, i].errorbar(atlas3d_sigma[mask], atlas3d_ssp[mask], 
			yerr=atlas3d_e_ssp[mask], c='k', fmt='.')

		if str_ssp == 'age':
			axs[j, i].set_yscale('log')
			axs[j, i].set_ylim(0.3, 20)

for a in axs.flatten():
	a.tick_params(top=True, bottom=True, left=True, right=True, direction='in',
		which='both' )
	a.set_xlim([1.6, 2.6])
fig.subplots_adjust(hspace=0, wspace=0)

fig.savefig('%s/Data/muse/analysis/ssp_sigma.png' % (cc.base_dir))
