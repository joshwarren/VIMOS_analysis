import numpy as np 
import matplotlib.pyplot as plt
from prefig import Prefig
Prefig(size=(20,20))
from checkcomp import checkcomp
cc = checkcomp()
from sauron_colormap import sauron#2 as sauron


fig, axs = plt.subplots(3, 3, sharex=True, sharey='row')

file = '%s/Data/vimos/analysis/pop_sigma.txt' % (cc.base_dir)
value = np.loadtxt(file, unpack=True, skiprows=2, usecols=range(1,25))
v_gals = np.loadtxt(file, usecols=(0,), unpack=True, dtype=str, skiprows=2)



file = '%s/Data/muse/analysis/pop_sigma.txt' % (cc.base_dir)
m_value = np.loadtxt(file, unpack=True, skiprows=2, usecols=range(1,25))
m_gals = np.loadtxt(file, usecols=(0,), unpack=True, dtype=str, skiprows=2)



Atlas_file = '%s/Data/atlas3d/XV_table1.dat' % (cc.base_dir)
atlas3d_sigma = np.loadtxt(Atlas_file, usecols=(1,), unpack=True)
atlas3d_gals = np.loadtxt(Atlas_file, usecols=(0,), unpack=True, dtype=str)
mask = np.array([g not in ['NGC4268, PGC170172', 'NGC1222', 'NGC4684', 
	'NGC5273', 'NGC3998', 'PGC029321'] for g in atlas3d_gals])


for i, str_ap in enumerate([r'R$_e$/8', r'R$_e$/2', r'R$_e$']):
	sigma = np.log10(value[i*9+0])
	e_sigma = np.abs(value[i*9+1]/value[i*9+0]/np.log(10))
	prob = value[i*9+8]
	print prob

	m_sigma = np.log10(m_value[i*8+0])
	e_m_sigma = np.abs(m_value[i*8+1]/m_value[i*8+0]/np.log(10))

	Atlas_file = '%s/Data/atlas3d/XXX_table%i.dat' % (cc.base_dir,i+1)

	for j, str_ssp in  enumerate(['age', 'metallicity', 'alpha']):
		ssp = value[i*9+j*2+2]
		e_ssp = value[i*9+j*2+3]

		m_ssp = m_value[i*8+j*2+2]
		e_m_ssp = m_value[i*8+j*2+3]

		atlas3d_ssp, atlas3d_e_ssp = np.loadtxt(Atlas_file, unpack=True, 
			usecols=(2*j+9, 2*j+10))

		axs[j, i].errorbar(sigma, ssp, xerr=e_sigma, yerr=e_ssp, c=prob, 
			fmt='o')
		axs[j, i].errorbar(m_sigma, m_ssp, xerr=e_m_sigma, yerr=e_m_ssp, c='r', 
			fmt=x)
		axs[j, i].text(0.1, 0.9, str_ap, va='center', 
			transform=axs[j,i].transAxes)

		axs[j, i].errorbar(atlas3d_sigma[mask], atlas3d_ssp[mask], 
			yerr=atlas3d_e_ssp[mask], c='k', fmt='.')

		if str_ssp == 'age':
			axs[j, i].set_yscale('log')
			axs[j, i].set_ylim(0.3, 20)

		for g in ['ic1459','ic4296','ngc1399']:
			m_m = m_gals == g
			v_m = v_gals == g
			axs[j,i].plot([sigma[v_m], m_sigma[m_m]], [ssp[v_m], m_ssp[m_m]], 'k--')

for a in axs.flatten():
	a.tick_params(top=True, bottom=True, left=True, right=True, direction='in',
		which='both' )
	a.set_xlim([1.6, 2.6])

axs[0,0].set_ylabel('Age (Gyr)')
axs[1,0].set_ylabel(r'Metallicity [Fe/H]$_\odot$')
axs[2,0].set_ylabel(r'Alpha-enhancement [$\alpha$/Fe]$_\odot$')
axs[2,1].set_xlabel(r'log $\sigma$ (km s$^{-1}$)')

# Turn off scientific notation for Age
axs[0,0].set_yticks([1, 10])
import matplotlib as mpl
axs[0,0].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())



fig.subplots_adjust(hspace=0, wspace=0)

fig.savefig('%s/Data/muse/analysis/ssp_sigma.png' % (cc.base_dir))
