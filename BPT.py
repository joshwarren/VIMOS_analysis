from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp
from plot_results import add_
from astropy.io import fits
from errors2 import get_dataCubeDirectory
import os
from prefig import Prefig
from sauron_colormap import sauron

bd = 2.86 # Balmer decrement 
lbd = np.log10(bd)

nd = 12.2300 # grad of bestfit to [NII] vs [NI]
lnd = np.log10(nd)

sd = 0.76 # grad continuum from Ha to Hb
lsd = np.log10(sd)
def NII_to_NI(v):
	a = nd
	b = 2.748e4
	return (v - b) / a

def log_NII_to_NI(v):
	return np.log10(NII_to_NI(10**v))

def e_NII_to_NI(v, s_v):
	a = nd
	b = 0

	s_aa = 0.1389**2 # ; s_ab = 429.6
	# s_bb = 1.417e6	

	f_a = (b - v) / a**2
	f_b = -1. / a
	f_v = -f_b

	# # Varience
	# out = f_a**2 * s_aa + f_b**2 * s_bb + f_v**2 * s_v**2
	# # Add covarience, assuming zero covarience between v and a, b and c
	# out += 2 * f_a * f_b * s_ab

	# Varience
	out = f_a**2 * s_aa + f_v**2 * s_v**2
	return np.sqrt(out) 

def OI_to_OIII(v):
	a = 8.017e-6
	b = 0.3879
	c = 1.369e-6
	# positive root
	return (-b + np.sqrt(b**2 - 4 * a * (c - v))) / 2 / a 

def log_OI_to_OIII(v):
	return np.log10(OI_to_OIII(10**v))

def e_OI_to_OIII(v, s_v):
	a = 8.017e-6
	b = 0.3879
	c = 1.369e-6

	s_aa = 4.842e-15 ; s_ab = -6.034e-11 ; s_ac = 3.262e-16
	s_bb = 1.661e-6 ; s_bc = -8.98e-12
	s_cc = 2.681e-16

	# Partial derrivatives
	f_a = b / 2 * a**2 + (-b**2 / 2 / a**3 - (c - v) / a**2) / 2
	f_b = -1. / 2 / a + b / 2 / a / np.sqrt(b**2 + 4 * a * (c - v))
	f_c = 1 / np.sqrt(b**2 + 4 * a * (c - v))
	f_v = -f_c

	# Variance (remember: s_aa = s_a^2)
	out = f_a**2 * s_aa + f_b**2 * s_bb + f_c**2 * s_cc + f_v**2 * s_v**2

	# Add covarience, assuming zero covarience between v and a, b and c
	out += 2 * f_a * f_b * s_ab + 2 * f_a * f_c * s_ac + 2 * f_b * f_c * s_bc
	return np.sqrt(out)

def OIII_to_OI(v):
	a = 8.017e-6
	b = 0.3879
	c = 1.369e-6
	# positive root
	return a * v**2 + b * v + c

# def e_OIII_to_OI(v, s_v):
# 	a = 8.017e-6
# 	b = 0.3879
# 	c = 1.369e-6

# 	s_aa = 4.842e-15 ; s_ab = -6.034e-11 ; s_ac = 3.262e-16
# 	s_bb = 1.661e-6 ; s_bc = -8.98e-12
# 	s_cc = 2.681e-16

# 	# Partial derrivatives
# 	f_a = b / 2 * a**2 + (-b**2 / 2 / a**3 - (c - v) / a**2) / 2
# 	f_b = -1. / 2 / a + b / 2 / a / np.sqrt(b**2 + 4 * a * (c - v))
# 	f_c = 1 / np.sqrt(b**2 + 4 * a * (c - v))
# 	f_v = -f_c

# 	# Variance (remember: s_aa = s_a^2)
# 	out = f_a**2 * s_aa + f_b**2 * s_bb + f_c**2 * s_cc + f_v**2 * s_v**2

# 	# Add covarience, assuming zero covarience between v and a, b and c
# 	out += 2 * f_a * f_b * s_ab + 2 * f_a * f_c * s_ac + 2 * f_b * f_c * s_bc
# 	return np.sqrt(out)


	

def add_grids(ax, xline, yline, x_Ha=False, y_Ha=False):
	# Starbursts - No tables in paper
	d = '%s/models/Dopita_starbursts' % (cc.home_dir)
	# Z = 0.2,0.4,1,2 solar by q = 0.5,1,2,4,8,15,30 X 10^7 cm/s
	NI = np.zeros((4,7))
	OIII = np.zeros((4,7))

	# AGN - I can't find any reference to N_e = 100 cm^-3 in paper
	d = '%s/models/Groves_AGN' % (cc.home_dir)
	# alpha = -1.7, -1.4, -1.2 by log U = -3.0, -2.6, -2.3, -2.0, -1.6, -1.3, -1.0 
	# by n_e = 100, 1000
	# NI = np.zeros((3, 7, 2))
	# OIII = np.zeros((3, 7, 2))

	"""
	#???
	# alpha = 1.2, -1.4, -1.7 by log U = -1, -2, -3, -4 by N_e = 100, 1000 cm^{-3}
	NI = np.zeros((3, 4, 1))
	OIII = np.zeros((3, 4, 1))
	
	for i, n in enumerate(['dusty']):#, 'dust_free', 'undeplete_dust_free']):
		with open(d+'/%s.txt' % (n), 'r') as f:
			rows = f.read().splitlines()
			new_table = np.where(['alpha = ' in r for r in rows])[0]
			# new_table = np.append(new_table, -1) # removes alpha = -2.0
			for j in range(len(new_table)-1):
				table = rows[new_table[j]:new_table[j+1]]
				OIII_row = np.where(['[O III] \\lambda 5007' in r for r in table]
					)[0][0]
				if n == 'dusty': 
					start_col = 2
				else:
					start_col = 1
				OIII[j, :, i] = np.array(table[OIII_row].split('\t')[start_col:]
					).astype(float)
				NI_row = np.where(['[N I] \\lambda 5200' in r for r in table])[0][0]
				NI[j, :, i] = np.array(table[NI_row].split('\t')[start_col:]
					).astype(float)

	for a in range(3):
		# for n in range(1):
		# 	ax.plot(np.log10(NI[a,:,n]), np.log10(OIII[a,:,n]), 'r:')
		# for u in range(4):
		# 	ax.plot(np.log10(NI[a,u,:]), np.log10(OIII[a,u,:]), 'r-')
		for n in range(1):
			ax.plot(NI[a,:,n], OIII[a,:,n], 'r:')
		for u in range(4):
			ax.plot(NI[a,u,:], OIII[a,u,:], 'r-')
	"""

	# LINER - from Allen 2008 ApJS 178 20
	d = '%s/models/Allen_shock_models' % (cc.home_dir)
	# V_s = 150, 200, 300, 500, 750, 1000 km/s by b = 0.5, 1, 2, 4
	x = np.zeros((6,4))
	y = np.zeros((6,4))

	for i, b in enumerate(['0_5', '1', '2', '4']):
		lines = np.genfromtxt('%s/M_n1_b%s_s_lines.txt' % (d, b), unpack=True, 
			usecols=(3,6,8,12,19,29,36), skip_footer=1)
		wav = np.array(lines[0,:])
		lines = lines[1:,:]

		atom, species = np.genfromtxt('%s/M_n1_b%s_s_lines.txt' % (d, b), 
			unpack=True, usecols=(1,2), dtype=str, skip_footer=1)

		if yline == '[OIII]':
			y_row = (atom=='O') * (species=='III') * (wav==5006.77)
		if xline == '[NI]':
			x_row = (atom=='N') * (species=='I') * (wav==5197.82)
			x_row2 = (atom=='N') * (species=='I') * (wav==5200.17)
		if xline == '[OI]':
			x_row = (atom=='O') * (species=='I') * (wav==6300.2)
		if xline == '[NII]':
			x_row = (atom=='N') * (species=='II') * (wav==6583.34)
		if xline == '[SII]':
			x_row = (atom=='S') * (species=='II') * (wav==6716.31)
			x_row2 = (atom=='S') * (species=='II') * (wav==6730.68)
			


		x[:, i] = lines[:, x_row].flatten()
		if 'x_row2' in locals().keys():
			x[:, i] += lines[:, x_row2].flatten()

		y[:, i] = lines[:, y_row].flatten()
		if 'y_row2' in locals().keys():
			y[:, i] += lines[:, y_row2].flatten()

		Ha_row = (atom=='H') * (species=='I') * (wav==6562.8)
		if y_Ha:
			y[:, i] /= lines[:, Ha_row].flatten()
		if x_Ha:
			x[:, i] /= lines[:, Ha_row].flatten()

	for i in range(6):
		ax.plot(np.log10(x[i, :]), np.log10(y[i, :]), 'g-')
	for i in range(4):
		ax.plot(np.log10(x[:, i]), np.log10(y[:, i]), 'g:')






def BPT(galaxy, D=None, opt='pop', return_sauron_diagnotics=False):
	print '   BPT'

	analysis_dir = "%s/Data/vimos/analysis" % (cc.base_dir)
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(galaxiesFile, unpack=True, 
		skiprows=1, usecols=(4,5), dtype=int)
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),
		dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)

	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()

	if 'Hbeta' not in D.e_line.keys(): 
		if return_sauron_diagnotics:
			return [np.nan], [np.nan]
		else:
			return D

# ------------=============== MEx diagram =================----------
	# from Atlas3D XXXI (Section 6.2.1)
	Prefig()
	fig, ax = plt.subplots()
	# Factor of 1.35 is to only use 5007 line, and not doublet counterpart
	y = np.log10(D.e_line['[OIII]5007d'].flux/1.35/D.e_line['Hbeta'].flux)
	return_y = np.array(y)

	y_err = y * np.sqrt((D.e_line['[OIII]5007d'].flux.uncert/
		D.e_line['[OIII]5007d'].flux)**2 + (D.e_line['Hbeta'].flux.uncert/
		D.e_line['Hbeta'].flux)**2)/np.log(10)

	x = D.components['stellar'].plot['sigma']

	large_err = y_err**2 > 1
	m = ~large_err * (D.e_line['[OIII]5007d'].equiv_width / 1.35 < 0.8)
	ax.errorbar(x[m], y[m], c='b', xerr = x.uncert[m], yerr=y_err[m], 
		fmt='.')
	ax.errorbar(np.nan, np.nan, xerr=np.nan, yerr=np.nan, 
		label='EW([OIII]) < 0.8')

	m = ~large_err * (D.e_line['[OIII]5007d'].equiv_width / 1.35 >= 0.8)
	ax.errorbar(x[m], y[m], c='r', xerr = x.uncert[m], yerr=y_err[m], 
		fmt='.')
	ax.errorbar(np.nan, np.nan, xerr=np.nan, yerr=np.nan, 
		label=r'EW([OIII]) $\geq 0.8$')

	ax.axvline(70, c='k')
	ax.axhline(np.log10(0.5), xmin=70./min(max(
		D.components['stellar'].plot['sigma'][m]), 500), c='k')
	ax.axhline(np.log10(1), xmin=70./min(max(
		D.components['stellar'].plot['sigma'][m]), 500), c='k')

	x_line = np.arange(70,1000,1)
	y_line = 1.6*10**-3 * x_line + 0.33
	ax.plot(x_line, y_line, c='k')

	ax.set_xlim([0, min(max(D.components['stellar'].plot['sigma'][m]), 500)])
	ax.set_ylim([-1.2, 1.5])

	ylim = ax.get_ylim()
	yrange = ylim[1] - ylim[0]
	ax.text(50, ylim[0] +0.96 * yrange, 'SF')
	ax.text(75, 0.55, 'Seyfert 2')
	ax.text(75, 0.15, 'LINER')
	ax.text(75, -0.23, 'Transition')
	ax.text(75, -0.5, 'Passive')

	ax.set_xlabel(r'$\sigma_\ast$')
	ax.set_ylabel(r'log [OIII]$\lambda$5007/H$\,\beta$')
	ax.set_title('Mass-excitation (MEx) diagnotics for %s' % (galaxy.upper()))
	ax.legend(facecolor='w')

	fig.savefig('%s/plots/MEx.png' % (output))
	plt.close(fig)

# ------------============= SAURON diagram =================----------
	# from SAURON XVI
	if '[NI]d' in D.e_components:
		Prefig()
		fig, ax = plt.subplots()

		# y as in MEx above
		x = np.log10(D.e_line['[NI]d'].flux/D.e_line['Hbeta'].flux)
		return_x = np.array(x)

		x_err = x * np.sqrt((D.e_line['[NI]d'].flux.uncert/
			D.e_line['[NI]d'].flux)**2 +
			(D.e_line['Hbeta'].flux.uncert/D.e_line['Hbeta'].flux)**2)/\
			np.log(10)

		large_err = (x_err > 0.5) + (y_err > 0.5)

		ax.errorbar(x, y, yerr=y_err, xerr=x_err, fmt='.')

		ax.set_xlim([-2.5, 0.5])
		ax.set_ylim([-1.5, 1.5])
		ax.set_xlabel(r'log [NI]$\lambda\lambda$5197,5200/H$\,\beta$')
		ax.set_ylabel(r'log [OIII]$\lambda$5007/H$\,\beta$')
		ax.set_title('SAURON diagnotics for %s' % (galaxy.upper()))

		# Add BPT line - transformed to NI
		xlim = ax.get_xlim()
		x_line = np.linspace(xlim[0], xlim[1], 100)
		y_line = 0.61/(x_line - 0.47) + 1.19

		m = y_line < 1
		plt.plot(x_line[m]-lnd+lbd, y_line[m], 'k')

		## Add MAPPINGS-III grids
		add_grids(ax, '[NI]', '[OIII]')

		fig.savefig('%s/plots/diagnotic.png' % (output))
		plt.close(fig)
	else:
		return_x = return_y * np.nan

# # ------------============== BPT diagram ==================----------
# 	# y is the same as MEx
# 	Prefig()
# 	fig, ax = plt.subplots()

# 	x = np.log10(OIII_to_OI(D.e_line['[OIII]5007d'].flux/1.35)/
# 		D.e_line['Hbeta'].flux/bd)

# 	ax.errorbar(x, y, yerr=y_err)

# 	x_line1 = np.arange(-2.2, 1, 0.001)
# 	y_line1 = 0.73/(x_line1 + 0.59) + 1.33
# 	m = y_line1 < 1
# 	ax.plot(x_line1[m], y_line1[m],'k')

# 	y_line2 = 1.18 * x_line1 + 1.30
# 	m = y_line2 > y_line1
# 	a = np.min(x_line1[m])
# 	x_line2 = np.arange(a, 0.60, 0.001)
# 	y_line2 = 1.18 * x_line2 + 1.30
# 	ax.plot(x_line2, y_line2, 'k')

# 	ax.set_xlim([-2.2, 0])
# 	ax.set_ylim([-1.2, 1.5])

# 	ax.set_ylabel(r'log [OIII]$\lambda$5007/H$\,\beta$')
# 	ax.set_ylabel(r'log [OI]$\lambda$6300/H$\,\alpha$')

# 	fig.savefig('%s/plots/BPT.png' % (output))

# ------------=========== Hbeta flux profile ===============----------
	Prefig()
	fig, ax = plt.subplots()

	Hb = D.e_line['Hbeta'].flux

	radius = np.sqrt((D.xBar - D.center[0])**2 + (D.yBar - D.center[1])**2)
	ax.scatter(radius, Hb,# cmap=sauron, 
		c=np.log10(D.e_line['[OIII]5007d'].flux/D.e_line['Hbeta'].flux).clip(
		-1,1))

	o = np.argsort(radius)

	# ax.plot(radius[o][1:], Hb[np.argmin(np.abs(radius-1))] * radius[o][1:]**-2, 'k')
	ax.plot(radius[o], np.nanmedian(radius[o][np.isfinite(Hb[o])])**2
		* np.nanmedian(Hb[o][np.isfinite(Hb[o])]) * radius[o]**-2, 'k')

	ax.set_yscale('log')
	ax.set_ylabel(r'H$_\beta$ Flux')
	ax.set_xlabel('Radius (pixels)')

	fig.savefig('%s/plots/Hbeta_profile.png' % (output))
	plt.close(fig)
	
# ------------============== WHaN2 diagram ================----------
# Use Hb and NI instead of Ha and NII.
	if all([l in D.e_components for l in ['[NI]d', 'Hbeta']]):
		Prefig(size=(16,12*3), transparent=False)
		fig, ax = plt.subplots(3)
		x = np.log10(D.e_line['[NI]d'].flux/D.e_line['Hbeta'].flux)
		y = np.log10(D.e_line['Hbeta'].equiv_width)

		Ha_Hapeak = D.e_line['Hbeta'].flux/np.nanmax(D.e_line['Hbeta'].flux)
		p1 = Ha_Hapeak <= 0.2

		p2 = Ha_Hapeak > 0.2

		c = np.array(np.sqrt((D.xBar-center[0])**2 +(D.yBar-center[1])**2))

		passive = (D.e_line['Hbeta'].equiv_width<0.5) * (
			D.e_line['[NI]d'].equiv_width<0.5)

		ax[0].scatter(x, y,c=c)
		ax[0].scatter(x[passive],y[passive],marker='x',c='r')


		xlim = ax[0].get_xlim()
		ylim = ax[0].get_ylim()
		xlim = [min(xlim[0], -2+lbd-lnd), max(xlim[1], 0.75+lbd-lnd)]
		ylim = [min(ylim[0], -1-lbd+lsd), max(ylim[1], 1.2-lbd+lsd)]

		ax[1].scatter(x[p1], y[p1],c=c[p1], vmin=np.min(c), vmax=np.max(c))
		ax[1].scatter(x[p1*passive],y[p1*passive],marker='x',c='r')
		ax[1].text(-1.9+lbd-lnd, -0.8-lbd,
			r'0 $<$ H$_\beta$/H$_\beta^{\mathrm{peak}}$ $\leq$ 0.2')

		ax[2].scatter(x[p2], y[p2],c=c[p2], vmin=np.min(c), vmax=np.max(c))
		ax[2].scatter(x[p2*passive],y[p2*passive],marker='x',c='r')
		ax[2].text(-1.9+lbd-lnd, -0.8-lbd, 
			r'0.2 $<$ H$_\beta$/H$_\beta^{\mathrm{peak}}$ $\leq$ 1')
		
		# Dianotics lines from R Cid Fernandes et.al. 2011 MNRAS 413 1687
		for a in ax:
			a.axhline(np.log10(3)-lbd+lsd, ls=':', c='k') # 
			a.plot([-0.4+lbd-lnd, -0.4+lbd-lnd],
				[np.log10(3)-lbd+lsd, ylim[1]], ls=':', c='k')
			a.plot([-0.4+lbd-lnd, xlim[1]], 
				[np.log10(6)-lbd+lsd, np.log10(6)-lbd+lsd], ls=':', c='k')
			a.set_xlim(xlim)
			a.set_ylim(ylim)
			a.set_ylabel(r'log(EW(H$_\beta$)/$\AA$)')
			a.set_xlabel(r'log [NI]$\lambda\lambda$5197,5200/H$\,\beta$')
			a.text(-1.9+lbd-lnd,1-lbd+lsd,'Star Forming')
			a.text(-0.35+lbd-lnd,1-lbd+lsd,'strong AGN')
			a.text(-0.35+lbd-lnd,0.55-lbd+lsd,'weak AGN')
			a.text(-1.9+lbd-lnd,0.25-lbd+lsd,'Retired Galaxies')
		# fig.suptitle('WHbN1: %s' % (galaxy))

		fig.savefig('%s/plots/WHbN1.png' % (output))
		plt.close(fig)






	if return_sauron_diagnotics:
		return return_x, return_y
	else:
		return D




if __name__=='__main__':
	galaxies = ['eso443-g024', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 
		'ngc1399', 'ngc3100', 'ngc3557', 'ngc7075', 'pks0718-34']
	# galaxies = ['ic1459']
	# fig, ax = plt.subplots()
	# for g in galaxies:
	# 	x, y = BPT(g, return_sauron_diagnotics=True)

	# 	if not all(np.isnan(x)) and not all(np.isnan(y)):
	# 		ax.scatter(x, y, label=g)
	# ax.legend(facecolor='w')

	# add_grids(ax)

	# ax.set_xlim([-2.5, 0.5])
	# ax.set_ylim([-1.5, 1.5])
	# ax.set_xlabel(r'log [NI]d/H$_\beta$')
	# ax.set_ylabel(r'log [OIII]d/H$_\beta$')
	# ax.set_title('SAURON diagnotics for all sample')
	# fig.savefig('%s/Data/vimos/analysis/diagnotic.png' % (cc.base_dir))

	BPT('ic1459')
