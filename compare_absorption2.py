# New routine to compare VIMOS observations with Rampazzo

import numpy as np 
from checkcomp import checkcomp
cc=checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
from prefig import Prefig
Prefig()
import cPickle as pickle
from astropy.io import fits
from classify import get_R_e
from errors2 import apply_range, get_dataCubeDirectory,run_ppxf, set_params
from glob import glob
from pop import get_absorption
from ppxf import create_plot
from spectools import spectrum
from scipy.interpolate import interp1d


c = 299792.458 # speed of light in km/s

def Lick_to_LIS(line, value, res=8.4):
	if line =='H_beta' or line=='Hb' or line=='hb':
		line='Hbeta'
	if line =='Mgb' or line == 'mgb' or line =='mg_b':
		line='Mg_b'


	file = '%s/Data/lit_absorption/Miles/lick_to_lis.txt' % (cc.base_dir)
	if res==5.0:
		a0, a1, a2, a3 = np.loadtxt(file, unpack=True, usecols=(1,2,3,4), skiprows=1)
	elif res==8.4:
		a0, a1, a2, a3 = np.loadtxt(file, unpack=True, usecols=(5,6,7,8), skiprows=1)
	elif res==14.0:
		a0, a1, a2, a3 = np.loadtxt(file, unpack=True, usecols=(9,10,11,12), 
			skiprows=1)
	line_name = np.loadtxt(file, unpack=True, usecols=(0), dtype=str, skiprows=1)

	i_line = np.where(line_name == line)[0][0]

	return a0[i_line] + a1[i_line]*value + a2[i_line]*value**2 + a3[i_line]*value**3

def compare_absortion(galaxy, R_sig=False, corr_lines='all'):
	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[0].header

	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
	color = ['r',        'b',       'g',     'c',    'purple',   'k',  'orange']



	R_e = get_R_e(galaxy)
	apertures = np.array([1.5, 2.5, 10, R_e/10, R_e/8, R_e/4, R_e/2]) # arcsec

	Ramp_sigma = {'ngc3557':[265, 247, 220], 'ic1459':[311, 269, 269], 
		'ic4296':[340, 310, 320]}

	R_sigma  = interp1d([R_e/8, R_e/4, R_e/2], Ramp_sigma[galaxy], 
		fill_value=(Ramp_sigma[galaxy][0], Ramp_sigma[galaxy][2]), 
		bounds_error=False)

	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals, x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,4,5), dtype='float,int,int')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	center = np.array([x_cent_gals[i_gal], y_cent_gals[i_gal]])

	index = np.zeros((40,40,2))
	for i in range(40):
		for j in range(40):
			index[i,j,:] = np.array([i,j]) - center

	fig, ax = plt.subplots()
	fig3, ax3 = plt.subplots()
	fig4, ax4 = plt.subplots()
	ax5 = ax4.twinx()
	my_values = {}
	my_errors = {}
	sigma = np.array([])
	t=[]
	e=[]
	g=[]
	h=[]
	j=[]
	r=[]
	w=[]
	for a in apertures:
		params = set_params(reps=0, opt='pop', gas=1, lines=corr_lines, 
			produce_plot=False)

		mask = np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) * header['CDELT3'] < a

		spec = np.nansum(f[0].data[:,mask], axis=1)
		noise = np.sqrt(np.nansum(f[1].data[:,mask]**2, axis=1))

		lam = np.arange(len(spec))*header['CDELT3'] + header['CRVAL3']
		spec, lam, cut = apply_range(spec, lam=lam, set_range=params.set_range, 
			return_cuts=True)
		lamRange = np.array([lam[0],lam[-1]])
		noise = noise[cut]

		pp = run_ppxf(galaxy, spec, noise, lamRange, header['CDELT3'], params)

		plot = create_plot(pp)
		plot.lam = pp.lam*(1+pp.z)/(1+pp.z+(pp.sol[0][0]/c))
		fig2, ax2, = plot.produce
		s = spectrum(lam=pp.lam, lamspec=pp.galaxy)
		for i, l in enumerate(lines):
		  if l=='H_beta' or l=='Hbeta':
			  l='hb'
		  elif l=='Mg_b':
			  l='mgb'
		  elif l=='NaD':
			  l='nad'
		  elif l=='TiO1':
			  l='tio1'
		  elif l=='TiO2':
			  l='tio2'
		  elif l=='Fe5270':
			  l='fe52'
		  elif l=='Fe5335':
			  l='fe53'
		  ax2.axvspan(*getattr(s,l),color='b', alpha=0.5)
		  ax2.axvspan(getattr(s,l+'cont')[0], getattr(s,l+'cont')[1], color='r', 
			  alpha=0.5)
		  ax2.axvspan(getattr(s,l+'cont')[2], getattr(s,l+'cont')[3], color='r', 
			  alpha=0.5)
		  lims = ax2.get_ylim()
		  if i%2==0:
			  ax2.text(np.mean(getattr(s,l)), lims[1] - 0.1*(lims[1]-lims[0]), l, 
				  size=8, ha='center')
		  else:
			  ax2.text(np.mean(getattr(s,l)), lims[1] - 0.15*(lims[1]-lims[0]), l, 
				  size=8, ha='center')

		fig2.savefig('/Data/lit_absorption/Rampazzo/%s_rad_%.2f.png'%(galaxy, a))
		plt.close(fig2)

		if R_sig:
			if isinstance(R_sig, bool):
				absorp, uncert = get_absorption(lines, pp=pp, sigma=R_sigma(a),
					instrument='vimos')
				sigma = np.append(sigma, R_sigma(a))
			else:
				absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos',
					sigma=R_sigma(a)+R_sig*(pp.sol[0][1]-R_sigma(a)))
				sigma = np.append(sigma, R_sigma(a)+R_sig*(pp.sol[0][1]-R_sigma(a)))
		else:
			absorp, uncert = get_absorption(lines, pp=pp, instrument='vimos')
			sigma = np.append(sigma, pp.sol[0][1])

		for l in lines:
			if a == min(apertures):
				my_values[l] = np.array([])
				my_errors[l] = np.array([])
			my_values[l] = np.append(my_values[l], absorp[l])
			my_errors[l] = np.append(my_errors[l], uncert[l])
		
		for i, l in enumerate(lines):
			ax.errorbar(a, absorp[l], yerr=uncert[l], color=color[i], fmt='x')
	for i, l in enumerate(lines):
	  ax.errorbar(np.nan, np.nan, color=color[i], fmt='x', label=l)
	ax.legend(facecolor='w')


	Rampazzo_file = '%s/Data/lit_absorption/Rampazzo_aperture.txt' % (cc.base_dir)
	file_headings = np.loadtxt(Rampazzo_file, dtype=str)[0]

	for i, l in enumerate(lines):
		col = np.where(file_headings==l)[0][0]
		try:
			col2 = np.where(file_headings==l)[0][1]
		except:
			try:
				col2 = np.where(file_headings=='_'+l)[0][0]
			except:
				col2 = np.where(file_headings=='e_'+l)[0][0]
		R_obs, R_err = np.loadtxt(Rampazzo_file, unpack=True, skiprows=1, 
			usecols=(col,col2))
		R_galaxies = np.loadtxt(Rampazzo_file, unpack=True, skiprows=1, usecols=(0,), 
			dtype=str)

		mask = R_galaxies==galaxy

		order = np.argsort(apertures)

		lit_value = Lick_to_LIS(l, R_obs[mask][order])
		err = np.mean([np.abs(Lick_to_LIS(l, R_obs[mask][order] + 
			R_err[mask][order]) - Lick_to_LIS(l, R_obs[mask][order])), 
			np.abs(Lick_to_LIS(l, R_obs[mask][order] - R_err[mask][order]) -
			Lick_to_LIS(l, R_obs[mask][order]))], axis=0) 

		ax.errorbar(apertures[order], lit_value, yerr=err, color=color[i])

		if l=='H_beta' or l=='Hbeta':
			l2='hb'
		elif l=='Mg_b':
			l2='mgb'
		elif l=='NaD':
			l2='nad'
		elif l=='TiO1':
			l2='tio1'
		elif l=='TiO2':
			l2='tio2'
		elif l=='Fe5270':
			l2='fe52'
		elif l=='Fe5335':
			l2='fe53'
		else:
			l2=l

		ax3.scatter(
		  np.abs(my_values[l][order] - lit_value)/my_values[l][order],
		  np.abs(my_values[l][order] - lit_value)/np.sqrt(err**2 + 
		  my_errors[l][order]**2), color=color[i], s=4*apertures[order]**2,
		  label=l)
		ax4.scatter(sigma[order], 
		  np.abs(my_values[l][order] - lit_value)/my_values[l][order],
		  color=color[i], s=4*apertures[order]**2, label=l)
		ax5.scatter(sigma[order], 
		  np.abs(my_values[l][order] - lit_value)/np.sqrt(err**2 + 
		  my_errors[l][order]**2), marker='x',
		  color=color[i], s=4*apertures[order]**2, label=l)
		t.extend(sigma[order])
		g.extend(my_values[l][order]) 
		h.extend(lit_value)
		j.extend(err)
		e.extend(my_errors[l][order])
		r.extend([lines[i]]*len(sigma))
		w.extend(4*apertures[order]**2)


	ax.set_ylabel(r'Index strength, $\AA$')
	ax.set_xlabel('Radius, arcsec')
	fig.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_%i.png' % (
	  cc.base_dir, galaxy, params.gas))
	plt.close(fig)


	ax3.set_ylabel(r'Sigma difference ((Mine - Ramp)/Combined Uncert)')
	ax3.set_xlabel('Fractional difference ((Mine - Ramp)/Mine)')
	ax3.legend()
	fig3.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_fractional.png' % (
	  cc.base_dir, galaxy))
	plt.close(fig3)


	ax4.set_xlabel('Vel dispersion')
	ax3.set_ylabel('Fractional difference ((Mine - Ramp)/Mine)')
	ax5.set_ylabel(r'Sigma difference ((Mine - Ramp)/Combined Uncert)')
	ax4.legend()
	fig4.savefig('%s/Data/lit_absorption/Rampazzo_aperture_%s_sigma.png' % (
	  cc.base_dir, galaxy))
	plt.close(fig4)


	return t, g, h, j, e, r, w

if __name__=='__main__':
	lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
	color = ['r',        'b',       'g',     'c',    'purple',   'k',  'orange']

	fig, ax = plt.subplots()
	ax2 = ax.twinx()
	s2 = [] 
	f2 = []
	g2 = []
	h2 = []
	j2 = []
	r2 = []
	z2 = []
	p2 = []
	for gal in ['ic1459','ic4296','ngc3557']:
		s, f, g, h, j, r, z = compare_absortion(gal, 
			corr_lines = ['Hbeta','[OIII]5007d'])
		s2.extend(s) # sigma
		f2.extend(f) # my_values
		g2.extend(g) # Ramp
		h2.extend(h) # Ramp err
		j2.extend(j) # my_error
		r2.extend(r) # index
		z2.extend(z) # apature size
		p2.extend([gal]*len(s)) # galaxy name



	s2=np.array(s2)
	f2=np.array(f2)
	g2=np.array(g2)
	h2=np.array(h2)
	j2=np.array(j2)
	r2=np.array(r2)
	z2=np.array(z2)
	p2=np.array(p2)

	fig, ax=plt.subplots(subplot_kw={'aspect':'equal'})

	for i in range(len(np.unique(r2))):
		for gal in ['ic1459','ic4296','ngc3557']:
			m = (r2 == np.unique(r2)[i]) * (p2 == gal)
			params = np.polyfit(f2[m], g2[m], 1, 
				w=1/np.sqrt(j2[m]**2 + h2[m]**2))
			i_gal = np.where(np.array(lines) == np.unique(r2)[i])[0][0]
			if gal == 'ic1459':
				ax.errorbar(f2[m], g2[m], xerr=j2[m], yerr=h2[m], fmt='x', 
					label=np.unique(r2)[i], color=color[i_gal])
			else:
				ax.errorbar(f2[m], g2[m], xerr=j2[m], yerr=h2[m], fmt='o', 
					color=color[i_gal])
	for i in range(len(np.unique(r2))):
			m = (r2 == np.unique(r2)[i])
			params = np.polyfit(f2[m], g2[m], 1, 
				w=1/np.sqrt(j2[m]**2 + h2[m]**2))
			print np.unique(r2)[i], 'Offset:', np.mean(f2[m] - g2[m]), \
				'dispersion:', np.std(f2[m] - g2[m])

	ax.plot([1,7],[1,7],'k')
	ax.legend(facecolor='w')
	ax.set_xlabel(r'Mine values $\AA$')
	ax.set_ylabel(r'Rampazzo values $\AA$')
	fig.savefig('%s/Data/lit_absorption/comparison_to_Rampazzo_vimos.png' % (
		cc.base_dir))

