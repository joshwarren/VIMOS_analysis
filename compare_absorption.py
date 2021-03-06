## ==================================================================
## 		Find absorption line values for given apature
## ==================================================================
## warrenj 20170206 A routine to return absorption line strengths for a 
##	given appature
# import matplotlib # 20160202 JP to stop lack-of X-windows error
# matplotlib.use('Agg')
import numpy as np
from checkcomp import checkcomp
cc = checkcomp()
from errors2 import determine_goodpixels, apply_range, get_stellar_templates, \
	get_emission_templates
import ppxf_util as util
from ppxf import ppxf
from scipy import ndimage # for gaussian blur
from absorption import absorption
from classify import get_R_e
import matplotlib.pyplot as plt 
from label_line import *
from prefig import Prefig 
Prefig()


c = 299792.458 # km/s

stellar_moments = 4 
gas_moments = 2
quiet = True
gas = 3


class compare_absorption(object):
	def __init__(self, galaxy='ngc3557', slit_h=4.5, slit_w=2, slit_pa=30, 
		method='Rampazzo_aperture', r1=0, r2=None, debug=False):
		print galaxy
		if 'Rampazzo' in method:
			from Rampazzo import rampazzo
			# Default is nuclear region: 0<r<R_e/16
			if r2 is None:
				r2 = get_R_e(galaxy)/16
			data = rampazzo(galaxy, method='aperture', slit_h=slit_h, r1=r1, r2=r2, 
				debug=debug)

			if 'gradient' in method:
				if '2' in method: data.method = 'gradient2'
				else: data.method = 'gradient'

		elif method == 'Ogando':
			from Ogando import ogando
			data = ogando(galaxy, debug=debug, slit_h=slit_h, slit_w=slit_w, 
				slit_pa=slit_pa)

		elif method == 'Miles':
			from Miles import miles
			data = miles(galaxy)

		gal_spec, gal_noise = data.get_spec()

		gal_spec, lam, cut = apply_range(gal_spec, window=201, repeats=3, 
			lam=data.lam, set_range=np.array([4200,10000]), return_cuts=True)
		gal_noise = gal_noise[cut]
		lamRange = np.array([lam[0],lam[-1]])/(1+data.z)
## ----------================= Templates ====================---------
		FWHM_gal = 2.5 # VIMOS documentation (and fits header)
		FWHM_gal = FWHM_gal/(1+data.z) # Adjust resolution in Angstrom

		stellar_templates = get_stellar_templates(galaxy, FWHM_gal)
		velscale = stellar_templates.velscale

		e_templates = get_emission_templates(gas, lamRange, 
			stellar_templates.logLam_template, FWHM_gal)

		if gas:
			templates = np.column_stack((stellar_templates.templates, e_templates.templates))
		else:
			templates = stellar_templates.templates
		component = [0]*len(stellar_templates.templatesToUse) + e_templates.component
		templatesToUse = np.append(stellar_templates.templatesToUse, 
			e_templates.templatesToUse)
		element = ['stellar'] + e_templates.element

		start = [[data.vel, data.sig]] * (max(component) + 1)
		moments = [stellar_moments] + [gas_moments] * max(component)
## ----------============== Final calibrations ==============---------
		## smooth spectrum to fit with templates resolution
		if FWHM_gal < stellar_templates.FWHM_tem:
			sigma = stellar_templates.FWHM_dif/2.355/data.f[0].header['CDELT3']
			gal_spec = ndimage.gaussian_filter1d(gal_spec, sigma)
			gal_noise = np.sqrt(ndimage.gaussian_filter1d(gal_noise**2, sigma))
		
		## rebin spectrum logarthmically
		gal_spec_log, logLam_bin, _ = util.log_rebin(lamRange, gal_spec, velscale=velscale)
		gal_noise_log, logLam_bin, _ = util.log_rebin(lamRange, gal_noise**2, 
			velscale=velscale)
		gal_noise_log = np.sqrt(gal_noise_log)

		gal_noise_log = gal_noise_log + 0.0000000000001



		dv = (stellar_templates.logLam_template[0]-logLam_bin[0])*c # km/s
		# Find the pixels to ignore to avoid being distracted by gas emission
		#; lines or atmospheric absorbsion line.  
		goodPixels = determine_goodpixels(logLam_bin,stellar_templates.lamRange_template,
			data.vel, data.z, gas=gas!=0) 
		lambdaq = np.exp(logLam_bin)
## ----------=================== pPXF =======================---------
		pp = ppxf(templates, gal_spec_log, gal_noise_log, velscale, start, 
			goodpixels=goodPixels, moments=moments, degree=-1, vsyst=dv, 
			component=component, lam=lambdaq, plot=not quiet, quiet=quiet, mdegree=10)
		self.pp = pp

		# Only use stellar templates (Remove emission lines)
		stellar_spec = gal_spec_log - \
			pp.matrix[:, -e_templates.ntemp:].dot(pp.weights[-e_templates.ntemp:])
		conv_spec = pp.matrix[:, :-e_templates.ntemp].dot(pp.weights[:-e_templates.ntemp])

		# Generate the unconvolved spectra ('0 km/s' resolution)
		unc_lam = stellar_templates.wav
		unc_spec = stellar_templates.lin_templates.dot(pp.weights[:stellar_templates.ntemp])
		unc_spec *= np.polynomial.legendre.legval(np.linspace(-1,1,
			len(unc_spec)), np.append(1, pp.mpolyweights))
## ----------============== Absorption Line =================---------

		lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
		self.result = {}
		self.uncert = {}
		for line in lines:		
			ab, ab_uncert = absorption(line, lambdaq, stellar_spec, noise=gal_noise_log, 
				unc_lam=unc_lam, unc_spec=unc_spec,	conv_spec=conv_spec)#, 
				#lick=True)

			if method == 'Ogando':
				# Aperture correction
				ab_file = '%s/Documents/useful_files/ab_linelist.dat' % (cc.home_dir)
				i1, i2, b1, b2, r1, r2, units = np.genfromtxt(ab_file, unpack=True, 
					usecols=(1,2,3,4,5,6,7), skip_header=2, skip_footer=2)
				ls = np.genfromtxt(ab_file, unpack=True, dtype=str, usecols=(8), 
					skip_header=2, skip_footer=2)
				l = np.where(ls == line)[0][0]
				index = [i1[l], i2[l]]

				# Convert to mag
				ab = -2.5 * np.log(1 - ab/(index[1]-index[0]))

				# Aperture Correction: beta values taken from paper
				beta = {'H_beta':0.002, 'Fe5015':-0.012, 'Mg_b':-0.031} 
				# If line is not observed by Ogando
				if line not in beta.keys():
					self.result[line] = np.nan
					self.uncert[line] = np.nan
					continue
				H = 70.0 # value used by Bolonga group.
				r_ab = np.degrees(1.19/1000 * H/(c * data.z)) * 60 * 60 # 1.19 kpc -> arcsec
				# Correction is def in eq (9) with r_ab and r_norm def in eq (1) - not 
				#	100% sure if I've got them correct.
				ab = ab - beta[line] * np.log(1.025 * np.sqrt(slit_w*r_ab/np.pi)/slit_h)

				# Back to Angstroms
				ab = (index[1]-index[0]) * (1 - np.exp(ab/-2.5))

			elif 'Rampazzo' in method:
				self.r = data.r
			self.result[line] = ab
			self.uncert[line] = ab_uncert

		if debug:
			for line in lines:
				print '%s:	%.3f +/- %.3f' % (line, self.result[line], self.uncert[line])









def run(galaxy='ic1459', method=None, debug=False):
	if 'Rampazzo' in method:
		pa = {'ic1459':40, 'ic4296':40, 'ngc3557':30}

		# lines measured by Rampazzo
		lines = ['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 'Mg_b']
		f, ax = plt.subplots(len(lines), sharex=True)
		result = {}
		uncert = {}
		for line in lines:
			result[line] = []
			uncert[line] = []
		r = []
		R_e = get_R_e(galaxy)
		
		if 'aperture' in method:
			h = [1.5,2.5,10.0, R_e/10,R_e/8,R_e/4,R_e/2]
			for i in h:
				print 'Aperture: %.3f arcsec' %(i)
				comp_ab = compare_absorption(galaxy, slit_h=R_e, slit_w=2, 
					slit_pa=pa[galaxy], method=method, r1=0, r2 = i, debug=debug)
				for line in lines:
					result[line].append(comp_ab.result[line])
					uncert[line].append(comp_ab.uncert[line])
				r.append(comp_ab.r)
		else: 
			h = [0,R_e/16, R_e/8, R_e/4, R_e/2]
			for i in xrange(len(h)-1):
				print 'Gradient between %.3f and %.3f' % (h[i], h[i+1])
				comp_ab = compare_absorption(galaxy, slit_h=R_e, slit_w=2, 
					slit_pa=pa[galaxy], method=method, r1=h[i], r2 = h[i+1], debug=debug)
				for line in lines:
					result[line].append(comp_ab.result[line])
					uncert[line].append(comp_ab.uncert[line])
				r.append(comp_ab.r)
		
		
		r = np.array(r)
		s = np.argsort(r)
		for line in lines:
			result[line] = np.array(result[line])
			uncert[line] = np.array(uncert[line])

		if 'aperture' in method:
			R_file = '%s/Data/lit_absorption/Rampazzo_aperture.txt' % (cc.base_dir)
		elif 'gradient' in method:
			R_file = '%s/Data/lit_absorption/Rampazzo_gradient.txt' % (cc.base_dir)
		gals = np.loadtxt(R_file, unpack=True, usecols=(0,), skiprows=1, dtype=str)
		i_gal = np.where(gals == galaxy)[0]
		lit_r, G4300, Fe4383, Ca4455, Fe4531, Fe4668, H_beta, Fe5015, Mg_b, e_G4300, \
			e_Fe4383, e_Ca4455, e_Fe4531, e_Fe4668, e_H_beta, e_Fe5015, e_Mg_b = \
			np.genfromtxt(R_file, unpack=True, 
			usecols=(4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21), 
			skip_header=min(i_gal)+1, max_rows=len(i_gal))

		lit_result = {'G4300':G4300, 'Fe4383':Fe4383, 'Ca4455':Ca4455, 'Fe4531':Fe4531, 
			'Fe4668':Fe4668, 'H_beta':H_beta, 'Fe5015':Fe5015, 'Mg_b':Mg_b}
		lit_uncert = {'G4300':e_G4300, 'Fe4383':e_Fe4383, 'Ca4455':e_Ca4455, 
			'Fe4531':e_Fe4531, 'Fe4668':e_Fe4668, 'H_beta':e_H_beta, 'Fe5015':e_Fe5015, 
			'Mg_b':e_Mg_b}
		lit_s = np.argsort(lit_r)

		# Used to covert literture results back to flux calibrated from LICK 'calibrated'
		lit_lick_corr_file = '%s/Data/lit_absorption/Rampazzo_lick_corrections.txt' % (
			cc.base_dir)
		lit_lick_corr_lines = np.loadtxt(lit_lick_corr_file, unpack=True, usecols=(0,), 
			skiprows=1, dtype=str)
		lit_lick_corr_alpha, lit_lick_corr_beta = np.loadtxt(lit_lick_corr_file, 
			unpack=True, usecols=(1,2), skiprows=1)
		from compare_absorption2 import Lick_to_LIS

		for i, line in enumerate(lines):
			# Plot my results
			# Waiting for clarification regarding <r_l> in paper. Use their values in 
			# meantime
			# ax[i].errorbar(r[s]/R_e, result[line][s], yerr=uncert[line][s], color='k')
			ax[i].errorbar(lit_r[lit_s], result[line][lit_s], yerr=uncert[line][lit_s], 
				color='k')
			# Plot Rampazzo results
			i_line = np.where(lit_lick_corr_lines == line)[0][0]
			# ax[i].errorbar(lit_r[lit_s], (lit_result[line][lit_s] - 
			# 	lit_lick_corr_beta[i_line])/lit_lick_corr_alpha[i_line], 
			# 	yerr=lit_uncert[line][lit_s], fmt='--', color='m')
			ax[i].errorbar(lit_r[lit_s], Lick_to_LIS(line, lit_result[line][lit_s]), 
				yerr=lit_uncert[line][lit_s], fmt='--', color='m')

			# Remove underscores as matplotlib is currently automatically calling latex...
			if line == 'H_beta':
				line = 'Hbeta'
			elif line =='Mg_b':
				line = 'Mgb'
			ax[i].set_title(line)

		f.savefig('%s/Data/lit_absorption/%s_%s.png' % (cc.base_dir, method, galaxy), 
			bbox_inches="tight")

		# plt.show()
	elif method == 'Ogando':
		####*********CHECK PAs*****************

		## NEED TO ADD PA TO KINEMATICS.PY - CURRENTLY ONLY SAVE THE MISALIGNMENTS


		gals = ['ic1459', 'ic4296','ngc1399','ngc3100','ngc3557','ngc7075']
		fig, ax = plt.subplots()

		for i, galaxy in enumerate(gals):
			comp_ab = compare_absorption(galaxy, slit_h=4.1, slit_w=2.5, slit_pa=30,
				method=method, debug=debug)
			data_file = '%s/Data/lit_absorption/Ogando.txt' % (cc.base_dir)
			galaxy_gals = np.loadtxt(data_file, skiprows=2, usecols=(0,),dtype=str,
				unpack=True)
			i_gal = np.where(galaxy_gals==galaxy)[0][0]
			Ogando_sig, Ogando_sig_unc, H_beta, H_beta_uncert, Fe5015, Fe5015_uncert, \
				Mg_b, Mg_b_uncert = np.loadtxt(data_file, skiprows=2, unpack=True, 
				usecols=(1,2,4,5,6,7,8,9))

			Ogando_obs = {'H_beta':H_beta[i_gal], 'Fe5015':Fe5015[i_gal], 
				'Mg_b':Mg_b[i_gal]}
			Ogando_unc = {'H_beta':H_beta_uncert[i_gal], 'Fe5015':Fe5015_uncert[i_gal], 
				'Mg_b':Mg_b_uncert[i_gal]}
			c = {'H_beta':'r', 'Fe5015':'b', 'Mg_b':'g'}
			for line in Ogando_obs.keys():
				ax.errorbar(np.log10(comp_ab.pp.sol[0][1]), comp_ab.result[line], 
					yerr=comp_ab.uncert[line], fmt='x', color=c[line])
				ax.errorbar(Ogando_sig[i_gal], Ogando_obs[line], 
					xerr=Ogando_sig_unc[i_gal], yerr=Ogando_unc[line], fmt='o', 
					color=c[line])
				l = ax.plot([np.log10(comp_ab.pp.sol[0][1]), Ogando_sig[i_gal]],
					[comp_ab.result[line], Ogando_obs[line]], '--', c=c[line])
				if line == 'H_beta':
					ax.annotate(galaxy, xy=(Ogando_sig[i_gal], Ogando_obs[line]), 
						textcoords='data')

			ax.legend(c)
			ax.set_xlabel('log(sigma (km/s))')
			ax.set_ylabel('Index Strength (A)')
		fig.savefig('/Data/lit_absorption/Ogando.png', bbox_inches='tight')
		# plt.show()

	elif method == 'Miles':
		fig, ax = plt.subplots()

		gals = ['ic1459', 'ic1531', 'ic4296','ngc0612', 'ngc1399','ngc3100','ngc3557',
			'ngc7075', 'pks0718-34', 'eso443-g024']
		lines = np.array(['G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'H_beta', 'Fe5015', 
			'Mg_b'])
		line_w = [4300, 4383, 4455,4531,4861,5015,5170]
		c = {'ic1459':'saddlebrown', 'ic1531':'k', 'ic4296':'g','ngc0612':'y', 'ngc1399':'m',
			'ngc3100':'r','ngc3557':'c','ngc7075':'orange', 'pks0718-34':'darkgreen', 
			'eso443-g024':'gray'}


		models_dir  = '%s/models/TMJ_SSPs/tmj.dat' % (cc.home_dir)

		for i, galaxy in enumerate(gals):
			comp_ab = compare_absorption(galaxy=galaxy, method=method)
			for j, l in enumerate(lines):
				ax.errorbar(line_w[j], comp_ab.result[l], yerr=comp_ab.uncert[l],
					color=c[galaxy], zorder=10, fmt='x')

		titles = np.loadtxt(models_dir, dtype=str)[0]
		import matplotlib.patches as patches

		for i, l in enumerate(titles):
			if l in lines:
				TMJ_index = np.loadtxt(models_dir, usecols=(i,), unpack=True, 
					skiprows=35)
				j = np.where(lines == l)[0][0]
				ax.add_patch(patches.Rectangle((line_w[j] - 20, np.min(TMJ_index)), 40,
					np.max(TMJ_index)-np.min(TMJ_index), zorder=5))

		ax.set_xlabel('Wavelength (A)')
		ax.set_ylabel('Index Indicies')
		ax.set_title('Global index measurements (with TMJ models ranges)')
#		ax.legend(c)

		fig.savefig('/Data/lit_absorption/Miles.png', bbox_inches='tight')


if __name__ == '__main__':
	# galaxy = 'ngc3557'
	galaxy = 'ic1459'
	galaxy = 'ic4296'
	import subprocess
	for galaxy in ['ic1459','ic4296','ngc3557']:

		run(galaxy = galaxy, method = 'Rampazzo_gradient2', debug=True)
		# run(galaxy = galaxy, method = 'Rampazzo_gradient2', debug=True)
	# 	subprocess.call(['/bin/bash', '-i', '-c', "push 'done %s Rampazzo'" % (galaxy)])

	# run(method = 'Ogando', debug=False)
	# subprocess.call(['/bin/bash', '-i', '-c', "push 'done Ogando comparison'"])
	# run(method = 'Miles', debug=False)
	# subprocess.call(['/bin/bash', '-i', '-c', "push 'done Miles comparison'"])