## Making use of Sam Vaughan's Thomas models code

from Thomas_models_index_model_variation import *
import numpy as np 
import matplotlib.pyplot as plt 
import cPickle as pickle
from checkcomp import checkcomp
cc=checkcomp()
from spectools import *
cm = plt.get_cmap('plasma') 

def abs_line_plots(galaxy, i1, i2):

	wav_range = '4200-'
	wav_range_dir = '4200-/'
	out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_pickle = '%s/pickled' % (output)
	pickleFile = open("%s/dataObj_%s_pop.pkl" % (out_pickle, wav_range), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	spectra=[]
	variance=[]
	index1_ob=[]
	index2_ob=[]
	index1_va=[]
	index2_va=[]

	for i in range(D.number_of_bins):
		spectra.append(spectrum(D.bin[i].lam, D.bin[i].continuum))
		variance.append(spectrum(D.bin[i].lam, D.bin[i].noise))
		index_value, index_va, continuum_def,feature_def = spectra[i].irindex(0.71,i1,varSED=variance[i])
		index1_ob.append(index_value)
		index1_va.append(index_va)
		index_value, index_va, continuum_def,feature_def = spectra[i].irindex(0.71,i2,varSED=variance[i])
		index2_ob.append(index_value)
		index2_va.append(index_va)

	index1_ob=np.array(index1_ob)
	index2_ob=np.array(index2_ob)
	index1_va=np.array(index1_va)
	index2_va=np.array(index2_va)
	r = np.sqrt((np.array(D.xBar)-16)**2 + (np.array(D.yBar)-16)**2)
	#r=np.vstack([r,r,r,r])


	f,ax = plt.subplots()
	gal_spec = spectrum(D.bin[5].lam,D.galaxy_continuum)
	gal_var = spectrum(D.bin[5].lam,D.galaxy_varience)
	index1_value_best, index1_va_best, continuum_def,feature_def = gal_spec.irindex(0.71,i1,varSED=gal_var)
	index2_value_best, index2_va_best, continuum_def,feature_def = gal_spec.irindex(0.71,i2,varSED=gal_var)

	#ax.errorbar(index1_value_best,index2_value_best, xerr=index1_va_best, yerr=index2_va_best,color='r')
	ax.scatter(index1_value_best,index2_value_best, color='black',marker='*',zorder=12)

	index1, index1_name = get_index(i1)
	index2, index2_name = get_index(i2)
	index1=np.swapaxes(np.transpose(index1.reshape(4, 6, 20), (0, 2, 1)), 1, 2)
	index2=np.swapaxes(np.transpose(index2.reshape(4, 6, 20), (0, 2, 1)), 1, 2)

	plot_index_vs_index(ax, index1, index1_name, index2, index2_name)
	c=np.array([cm(1.0*a/D.number_of_bins) for a in np.argsort(r)])#,axis=1)###################<-----------not working
	mask = (index1_va<index1_ob) + (index2_va<index2_ob)
	for i in range(D.number_of_bins):
		ax.scatter(index1_ob[i],index2_ob[i],2,c=c[i])
		#if index2_va[i] < 0.2*index2_ob[i] and index1_va[i] < 0.2*index1_ob[i]:
		#	ax.errorbar(index1_ob[i],index2_ob[i],xerr=index1_va[i],yerr=index2_va[i],c=c[i])


	a = np.where((index1_va.flatten()/index1_ob.flatten() < 0.2)*(index1_ob.flatten() > 0))[0]
	b = np.where((index2_va.flatten()/index2_ob.flatten() < 0.2)*(index2_ob.flatten() > 0))[0]

	# for i in a:
	# 	ax.errorbar(index1_ob[i],index2_ob[i],xerr=index1_va[i],yerr=index2_va[i],c=c[i])
	# for i in b:
	# 	ax.errorbar(index1_ob[i],index2_ob[i],xerr=index1_va[i],yerr=index2_va[i],c=c[i])

	legend2=ax.legend(fancybox=True, fontsize=10, loc="best")


if __name__=="__main__":

	galaxy = 'ngc3100'
	i1 = 'Hb'
	i2 = 'Mgb'
	i3='Fe5015'
	absorption(galaxy, i1, i2)
	absorption(galaxy, i1, i3)
	absorption(galaxy, i2, i3)

	plt.show()
