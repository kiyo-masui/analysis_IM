#! /usr/bin/env python

import scipy as sp
import numpy as np
#from numpy.fft import *
import scipy.linalg as linalg

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import MakePower
import matplotlib.pyplot as plt
import fftw3 as FFTW
from scipy.interpolate import interp1d

#dirty_map = algebra.load("15hr_41_dirty_map_I.npy")
#dirty_map = algebra.make_vect(dirty_map)

#print dirty_map.mean()
#print dirty_map.axes
#
#fq = dirty_map.get_axis('freq')/1.0e6
#print fq.size, fq[0]

pi = 3.1415926
deg2rad = pi/180.

params_init = {
	# the process number used in this module, now we just support processes=1
	'processes' : 1,

	# indicated the location of the data to be plot
	'input_root' : '../newmaps/',
	# the dirction where the figure will be saved in
	'output_root' : './',

	# parameters that used to make up the file name
	# resultf == None:
	# 		In this case, the result data must be the cross(auto) power spectrum 
	# 		of two(two same) maps. The filename will make up by the parameter
	# 		"hr" and "last".
	# 		resultf = hr[0] + last[0] + '-' + hr[0] + last[0]
	# 		If len(last)==0, it will ignore the 'last[]'
	#
	# 		for power spectrum, the data should be 
	# 		'k_' + params['resultf'] + '.npy'
	#  	'PK_' + params['resultf'] + '.npy'
	# result == 'XXXX'
	# 		In this case, the result data should be the result of more than two
	#		maps. 
	# 		for power spectrum, the data should be 
	# 		'k_combined' + params['resultf'] + '.npy'
	#  	'PK_combined' + params['resultf'] + '.npy'
	'resultf': None,
	'hr' : (),
	'last' : (),

	'jknumber' : 100.,
	# set to import the FKP weighted power spectrum
	'FKPweight' : False,

	# set the parameters for the theoretical power spectrum from CAMB
	'OmegaHI' : 1.e-3,
	'Omegam' : 0.24,
	'OmegaL' : 0.76,
	'z' : 1.,

	# set the lower limit of the y-axis in the figure
	'ymin' : 1.e3,
	# set the upper and lower limit of the x-axis in the figure
	'kmin' : None,
	'kmax' : None,
	# set the unit of the power spectrum
	'PKunit' : 'mk',
	# set to plot the cross power spectrum. 
	# Then the unit will be changed to square root of the auto power spectrum
	'cross' : False,
	# set to calibrate with transfer function
	'biascal' : False,
	# set the location of the transfer function file
	'biascaldir_k' : '', 
	'biascaldir_b' : '',
}
prefix = 'pt_'

class PowerSpectrumPlot(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

	def execute(self, nprocesses=1):
		params = self.params

		resultf = params['resultf']
		if resultf==None:
			resultf = params['hr'][0]
			if len(params['last']) != 0:
				resultf = resultf + params['last'][0]
			resultf = resultf + '-' + params['hr'][1]
			if len(params['last']) != 0:
				resultf = resultf + params['last'][1]

		FKPweight = params['FKPweight']
		kmin = params['kmin']
		kmax = params['kmax']

		#power_th = np.load('../maps/test/power_yuebin.npy')
		#power_th[0] = power_th[0]/0.705
		#power_th[1] = power_th[1]*0.705**3




		# import the power spectrum result
		if params['resultf']==None:
			k = sp.load(params['input_root'] + 'k_'+resultf+'.npy')
			PK = sp.load(params['input_root'] + 'PK_'+resultf+'.npy')
		else:
			k = sp.load(params['input_root'] + 'k_combined_'+resultf+'.npy')
			PK = sp.load(params['input_root'] + 'PK_combined_'+resultf+'.npy')

		PK[PK<0]=0

		non0 = PK.nonzero()

		PK = PK.take(non0)[0]
		k = k.take(non0)[0]

		## Calibrate the bias
		#if params['biascal']:
		#	biask = sp.load(params['biascaldir_k']).take(non0)[0]
		#	bias = sp.load(params['biascaldir_b']).take(non0)[0]
		#	if (biask-k).any():
		#		print "biask and k not match!!"
		#		return 0
		#	#B = interp1d(biask, bias, kind='cubic')
		#	#non0 = PK.nonzero()
		#	#k_biascal = k.take(non0)[0]
		#	print biask
		#	print bias
		#	#if params['cross']:
		#	#	biascal = biascal.__pow__(0.5)
		#	#	print biascal
		#	PK_biascal = PK * bias
		#	k_biascal = biask
		#	Averageplot = True
		#	try:
		#		PKvar = sp.load(params['input_root'] + 
		#			'PKvar_combined_'+resultf+'.npy')
		#		PKerr_average_biascal = np.ndarray(shape=(2,len(non0[0])))
		#		PKerr_average_biascal[0] = PKvar.copy().take(non0)[0] * bias
		#		PKerr_average_biascal[1] = PKvar.copy().take(non0)[0] * bias
		#		print PKerr_average_biascal
		#		for i in range(len(non0[0])):
		#			if PKerr_average_biascal[0][i] >=PK_biascal[i]:
		#				PKerr_average_biascal[0][i] = PK_biascal[i]-1.e-15
		#	except IOError:
		#		print '\t::No Average Error!!'
		#		Averageplot = False

		#	#PK_biascal[PK_biascal<0]=0
		#	#non0 = PK_biascal.nonzero()
		#	#PK_biascal = PK_biascal.take(non0)[0]
		#	#k_biascal = k_biascal.take(non0)[0]
		#	#if Averageplot:
		#	#	PKerr_average_biascal = np.ndarray(shape=(2,len(non0[0])))
		#	#	PKerr_average_biascal[0]= \
		#	#		PKerr_average_biascal_full[0].take(non0)[0]
		#	#	PKerr_average_biascal[1]= \
		#	#		PKerr_average_biascal_full[1].take(non0)[0]
		#	#	for i in range(len(PKerr_average_biascal[0])):
		#	#		if PKerr_average_biascal[0][i] >=PK_biascal[i]:
		#	#			PKerr_average_biascal[0][i] = PK_biascal[i]-1.e-15

#		PKeach = sp.load(params['input_root'] + 'PKeach_' + resultf + '.npy')
#		for ii in range(PKeach.shape[0]):
#			plt.scatter(k, PKeach[ii])
#		plt.loglog()
#		#plt.semilogx()
#		ymin = params['ymin']
#		plt.ylim(ymin=ymin)	
#		plt.xlim(xmin=0.01, xmax=0.9)
#		plt.show()
#
#		return 0

		# import the jk error
		#PKjk = sp.load(params['input_root'] + 'PKjk_'+resultf+'.npy')
		JKplot = True
		try:
			PKvar = sp.load(params['input_root'] + 'PKvar_'+resultf+'.npy')
			PKerr = np.ndarray(shape=(2,len(non0[0])))
			PKerr[0] = PKvar.take(non0)[0]
			PKerr[1] = PKvar.take(non0)[0]

			for i in range(len(PKerr[0])):
				if PKerr[0][i] >=PK[i]:
					PKerr[0][i] = PK[i]-1.e-10
		except IOError:
			print '\t::No Jeckknife Error!!'
			JKplot = False

		# import the average error
		Averageplot = True
		try:
			PKvar = sp.load(params['input_root'] + 
				'PKvar_combined_'+resultf+'.npy')
			PKerr_average = np.ndarray(shape=(2,len(non0[0])))
			PKerr_average[0] = PKvar.take(non0)[0]
			PKerr_average[1] = PKvar.take(non0)[0]

			for i in range(len(PKerr_average[0])):
				if PKerr_average[0][i] >=PK[i]:
					PKerr_average[0][i] = PK[i]-1.e-15
		except IOError:
			print '\t::No Average Error!!'
			Averageplot = False

		# import the fkp result
		if FKPweight:
			kfkp = sp.load(params['input_root'] + \
				'k_fkp_'+resultf+'.npy')
			PKfkp = sp.load(params['input_root'] + \
				'PK_fkp_'+resultf+'.npy')
			PKjkfkp = sp.load(params['input_root'] + \
				'PKjk_fkp_'+resultf+'.npy')
			PKvarfkp = sp.load(params['input_root'] + \
				'PKvar_fkp_'+resultf+'.npy')
			non0fkp = PKfkp.nonzero()
			kfkp = kfkp.take(non0fkp)[0]
			PKfkp = PKfkp.take(non0fkp)[0]
			PKerrfkp = np.ndarray(shape=(2,len(non0fkp[0])))
			PKerrfkp[0] = PKvarfkp.take(non0fkp)[0]
			PKerrfkp[1] = PKvarfkp.take(non0fkp)[0]
			for i in range(len(PKerrfkp[0])):
				if PKerrfkp[0][i]>= PKfkp[i]:
					PKerrfkp[0][i] = PKfkp[i]-1.e-10

		# import the theoretical power spectrum 
		PKcamb = sp.load(params['input_root']+'PKcamb.npy')
		#PKnonl = sp.load(params['input_root']+'nonlPK_'+resultf+'.npy')
		#PKnonl_k = sp.load(params['input_root']+'k_nonlPK_'+resultf+'.npy')
		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		z = params['z']
		a3 = (1+z)**(-3)
		Tb = (OmegaHI) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
			* ((1.+z)/2.5)**0.5
		if params['PKunit']=='mK':
			Tb = Tb/1.e-3

		xx = ((1.0/Omegam)-1.0)/(1.0+z)**3
		num = 1.0 + 1.175*xx + 0.3046*xx**2 + 0.005335*xx**3
		den = 1.0 + 1.875*xx + 1.021 *xx**2 + 0.1530  *xx**3

		G = (1.0 + xx)**0.5/(1.0+z)*num/den
		#print G**2*(1+z)**(-2)

		PKcamb[1] = PKcamb[1]*(G**2*(1+z)**(-2))
		if params['cross']: 
			PKcamb[1] = PKcamb[1]*T
		else: 
			PKcamb[1] = PKcamb[1]*(Tb**2)
			
#		PKnonl = PKnonl*(G**2*(1+z)**(-2))
#		PKnonl = PKnonl*(Tb**2)
#
		#power_th[1] = power_th[1]*(G**2*(1+z)**(-2))
		#power_th[1] = power_th[1]*(Tb**2)

		
		# plot
		plt.figure(figsize=(8,4))
		
		#for i in range(0,50):
		#	plt.scatter(k, PKjk[i].take(non0), alpha=0.5, c='b')
		#for i in range(200,250):
		#	plt.scatter(k, PKjk[i].take(non0), alpha=0.5, c='r')
		#plt.loglog()
		#plt.ylim(ymin=1.e-17)	
		#plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		#plt.show()
		#return 

		plt.subplot('111')
		plt.plot(PKcamb[0], PKcamb[1], 'g-', linewidth=2, 
			label='Theoretical Power Spectrum')
		#plt.plot(PKnonl_k, PKnonl, 'k-', linewidth=2,
		#	label='Theoretical Convolved wiht Window')
		if JKplot:
			plt.errorbar(k, PK, PKerr, fmt='o', c='b', 
				label='Noise Inv Weight', capsize=4.5, elinewidth=1)
		elif Averageplot:
			plt.errorbar(k, PK, PKerr_average, fmt='o', c='r',
				label='Noise Inv Weight', capsize=4.5, elinewidth=1)
			#if params['biascal']:
			#	print PK_biascal
			#	plt.errorbar(k_biascal, PK_biascal, PKerr_average_biascal,
			#		fmt='o', c='r',
			#		#label='Noise Inverse Weighted', 
			#		label='21cm Auto Power Spectrum', 
			#		capsize=4.5, elinewidth=1)
		else:
			plt.scatter(k, PK, c='b', label='Noise Inv Weight')
		if FKPweight:
			plt.errorbar(kfkp, PKfkp, PKerrfkp, fmt='o', c='r',
				label='FPK Weight', capsize=4)
#		for ii in range(PKeach.shape[0]):
#			plt.scatter(k, PKeach[ii].take(non0)[0])
		plt.loglog()
		#plt.semilogx()
		ymin = params['ymin']
		plt.ylim(ymin=ymin)	
		if (kmin==None) or (kmax==None):
			plt.xlim(xmin=0.01, xmax=0.9)
		else:
			plt.xlim(xmin=kmin-0.01*kmin, xmax=kmax+0.01*kmax)

		#plt.xlim(xmin=0.01, xmax=k.max()+0.1*k.max())
		#plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		plt.title('Power Spectrum')
		plt.xlabel('$k (hMpc^{-1})$')
		if params['cross']:
			ylabel = '$P(k) (%(PKunit)s(h^{-1}Mpc)^3)$' %params
		else:
			ylabel = '$P(k) (%(PKunit)s^{2}(h^{-1}Mpc)^3)$' %params
		plt.ylabel(ylabel)
		plt.legend(loc=0, scatterpoints=1)

		plt.savefig(params['output_root']+'power_err_'+resultf+'.eps', 
			format='eps')

		plt.show()

#	def xyz(self, ra, de, r, ra0=0.):
#		x = r*sin(0.5*pi-de)*cos(ra-ra0)
#		y = r*sin(0.5*pi-de)*sin(ra-ra0)
#		z = r*cos(0.5*pi-de)
#		return x, y, z
#
#	def fq2r(self, freq, freq0=1.4e9 , c_H0 = 2.99e3, Omegam=0.27, Omegal=0.73):
#		"""change the freq to distence"""
#		zz =  freq0/freq - 1.
#		for i in range(0, zz.shape[0]):
#			zz[i] = c_H0*self.funcdl(zz[i], Omegam, Omegal)
#		return zz
#	
#	def discrete(self, array):
#		"""discrete the data pixel into small size"""
#		newarray = sp.zeros(self.params['discrete']*(array.shape[0]-1)+1)
#		for i in range(0, array.shape[0]-1):
#			delta = (array[i+1]-array[i])/float(self.params['discrete'])
#			for j in range(0, self.params['discrete']):
#				newarray[i*self.params['discrete']+j] = array[i] + j*delta
#		newarray[-1] = array[-1]
#		return newarray
#
#	def funcdl(self, z, omegam, omegal):
#		func = lambda z, omegam, omegal: \
#			((1.+z)**2*(1.+omegam*z)-z*(2.+z)*omegal)**(-0.5)
#		dl, dlerr = integrate.quad(func, 0, z, args=(omegam, omegal))
#	
#		if omegam+omegal>1. :
#			k = (omegam+omegal-1.)**(0.5)
#			return sin(k*dl)/k
#		elif omegam+omegal<1.:
#			k = (1.-omegam-omegal)**(0.5)
#			return sinh(k*dl)/k
#		elif omegam+omegal==1.:
#			return dl



if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
