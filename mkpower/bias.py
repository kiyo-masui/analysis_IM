#! /usr/bin/env python
'''
	This module is used to calculate the bias(Transfer Function)
	caused by the mode subtraction. 

	The bias is defined as the ratio of the power spectrum 
	between the simulation maps with and without mode subtraction.

'''

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
from scipy.optimize import leastsq

params_init = {
	'processes' : 1,
	'plot' : True,
	'input_root' : '../newmaps/',
	'resultf' : '',
	'output_root' : './',

	'simmap_root' : './',

	'kbinNum' : 20,
	'kmin' : None,
	'kmax' : None,

	'FKPweight' : False,

	'OmegaHI' : 1.e-3,
	'Omegam' : 0.23,
	'OmegaL' : 0.74,
	'z' : 1,
	'PKunit' : 'K',
	'cross' : False,
}

pi = 3.1415926
deg2rad = pi/180.
prefix = 'bc_'

class BiasCalibrate(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict,
			params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):

		params = self.params
		out_root = params['output_root']
		in_root = params['input_root']
		resultf = params['resultf']
		#resultf = params['hr'][0]
		#if len(params['last']) != 0:
		#	resultf = resultf + params['last'][0]
		#resultf = resultf + '-' + params['hr'][1]
		#if len(params['last']) != 0:
		#	resultf = resultf + params['last'][1]

		FKPweight = params['FKPweight']


		# Read the Power Spectrum Result
		k = sp.load(in_root + 'k_combined_' + resultf + '.npy')
		if FKPweight:
			pkvar = sp.load(in_root + 'PKvar_combined_fkp_' + resultf + '.npy')
			pk = sp.load(in_root + 'PK_combined_fkp_' + resultf + '.npy')
			pke = sp.load(in_root + 'PKeach_fkp_' + resultf + '.npy')
		else:
			pkvar = sp.load(in_root + 'PKvar_combined_' + resultf + '.npy')
			pk = sp.load(in_root + 'PK_combined_' + resultf + '.npy')
			pke = sp.load(in_root + 'PKeach_' + resultf + '.npy')

		non0 = pk.nonzero()
		pk = pk.take(non0)[0]
		k = k.take(non0)[0]
		pkvar = pkvar.take(non0)[0]
		pkerr = np.ndarray(shape=(2, len(non0[0])))
		pkerr[0] = pkvar
		pkerr[1] = pkvar

		# Read the Pwer Spectrum without mode subtraction
		simmap_root = params['simmap_root']
		k0 = sp.load(simmap_root + 'k_combined_simmaps_with_beam.npy')
		pk0= sp.load(simmap_root + 'PK_combined_simmaps_with_beam.npy')

		# Test if k and k0 are match
		if (k-k0).any():
			print "k and k0 are not match!!"
			return 0


		# Read the Theoretical Power Spectrum
		#PKcamb = sp.load(in_root + 'PKcamb.npy')

		#OmegaHI = params['OmegaHI']
		#Omegam = params['Omegam']
		#OmegaL = params['OmegaL']
		#z = params['z']
		#a3 = (1+z)**(-3)
		#Tb = (OmegaHI) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
		#	* ((1.+z)/2.5)**0.5
		#if params['PKunit']=='mK':
		#	Tb = Tb/1.e-3

		#xx = ((1.0/Omegam)-1.0)/(1.0+z)**3
		#num = 1.0 + 1.175*xx + 0.3046*xx**2 + 0.005335*xx**3
		#den = 1.0 + 1.875*xx + 1.021 *xx**2 + 0.1530  *xx**3

		#G = (1.0 + xx)**0.5/(1.0+z)*num/den
		##print G**2*(1+z)**(-2)

		#PKcamb[1] = PKcamb[1]*(G**2*(1+z)**(-2))
		#if params['cross']: 
		#	PKcamb[1] = PKcamb[1]*Tb
		#else: 
		#	PKcamb[1] = PKcamb[1]*(Tb**2)

		#P = interp1d(PKcamb[0], PKcamb[1], kind='cubic')
		#pk_th = P(k)

		dpke = pke.copy()
		dpke[:] = pk0/dpke[:]
		print dpke

		dpk = pk0/pk
		print '\tFor A maps:'
		pkA, pkerrA, dpkA = self.getbias(pke, 0, 3, pk0)
		print '\tFor B maps:'
		pkB, pkerrB, dpkB = self.getbias(pke, 3, 6, pk0)
		print '\tFor C maps:'
		pkC, pkerrC, dpkC = self.getbias(pke, 6, 9, pk0)
		print '\tFor D maps:'
		pkD, pkerrD, dpkD = self.getbias(pke, 9, 12, pk0)

		B = interp1d(k, dpk, kind='cubic')

		ki = np.logspace(log10(k.min()+0.001), log10(k.max()-0.001), num=300)
		bias = B(ki)

		sp.save(out_root + 'Bias_k_' + resultf, k)
		sp.save(out_root + 'Bias_B_' + resultf, dpk)
		sp.save(out_root + 'Bias_Be_' + resultf, dpke)


		if self.plot==True:
			plt.figure(figsize=(8,4))
			plt.subplot('111')
			plt.errorbar(k, dpk, pkerr, fmt='o', c='k', capsize=4.5)
			for i in range(12):
				plt.scatter(k, dpke[i], c='r' )
			#plt.errorbar(k, dpkA, pkerrA, fmt='o', c='r', capsize=4.5)
			#plt.errorbar(k, dpkB, pkerrB, fmt='o', c='g', capsize=4.5)
			#plt.errorbar(k, dpkC, pkerrC, fmt='o', c='b', capsize=4.5)
			#plt.errorbar(k, dpkD, pkerrD, fmt='o', c='y', capsize=4.5)
			plt.plot(ki, bias)
			#plt.scatter(k, dpk)
			plt.loglog()
			#plt.semilogx()
			plt.ylim(ymin=dpk.min())
			plt.xlim(xmin=k.min(), xmax=k.max())
			plt.title('Power Spectrum Bias')
			plt.xlabel('$k$')
			plt.ylabel('$dP(k) (mk^{2}(h^{-1}Mpc)^3)$')
			plt.show()

	def getbias(self, pke, begin, end, pk0):
		print pke[begin:end]
		pkA = pke[begin:end].mean(axis=0)
		pke[begin:end] = (pke[begin:end][:]-pkA)**2
		pkvarA = np.sum(pke[begin:end], axis=0)
		pkvarA = pkvarA/(end-begin)
		pkvarA = np.sqrt(pkvarA)
		print pk0
		print pkA
		print pkvarA
		pkerrA = np.ndarray(shape=(2, len(pkA)))
		pkerrA[0] = pkvarA
		pkerrA[1] = pkvarA
		dpkA = pk0/pkA
		print dpkA
		return pkA, pkerrA, dpkA

if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		BiasCalibrate(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		BiasCalibrate().execute()
