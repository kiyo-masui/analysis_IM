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
	'processes' : 1,
	'input_root' : '../newmaps/',
	'output_root' : './',
	'hr' : (),
	'last' : (),
	'jknumber' : 100.,
	'FKPweight' : False,
	'OmegaHI' : 1.e-3,
	'Omegam' : 0.24,
	'OmegaL' : 0.76,
	'z' : 1.
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

		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		FKPweight = params['FKPweight']

		# Make parent directory and write parameter file.
		k = sp.load(params['input_root'] + 'k_'+resultf+'.npy')
		PK = sp.load(params['input_root'] + 'PK_'+resultf+'.npy')
		PKjk = sp.load(params['input_root'] + 'PKjk_'+resultf+'.npy')
		PKvar = sp.load(params['input_root'] + 'PKvar_'+resultf+'.npy')
		non0 = PK.nonzero()

		PK = PK.take(non0)[0]
		k = k.take(non0)[0]
		PKmean = PKjk.mean(axis=0)
		#PKjk[:] = (PKjk[:]-PKmean)**2
		#PKvar = np.sum(PKjk[:-1], axis=0).take(non0)[0]
		#PKvar = PKvar*(params['jknumber']-1)/params['jknumber']
		#PKvar = np.sqrt(PKvar)
		PKerr = np.ndarray(shape=(2,len(non0[0])))
		PKerr[0] = PKvar.take(non0)[0]
		PKerr[1] = PKvar.take(non0)[0]

		for i in range(len(PKerr[0])):
			if PKerr[0][i]>=PK[i]:
				PKerr[0][i] = PK[i]-1.e-10
		#PKerr[0][PKerr[0]>=PK] = PK-1.e-15
		#print PKerr
		#PKerr = np.ndarray(shape=(2,len(non0[0])))
		#PKerr[0] = PK - PKvar.take(non0)[0]
		#PKerr[1] = PK + PKvar.take(non0)[0]

		#print PK
		#print PK.take(non00)
		#print k[non0]
		#return 0

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

		PKcamb = sp.load(params['input_root']+'PKcamb.npy')
		PKnonl = sp.load(params['input_root']+'nonlPK_'+resultf+'.npy')
		PKnonl_k = sp.load(params['input_root']+'k_nonlPK_'+resultf+'.npy')
		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		z = params['z']
		a3 = (1+z)**(-3)
		Tb = 0.3e-3 * (OmegaHI/1.e-3) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
			* ((1.+z)/2.5)**0.5

		PKcamb[1] = PKcamb[1]*(Tb**2)

		plt.figure(figsize=(8,8))
		#print k
		#print PK
		plt.subplot('211')
		
		#for i in range(0,50):
		#	plt.scatter(k, PKjk[i].take(non0), alpha=0.5, c='b')
		#for i in range(200,250):
		#	plt.scatter(k, PKjk[i].take(non0), alpha=0.5, c='r')
		#plt.loglog()
		#plt.ylim(ymin=1.e-17)	
		#plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		#plt.show()
		#return 

		#plt.scatter(k, PK, c='b')
		plt.plot(PKcamb[0], PKcamb[1], 'g-', linewidth=2, label='Theoretical')
		plt.plot(PKnonl_k, PKnonl, 'k-', linewidth=2,
			label='Theoretical Convolved wiht Window')
		plt.errorbar(k, PK, PKerr, fmt='o', c='b', label='Noise Inv Weight')
		if FKPweight:
			plt.errorbar(kfkp, PKfkp, PKerrfkp, fmt='o', c='r',
				label='FPK Weight', capsize=4)
		plt.loglog()
		plt.ylim(ymin=1.e-6)	
		plt.xlim(xmin=0.01, xmax=0.9)
		#plt.xlim(xmin=0.01, xmax=k.max()+0.1*k.max())
		#plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		plt.title('Power Spectrum')
		#plt.xlabel('$k$')
		plt.ylabel('$P(k) (Kelvin^{2}(h^{-1}Mpc)^3)$')
		#plt.legend()

		PK = PK*k*k*k/2./pi/pi
		PKerr = PKerr*k*k*k/2./pi/pi
		if FKPweight:
			PKfkp = PKfkp*kfkp*kfkp*kfkp/2./pi/pi
			PKerrfkp = PKerrfkp*kfkp*kfkp*kfkp/2./pi/pi
		PKcamb[1] = PKcamb[1]*PKcamb[0]*PKcamb[0]*PKcamb[0]/2./pi/pi
		PKnonl = PKnonl*PKnonl_k**3/2./pi/pi
		#print PK
		plt.legend(loc=8)
		plt.subplot('212')
		#plt.scatter(k, PK, c='b')
		plt.plot(PKcamb[0], PKcamb[1], 'g-', linewidth=2, label='Theoretical')
		plt.plot(PKnonl_k, PKnonl, c='k', linewidth=2,
			label='Theoretical Convolved wiht Window')
		plt.errorbar(k, PK, PKerr, fmt='o', c='b', 
			label='Noise Inv Weight', capsize=4, elinewidth=2)
		if FKPweight:
			plt.errorbar(kfkp, PKfkp, PKerrfkp, fmt='o', c='r',
				label='FPK Weight',capsize=4,elinewidth=2)
		plt.loglog()
		plt.ylim(ymin=1.e-12)	
		plt.xlim(xmin=0.01, xmax=0.9)
		#plt.xlim(xmin=0.01, xmax=k.max()+0.1*k.max())
		#plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		#plt.xlim(xmin=k.min(), xmax=k.max())
		plt.xlabel('$k (h Mpc^{-1})$')
		plt.ylabel('$\Delta^2 (Kelvin^{2})$')
		#plt.legend()
		#plt.show()
		plt.legend(loc=0)
		plt.savefig(params['output_root']+'power_err_'+resultf+'.eps', 
			format='eps')

		plt.show()

	def xyz(self, ra, de, r, ra0=0.):
		x = r*sin(0.5*pi-de)*cos(ra-ra0)
		y = r*sin(0.5*pi-de)*sin(ra-ra0)
		z = r*cos(0.5*pi-de)
		return x, y, z

	def fq2r(self, freq, freq0=1.4e9 , c_H0 = 2.99e3, Omegam=0.27, Omegal=0.73):
		"""change the freq to distence"""
		zz =  freq0/freq - 1.
		for i in range(0, zz.shape[0]):
			zz[i] = c_H0*self.funcdl(zz[i], Omegam, Omegal)
		return zz
	
	def discrete(self, array):
		"""discrete the data pixel into small size"""
		newarray = sp.zeros(self.params['discrete']*(array.shape[0]-1)+1)
		for i in range(0, array.shape[0]-1):
			delta = (array[i+1]-array[i])/float(self.params['discrete'])
			for j in range(0, self.params['discrete']):
				newarray[i*self.params['discrete']+j] = array[i] + j*delta
		newarray[-1] = array[-1]
		return newarray

	def funcdl(self, z, omegam, omegal):
		func = lambda z, omegam, omegal: \
			((1.+z)**2*(1.+omegam*z)-z*(2.+z)*omegal)**(-0.5)
		dl, dlerr = integrate.quad(func, 0, z, args=(omegam, omegal))
	
		if omegam+omegal>1. :
			k = (omegam+omegal-1.)**(0.5)
			return sin(k*dl)/k
		elif omegam+omegal<1.:
			k = (1.-omegam-omegal)**(0.5)
			return sinh(k*dl)/k
		elif omegam+omegal==1.:
			return dl



if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
