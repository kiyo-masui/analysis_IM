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
import windowf
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.integrate import romberg
from scipy.integrate import quad

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
	'plot' : True,
	'input_root' : '../newmaps/',
	'output_root' : './',
	'camb_input_root' : './',

	'hr' : (),
	'mid' : (),
	'polarizations' : (),
	'last' : (),

	'FKPweight' : False,
	'FKPpk' : 1.e-3,
	'boxshapelist' : ((),),
	'boxshape' : (60,12,6),
	'boxunitlist' : (),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),

	'OmegaHI' : 1.e-3,
	'Omegam' : 0.24,
	'OmegaL' : 0.76,
	'z' : 1.
}
prefix = 'nl_'
def windowfunction(k, A):
	return 1.e13*A[0]/(1. + (k/(1.e-3*A[1]))**2 + (k/(1.e-3*A[2]))**4)

class TheoryPowerSpectrumMaker(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params
		boxshape = params['boxshape']
		boxunit = params['boxunit']
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='nl_' )
		in_root = params['input_root']
		out_root = params['output_root']
		cambin_root = params['camb_input_root']
		all_out_fname_list = []
		all_in_fname_list = []
		
		#### Process ####
		kiyopy.utils.mkparents(params['output_root'])

		PKcamb_fname = cambin_root + 'PKcamb.npy'
		PKcamb = algebra.load(PKcamb_fname)

		N = len(params['boxunitlist'])
		yy = np.ndarray(shape=(N, 10))
		xx = np.ndarray(shape=(N, 10))
		for params['boxshape'], params['boxunit'], i\
			in zip(params['boxshapelist'], params['boxunitlist'], range(N)):
			params['plot'] = False
			parse_ini.write_params(params, 
				params['output_root']+'params.ini',prefix='wd_' )
			WindowF, k = \
				windowf.WindowFunctionMaker(params['output_root']+'params.ini',
				feedback=self.feedback).execute()
			if yy.shape[1] != WindowF.shape[0]:
				yy.resize((N,WindowF.shape[0]))
				xx.resize((N,WindowF.shape[0]))
			yy[i] =  WindowF.copy()
			xx[i] =  k.copy()

		def chisq(A, y, x, e):
			err = (y - windowfunction(x, A))**2/e**2
			return err

		non0 = yy[0].nonzero()
		y = yy[0].take(non0)[0][:-10]
		x = xx[0].take(non0)[0][:-10]
		non0 = yy[-1].nonzero()
		y = np.append(y, yy[-1].take(non0)[0][10:-4])
		x = np.append(x, xx[-1].take(non0)[0][10:-4])
		err = y.copy()*10.
		err[5:] = err[5:]*1.e-8

		print x.min(), x.max()
		ki = np.logspace(log10(0.01), log10(1.5), num=300)

		A1 = 1.
		A2 = 1.
		A3 = 1.8
		A0 = np.array([A1, A2, A3])
		A, status = leastsq(chisq, A0, args=(y, x, err), maxfev=20000)

		window = windowfunction(PKcamb[0], A)
		#boxinf = str(boxshape[0])+'x'\
		#	+str(boxshape[1])+'x'+str(boxshape[2])+'x'+str(boxunit)
		sp.save(out_root+'window_fit_'+resultf, window)

		CC = 1.
	#	CC = romberg(lambda k2: K(ki,k2)*k2*k2, PKcamb[0].min(), PKcamb[0].max())
#	#	CC = romberg(lambda k2: K(ki,k2)*k2*k2, 1.e-10, 1.e10)
		
		print A
		aaa = A[1]*1.e-3
		bbb = A[2]*1.e-3
		if bbb**4<4*aaa**4:
			CC = 1./(pi*bbb*(2.-(bbb/aaa)**2)**(0.5))
			def g(x):
				return atan((bbb**4 + 2.*aaa**2*x**2)/
					(bbb**2*(4.*aaa**4-bbb**4)**0.5))
			def K(k1, k2):
				return CC/(k1*k2)*(g(k1+k2)-g(k1-k2))
		else:
			mu = bbb**2*(bbb**4-4.*aaa**4)**0.5
			CC = aaa/(pi*2**0.5*((bbb**4+mu)**0.5-(bbb**4-mu)**0.5))
			def g(x):
				return (mu+bbb**4+2*aaa**2*x**2)/(mu-bbb**4-2*aaa**2*x**2)
			def K(k1, k2):
				return CC/(k1*k2)*log(g(k1-k2)/g(k1+k2))

		#def K(k1,k2):
		#	uplim = k1+k2
		#	downlim = np.fabs(k1-k2)
		#	C = 8*pi**2/(k1*k2)*CC
		#	return C*romberg(lambda Q: windowfunction(Q,A)*Q, downlim, uplim)


	#	print CC

		P = interp1d(PKcamb[0], PKcamb[1], kind='cubic')

		#print PKcamb[0].min(), PKcamb[0].max()

		Pnl = np.zeros(len(ki))
		Pnl_err = np.zeros(len(ki))
		for i in range(len(ki)):
			#Pnl[i] = romberg(lambda k1: k1**2*P(k1)*K(k1,ki[i]),
			Pnl[i], Pnl_err = quad(lambda k1: k1**2*P(k1)*K(k1,ki[i]),
				PKcamb[0].min(), PKcamb[0].max(), limit=200)
		#Pnl = sp.load(out_root+'nonlPK_'+resultf+'.npy')	

		CCC = romberg(lambda k1: k1**2*K(k1, 0.01), ki.min(), ki.max())
		print CCC
		#Pnl = Pnl/CCC

		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		z = params['z']
		a3 = (1+z)**(-3)
		Tb = 0.3e-3 * (OmegaHI/1.e-3) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
			* ((1.+z)/2.5)**0.5

		#Pnl = Pnl*(Tb**2)
		#PKcamb[1] = PKcamb[1]*(Tb**2)

		#print Pnl

		sp.save(out_root+'nonlPK_'+resultf, Pnl)
		sp.save(out_root+'k_nonlPK_'+resultf, ki)

		if self.plot==True:
			#plt.figure(figsize=(6,6))
			##print k
			#plt.subplot('111')

			##kj = sp.linspace(0,PKcamb[0][-1], num=500)
			##KI = np.zeros(500)
			##for j in sp.linspace(ki.min(),ki.max(), num=20):
			##	for i in range(500):
			##		KI[i] = K(j, kj[i])
			##	#plt.plot(kj, KI, label=str(ki))
			##	plt.plot(kj, KI, 'r-', linewidth=1)

			#plt.semilogy()
			##plt.loglog()
			##plt.ylim(ymin=1.e-0)	
			#plt.xlim(xmin=0, xmax=ki.max())
			#plt.title('Coupling Kernels')
			#plt.xlabel('$k$')
			#plt.ylabel('$K(k, k_i)$')
			##plt.legend()


			#plt.savefig(out_root+'Ki.eps', format='eps')

			plt.figure(figsize=(8,4))
			plt.subplot('111')
			#plt.plot(PKcamb[0], window, 'b--', linewidth=1,
			#	label='Fitted Window Function')
			plt.plot(PKcamb[0], PKcamb[1], 'g-', linewidth=1,
				label='Camb Power Spectrum')
			plt.plot(ki, Pnl, 'r-', linewidth=1, 
				label='Power Spectrum')
			plt.loglog()
			plt.xlim(xmin=ki.min(), xmax=ki.max())

			plt.legend()
			plt.savefig(out_root+'nonlPK.eps', format='eps')

			plt.show()
			#print 'Finished @_@ '
		return PKcamb

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
		TheoryPowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		TheoryPowerSpectrumMaker().execute()
