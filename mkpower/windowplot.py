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
	'boxshape' :((),),
	'boxunit' : 4,
}
prefix = 'wpt_'

class WindowFunctionPlot(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

	def execute(self, nprocesses=1):
		params = self.params

		boxshape = params['boxshape']
		boxunit = params['boxunit']
		input_root = params['input_root']

		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		FKPweight = params['FKPweight']

		# Make parent directory and write parameter file.
		#windowN = len(boxshape)
		#boxinf = str(boxshape[0][0])+'x'+\
		#	str(boxshape[0][1])+'x'+str(boxshape[0][2])+'x'+str(boxunit[0])
		#window_fname = input_root+'window_'+boxinf+'_'+resultf+'.npy'
		#w = sp.load(window_fname)
		#if windowN>1: 
		#	window = np.ndarray(shape=[windowN,len(w)])
		#	window[0] = w
		#	for i in range(1,windowN):
		#		boxinf = str(boxshape[i][0])+'x'+\
		#			str(boxshape[i][1])+'x'+str(boxshape[i][2])+'x'+str(boxunit[i])
		#		window_fname = input_root+'window_'+boxinf+'_'+resultf+'.npy'
		#		w = sp.load(window_fname)
		#		window[i] = w
		#else: window = w

		window_fname = input_root+'window_fit_'+resultf+'.npy'
		window = sp.load(window_fname)

		PKcamb = sp.load(params['input_root']+'PKcamb.npy')

		I = 0
		boxinf = str(boxshape[I][0])+'x'+\
			str(boxshape[I][1])+'x'+str(boxshape[I][2])+'x'+str(boxunit[I])
		window_fname = input_root+'WindowF_'+boxinf+'_'+resultf+'.npy'
		window_data0 = sp.load(window_fname)
		window_fname = input_root+'k_WindowF_'+boxinf+'_'+resultf+'.npy'
		window_data0_k = sp.load(window_fname)
		label0 = boxinf

		I = 1
		boxinf = str(boxshape[I][0])+'x'+\
			str(boxshape[I][1])+'x'+str(boxshape[I][2])+'x'+str(boxunit[I])
		window_fname = input_root+'WindowF_'+boxinf+'_'+resultf+'.npy'
		window_data1 = sp.load(window_fname)
		window_fname = input_root+'k_WindowF_'+boxinf+'_'+resultf+'.npy'
		window_data1_k = sp.load(window_fname)
		label1 = boxinf


		plt.figure(figsize=(8,4))
		plt.subplot('111')

		k = PKcamb[0]

		
		#boxinf = str(boxshape[0][0])+'x'+\
		#	str(boxshape[0][1])+'x'+str(boxshape[0][2])+'x'+str(boxunit[0])
		#if windowN>1:
		#	for i in range(windowN):
		#		boxinf = str(boxshape[i][0])+'x'+\
		#			str(boxshape[i][1])+'x'+str(boxshape[i][2])+'x'+str(boxunit[i])
		#		plt.plot(PKcamb[0], window[i], linewidth=2, label=boxinf)
		#else: plt.plot(PKcamb[0], window, linewidth=2, label=boxinf)
		
		plt.plot(PKcamb[0], window, 'r-', linewidth=2.5, label='Fit Line')

		plt.plot(window_data0_k, window_data0, 'go', label=label0)
		plt.plot(window_data1_k, window_data1, 'w^', label=label1)
		plt.loglog()
		#plt.ylim(ymin=1.e-7)	
		plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		plt.title('Window Function')
		#plt.xlabel('$k$')
		plt.ylabel('$|W(k)|^2$')

		plt.legend(loc=0)
		plt.savefig(params['output_root']+'windowfunction_'+resultf+'.eps', 
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
