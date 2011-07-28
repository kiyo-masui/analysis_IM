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
	'plot' : True,
	'input_root' : '../newmaps/',
	'output_root' : './',

	'FKPweight' : False,
	'boxshape' : (60,12,6),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),

	'hr' : (),
	'mid' : (),
	'polarizations' : (),
	'last' : (),
}
prefix = 'wd_'


class WindowFunctionMaker(object):
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
		
		FKPweight = params['FKPweight']
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']


		WindowF_fname = out_root+'WindowF_'+\
			str(boxshape[0])+'x'+str(boxshape[1])+'x'+\
			str(boxshape[2])+'x'+str(boxunit)+'_'+resultf
		kWindowF_fname = out_root+'k_WindowF_'+\
			str(boxshape[0])+'x'+str(boxshape[1])+'x'+\
			str(boxshape[2])+'x'+str(boxunit)+'_'+resultf

		print WindowF_fname

		try:
			WindowF = sp.load(WindowF_fname+'.npy')
			k = sp.load(kWindowF_fname+'.npy')
		except IOError:
			print '\tWindow Functin ReMake'
			WindowF, k = self.GetWindowFunctionData()

		non0 = WindowF.nonzero()
		sp.save(WindowF_fname, WindowF)
		sp.save(kWindowF_fname, k)

		#txtf = open(out_root+'window_for_idl.txt', 'w')
		#try:
		#	for i in range(len(WindowF)):
		#		if WindowF[i]==0: continue
		#		print >>txtf, '{0} {1}'.format(k[i], WindowF[i])
		#finally:
		#	txtf.close()

		return WindowF, k


	def GetWindowFunctionData(self):
		params = self.params
		boxshape = params['boxshape']
		boxunit = params['boxunit']
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]
		
		FKPweight = params['FKPweight']
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, 
			params['output_root']+'params.ini',prefix='pk_' )
		all_out_fname_list = []
		all_in_fname_list = []
		
		#### Process ####
		pol_str = params['polarizations'][0]
		hr_str = params['hr'][0]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][0]
		imap_fname = in_root + hr_str + mid[0] + end + '.npy'
		imap = algebra.load(imap_fname)
		imap = algebra.make_vect(imap)
		if imap.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')
		nmap_fname = in_root + hr_str + mid[1] + end + '.npy'
		try:
			nmap = algebra.load(nmap_fname)
			nmap = algebra.make_vect(nmap)

			bad = nmap<1.e-5*nmap.flatten().max()
			nmap[bad] = 0.
			non0 = nmap.nonzero()
			#imap[non0] = imap[non0]/nmap[non0]
		except IOError:
			print 'NO Noise File :: Set Noise to One'
			nmap = algebra.info_array(sp.ones(imap.shape))
			nmap.axes = imap.axes
			nmap = algebra.make_vect(nmap)
		if FKPweight:
			for i in range(nmap.shape[0]):
				nmap[i] = nmap[i]/(1.+nmap[i]*1.e4)


		hr_str = params['hr'][1]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][1]
		nmap_fname = in_root + hr_str + mid[1] + end + '.npy'
		try:
			nmap2 = algebra.load(nmap_fname)
			nmap2 = algebra.make_vect(nmap2)
	
			bad = nmap2<1.e-5*nmap2.flatten().max()
			nmap2[bad] = 0.
			non0 = nmap2.nonzero()
		except IOError:
			print 'NO Noise File :: Set Noise to One'
			nmap2 = algebra.info_array(sp.ones(imap.shape))
			nmap2.axes = imap.axes
			nmap2 = algebra.make_vect(nmap2)
		if FKPweight:
			for i in range(nmap.shape[0]):
				nmap2[i] = nmap2[i]/(1.+nmap2[i]*1.e4)

		r  = self.discrete(self.fq2r(imap.get_axis('freq')))
		ra = self.discrete(imap.get_axis('ra'))*deg2rad
		de = self.discrete(imap.get_axis('dec'))*deg2rad
		ra0= ra[int(ra.shape[0]/2)]
		ra = ra - ra0
		dr = r.ptp()/r.shape[0]
		dra= ra.ptp()/ra.shape[0]
		dde= de.ptp()/de.shape[0]
		disc_n = params['discrete']
		mapinf = [dr, dra, dde, disc_n]
		mapinf = np.array(mapinf)

		xrange0 = params['Xrange'][0]
		yrange0 = params['Yrange'][0]
		zrange0 = params['Zrange'][0]
		boxunit = params['boxunit']
		shapex = params['boxshape'][2]
		shapera = ra.shape[0]
		V = params['boxunit']**3
		boxinf = [xrange0, yrange0, zrange0, boxunit]
		boxinf = np.array(boxinf)

		weight = algebra.info_array(sp.zeros(params['boxshape']))
		weight.axes = ('x','y','z')
		weight = algebra.make_vect(weight)
		weight2 = algebra.info_array(sp.zeros(params['boxshape']))
		weight2.axes = ('x','y','z')
		weight2 = algebra.make_vect(weight2)

		MakePower.Filling(nmap, nmap2, 
			weight, weight2, r, ra, de, boxinf, mapinf)
		#weight_fname = in_root + 'Weight_' + resultf + '.npy'
		#weight = algebra.load(weight_fname)
		#weight = algebra.make_vect(weight)

		#weight_fname = in_root + 'Weight2_' + resultf + '.npy'
		#weight2 = algebra.load(weight_fname)
		#weight2 = algebra.make_vect(weight2)
		
		normal = (weight**2).flatten().sum()
		normal2 = (weight2**2).flatten().sum()
		normal = sqrt(normal)*sqrt(normal2)


		print "WindowFunctionMaker: FFTing "
		inputa = np.zeros(params['boxshape'], dtype=complex)
		inputa.real = weight
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		FFTW.execute(fft)
		weight = outputa.real**2 + outputa.imag**2
		inputa = np.zeros(params['boxshape'], dtype=complex)
		inputa.real = weight2
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		FFTW.execute(fft)
		weight2 = outputa.real**2 + outputa.imag**2
		fftbox = weight.__pow__(0.5)*weight2.__pow__(0.5)
		fftbox = fftbox*V*V#/normal

		WindowF = np.zeros(40)
		k = np.zeros(40)
		WindowF2 = np.zeros(shape=(10, 10))
		k2 = np.zeros(shape=(2, 10))
		MakePower.Make(fftbox, WindowF, k, WindowF2, k2)
		kunit = 2.*pi/(params['boxunit'])
		k = k*kunit
		k2 = k2*kunit
		WindowF = WindowF#/V
		WindowF2 = WindowF2#/V

		return WindowF, k

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
		WindowFunctionMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		WindowFunctionMaker().execute()
