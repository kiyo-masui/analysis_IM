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
	'FKPpk' : 1.e-3,
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
		fkpp = params['FKPpk']


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
		fkpp = params['FKPpk']

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, 
			params['output_root']+'params.ini',prefix='wd_' )
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
		nmap.info = imap.info
		if FKPweight:
			for i in range(nmap.shape[0]):
				nmap[i] = 1./(1.+nmap[i]*fkpp)
				#nmap[i] = nmap[i]/(1.+nmap[i]*fkpp)


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
				nmap2[i] = 1./(1.+nmap2[i]*fkpp)
				#nmap2[i] = nmap2[i]/(1.+nmap2[i]*fkpp)


		weight, weight2 = self.fill(nmap, nmap2)

		#MakePower.Filling(nmap, weight, r, ra, de, boxinf, mapinf)
		#MakePower.Filling(nmap2, weight2, r, ra, de, boxinf, mapinf)
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
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		inputa.imag = 0.
		inputa.real = weight
		FFTW.execute(fft)

		#inputa = np.zeros(params['boxshape'], dtype=complex)
		outputb = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputb, direction='forward', flags=['measure'])
		inputa.imag = 0.
		inputa.real = weight2
		FFTW.execute(fft)

		V = params['boxunit']**3

		fftbox = (outputa*(outputb.conjugate())).real
		fftbox = fftbox*V/normal

		WindowF = np.zeros(40)
		k = np.zeros(40)
		WindowF2 = np.zeros(shape=(10, 10))
		k2 = np.zeros(shape=(2, 10))
		MakePower.Make(fftbox, WindowF, k, WindowF2, k2)
		kunit = 2.*pi/(params['boxunit'])
		k = k*kunit
		k2 = k2*kunit
		#WindowF = WindowF#/V
		#WindowF2 = WindowF2#/V

		return WindowF, k

	def fill(self, imap, nmap):
		params = self.params
		
		mapshape = np.array(imap.shape)

		r  = self.fq2r(imap.get_axis('freq'))
		ra = imap.get_axis('ra')*deg2rad
		de = imap.get_axis('dec')*deg2rad
		ra0= ra[int(ra.shape[0]/2)]
		ra = ra - ra0
		dra= ra.ptp()/ra.shape[0]
		dde= de.ptp()/de.shape[0]


		#print r.min(), r.max()
		#print self.xyz(ra.min(), de.min(), r.min())
		#print self.xyz(ra.max(), de.min(), r.min())
		#print self.xyz(ra.min(), de.max(), r.min())
		#print self.xyz(ra.max(), de.max(), r.min())
		#print self.xyz(ra.min(), de.min(), r.max())
		#print self.xyz(ra.max(), de.min(), r.max())
		#print self.xyz(ra.min(), de.max(), r.max())
		#print self.xyz(ra.max(), de.max(), r.max())

		###return 0

		mapinf = [ra.min(), dra, de.min(), dde]
		mapinf = np.array(mapinf)

		box = algebra.info_array(sp.zeros(params['boxshape']))
		box.axes = ('x','y','z')
		box = algebra.make_vect(box)
		
		box_xrange = params['Xrange']
		box_yrange = params['Yrange']
		box_zrange = params['Zrange']
		box_unit = params['boxunit']
		box_disc = params['discrete']

		box_x = np.arange(box_xrange[0], box_xrange[1], box_unit/box_disc)
		box_y = np.arange(box_yrange[0], box_yrange[1], box_unit/box_disc)
		box_z = np.arange(box_zrange[0], box_zrange[1], box_unit/box_disc)

		#print box_x.shape
		#print box_y.shape
		#print box_z.shape

		boxshape = np.array(box.shape)*box_disc

		boxinf0 = [0, 0, 0]
		boxinf0 = np.array(boxinf0)
		boxinf1 = [boxshape[0], boxshape[1], boxshape[2]]
		boxinf1 = np.array(boxinf1)

		print "MapPrepare: Filling the FFT BOX"
		MakePower.Filling(
			imap, r, mapinf, box, boxinf0, boxinf1, box_x, box_y, box_z)


		nbox = algebra.info_array(sp.zeros(params['boxshape']))
		nbox.axes = ('x','y','z')
		nbox = algebra.make_vect(nbox)


		MakePower.Filling(
			nmap, r, mapinf, nbox, boxinf0, boxinf1, box_x, box_y, box_z)

		#return imap, nmap
		return box, nbox

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
