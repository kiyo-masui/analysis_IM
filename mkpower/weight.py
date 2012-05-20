#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
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
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : 'cleaned_clean_map_',
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : './',
}
prefix = 'wt_'

class NoiseInversWeight(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='pk_')
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']
		all_out_fname_list = []
		all_in_fname_list = []
		
		#### Process ####
		pol_str = params['polarizations'][0]
		#hr_str = params['hr'][0]
		for hr_str, ii in zip(params['hr'],range(len(params['hr']))):
			end = pol_str
			if len(last)!=0:
				end = end + last[ii]
			imap_fname = in_root + hr_str + mid[0] + end + '.npy'
			imap = algebra.load(imap_fname)
			imap = algebra.make_vect(imap)
			if imap.axes != ('freq', 'ra', 'dec') :
				raise ce.DataError('AXES ERROR!')

			nmap_fname = in_root + hr_str + mid[1] + end + '.npy'
			nmap = algebra.load(nmap_fname)
			nmap = algebra.make_vect(nmap)

			#invers noise weight
			print 'Inverse Noise Weight... Map:' + hr_str[:-1]
			self.weight(imap, nmap, 
				out_root+hr_str+'wt_cleaned_clean_map_'+end+'.png')

			dmap_fname = out_root + 'wt_' + hr_str + mid[0] + end + '.npy'
			algebra.save(dmap_fname, imap)
			all_out_fname_list.append(
				kiyopy.utils.abbreviate_file_path(dmap_fname))

			nmap_fname = out_root + 'wt_' + hr_str + mid[1] + end + '.npy'
			algebra.save(nmap_fname, nmap)
			all_out_fname_list.append(
				kiyopy.utils.abbreviate_file_path(nmap_fname))

		return 0



	def weight(self, imap, nmap, fname):
		freq = imap.get_axis('freq')/1.e6
		if self.plot==True:
			plt.figure(figsize=(8,8))
			plt.subplot(211)
			plt.title('Map with inverse noise weighted')
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) (clean map)')
			for i in range(0,imap.shape[2]):
				#for j in range(1, 2):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap.swapaxes(0,2)[i][j])
		for i in range(0, imap.shape[0]):
			for j in range(0, imap.shape[1]):
				for k in range(0, imap.shape[2]):
					imap[i][j][k] = (imap[i][j][k]*nmap[i][j][k])
			
		if self.plot==True:
			plt.subplot(212)
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) (dirty map)')
			for i in range(0,imap.shape[2]):
				#for j in range(1, 2):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap.swapaxes(0,2)[i][j])
			plt.savefig(fname, format='png')
			#plt.show()

	def xyzv(self, ra, de, r, ra0=0.):
		x = r*sin(0.5*pi-de)*cos(ra-ra0)
		y = r*sin(0.5*pi-de)*sin(ra-ra0)
		z = r*cos(0.5*pi-de)
		v = r**2*sin(0.5*pi-de)
		return x, y, z, v

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
		NoiseInversWeight(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		NoiseInversWeight().execute()
