#! /usr/bin/env python

import ctypes
import scipy as sp
import numpy as np
from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp
import random

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import matplotlib.pyplot as plt
import MakePower


pi = 3.1415926
deg2rad = pi/180.

params_init = {
	'processes' : 1,
	'plot' : True,
	'input_root' : '../newmaps/',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_','noise_inv_diag_'),
	'polarizations' : ('I',),
	'output_root' : './',
	'last' : (),

	'boxshape' : (60,12,6),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),
}
prefix = 'pre_'

class Prepare(object):
	"""Remove the Big Peak in the map"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(
			parameter_file_or_dict, params_init, prefix=prefix)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='pk_')
		hr = params['hr']
		mid = params['mid']
		last = params['last']
		all_out_fname_list = []
		all_in_fname_list = []
		pol_str = params['polarizations'][0]
		n_processes = params['processes']
		
		#### Process ####
		n_new = n_processes -1
		n_map = len(hr)

		if n_new <=0:
			for hr_str, ii in zip(params['hr'],range(len(params['hr']))):
				end = pol_str
				if len(last)!=0:
					end = end + last[ii]
				#imap_fname = in_root + hr_str + 'dirty_map_' + pol_str + '.npy'
				#imap_fname = in_root + hr_str + mid + pol_str + '.npy'
				imap_fname = hr_str + mid[0] + end + '.npy'
				nmap_fname = hr_str + mid[1] + end + '.npy'

				self.process_map(imap_fname, nmap_fname)
		elif n_new >32:
			raise ValueError("Processes limit is 32")
		else: 
			process_list = range(n_new)
			for ii in xrange(n_map + n_new):
				if ii >= n_new:
					process_list[ii%n_new].join()
					if process_list[ii%n_new].exitcode != 0:
						raise RuntimeError("A thred faild with exit code"
							+ str(process_list[ii%n_new].exitcode))
				if ii < n_map:
					end = pol_str
					if len(last)!=0:
						end = end + last[ii]
					imap_fname = hr[ii] + mid[0] + end + '.npy'
					nmap_fname = hr[ii] + mid[1] + end + '.npy'
					#mock_fname = hr[ii] + 'mock_map_' + end + '.npy'
					process_list[ii%n_new] = mp.Process(
						target=self.process_map, 
						args=(imap_fname, nmap_fname))
						#args=(imap_fname, nmap_fname, mock_fname))
					process_list[ii%n_new].start()
		return 0

	def process_map(self, imap_fname, nmap_fname, mock_fname=None):
		params = self.params
		out_root = params['output_root']
		in_root = params['input_root']
		
		imap = algebra.load(in_root + imap_fname)
		imap = algebra.make_vect(imap)
		print imap.flatten().mean()
		imap = imap - imap.flatten().mean()
		if imap.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		try:
			nmap = algebra.load(in_root + nmap_fname)
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
		if nmap.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		if mock_fname != None:
			mmap = algebra.info_array(
				2.*np.random.rand(imap.shape[0],imap.shape[1], imap.shape[2])-0.5)
			mmap.axes = imap.axes
			mmap = algebra.make_vect(mmap)
			box, nbox, mbox = self.fill(imap, nmap, mmap)
			pkrm_nfname = out_root + 'fftbox_' +  mock_fname
			algebra.save(pkrm_nfname, mbox)
		else:
			box, nbox = self.fill(imap, nmap)

		pkrm_fname = out_root + 'fftbox_' + imap_fname
		algebra.save(pkrm_fname, box)

		pkrm_nfname = out_root + 'fftbox_' +  nmap_fname
		algebra.save(pkrm_nfname, nbox)

#		print mmap.shape

#		print 'Removing Peak... Map:' + hr_str[:-1]
#		self.pkrm(imap,nmap,
#			out_root+'pkrm'+hr_str+'dirty_map_'+pol_str+'.png', threshold=2.5)

		#print imap.flatten().max()
		#print imap.flatten().min()
		#print box.flatten().max()
		#print box.flatten().min()


	def fill(self, imap, nmap, mmap=None):
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

		#nbox = algebra.info_array(sp.ones(params['boxshape']))
		#nbox.axes = ('x','y','z')
		#nbox = algebra.make_vect(nbox)

		MakePower.Filling(
			nmap, r, mapinf, nbox, boxinf0, boxinf1, box_x, box_y, box_z)


		if mmap != None:
			mbox = algebra.info_array(sp.zeros(params['boxshape']))
			mbox.axes = ('x','y','z')
			mbox = algebra.make_vect(mbox)
			MakePower.Filling(
				mmap, r, mapinf, mbox, boxinf0, boxinf1, box_x, box_y, box_z)
			return box, nbox, mbox
		else:
			return box, nbox

		#return imap, nmap

	def pkrm(self, imap, nmap, fname, threshold=2.0):
		freq = imap.get_axis('freq')/1.e6
		if self.plot==True:
			plt.figure(figsize=(8,8))
			plt.subplot(211)
			plt.title('Map with peak remove')
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) Without Foreground')
			#for i in range(13,14):
			#	for j in range(25, 26):
			for i in range(0,imap.shape[2]):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap.swapaxes(0,2)[i][j])

		dsigma = 50
		while(dsigma>0.01):
			sigma = imap.std(0)
			for i in range(0, imap.shape[0]):
				good = [np.fabs(imap[i]).__lt__(threshold*sigma)]
				choicelist = [imap[i]]
				imap[i] = np.select(good, choicelist)
				choicelist = [nmap[i]]
				nmap[i] = np.select(good, choicelist)
			dsigma = (sigma.__sub__(imap.std(0))).max()
			#print dsigma
		#print '\n'

		if self.plot==True:
			plt.subplot(212)
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) Without Foreground')
			#for i in range(13,14):
			#	for j in range(25, 26):
			for i in range(0,imap.shape[2]):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap.swapaxes(0,2)[i][j])

			plt.savefig(fname, format='png')
			plt.show()
			#plt.ylim(-0.0001,0.0001)

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
		Prepare(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		Prepare().execute()
