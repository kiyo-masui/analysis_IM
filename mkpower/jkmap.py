#! /usr/bin/env python

import scipy as sp
import numpy as np
from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import matplotlib.pyplot as plt
import multiprocessing as mp
import MakePower

pi = 3.1415926
deg2rad = pi/180.

params_init = {
	'processes' : 1,
	'plot' : False,
	'input_root' : '../newmaps/',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : './map/',

	'boxshape' : (60,12,6),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),

	'jkn0' : 10,
	'jkn1' : 5,
	'jkn2' : 3,

}
prefix = 'jm_'


class JackKnifeErrorMap(object):
	"""Make maps for JackKnife error calculation"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		#Read in the parameter.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix)
		self.feedback = feedback

	def execute(self, nprocesses=1):
		params = self.params

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='jk_')
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']
		all_out_fname_list = []
		all_in_fname_list = []
		n_processes = params['processes']

		#### Process ####
		pol_str = params['polarizations'][0]
		#hr_str = params['hr'][0]
		for hr_str, ii in zip(params['hr'], range(len(params['hr']))):
			end = pol_str
			if (len(params['last'])!=0):
				end = pol_str + params['last'][ii]

			print 'Making JK Map for:' + hr_str[:-1]
			#imap_fname = in_root + hr_str + 'dirty_map_' + pol_str + '.npy'
			imap_fname = in_root + hr_str + mid[0] + end + '.npy'
			imap = algebra.load(imap_fname)
			imap = algebra.make_vect(imap)
			imap = imap - imap.flatten().mean()
			if imap.axes != ('freq', 'ra', 'dec') :
				raise ce.DataError('AXES ERROR!')

			imap_fname = in_root + hr_str + mid[1] + end + '.npy'
			try:
				nmap = algebra.load(imap_fname)
				nmap = algebra.make_vect(nmap)
				if nmap.axes != ('freq', 'ra', 'dec') :
					raise ce.DataError('AXES ERROR!')
			except IOError:
				print 'NO Noise File :: Set Noise to One'
				nmap = algebra.info_array(sp.ones(imap.shape))
				nmap.axes = imap.axes
				nmap = algebra.make_vect(nmap)
			nmap.info = imap.info
			if nmap.axes != ('freq', 'ra', 'dec') :
				raise ce.DataError('AXES ERROR!')

			shape = np.array(imap.shape)
			print shape
			begin0=0
			begin1=0
			begin2=0
			stop0 =0
			stop1 =0
			stop2 =0
			r  = self.fq2r(imap.get_axis('freq'))
			dr = (r.max()-r.min())/params['jkn0']
			ranger = np.searchsorted(r,np.arange(r[0], r[-1], dr))
			num = 0

			#print range(0, shape[1], shape[1]/params['jkn1'])
			#print range(0, shape[2], shape[2]/params['jkn2'])
			#return 0
			#jkmap = algebra.info_array(sp.zeros(imap.shape))
			#jkmap.info = dict(imap.info)
			#jkmap = algebra.make_vect(jkmap)

			#njkmap = algebra.info_array(sp.zeros(nmap.shape))
			#njkmap.info = dict(nmap.info)
			#njkmap = algebra.make_vect(njkmap)

			list_abandon = []

			for ii in range(params['jkn0']):
				begin0 = ranger[ii]
				if begin0==ranger[-1]:
					stop0 = shape[0]
				else:
					stop0 = ranger[ii+1]

				for begin1 in range(0, shape[1], shape[1]/params['jkn1']):
					stop1 = begin1 + shape[1]/params['jkn1']

					for begin2 in range(0, shape[2], shape[2]/params['jkn2']):
						stop2 = begin2 + shape[2]/params['jkn2']

						list_abandon.append(
							[begin0, stop0, begin1, stop1, begin2, stop2])
					num = num + 1

			n_new = n_processes
			n_map = len(list_abandon)
			if n_new <=1:
				for ii in range(n_map):
					jkmap_fname = \
						out_root+hr_str+'jk'+str(ii)+mid[0]+end+'.npy'
					jknmap_fname = \
						out_root+hr_str+'jk'+str(num)+mid[1]+end+'.npy'
					
					self.process_map(
						list_abandon[ii], imap, nmap, jkmap_fname, jknmap_fname)
			elif n_new >32:
				raise ValueError("Processes limit is 32")
			else:
				process_list = range(n_new)
				for ii in xrange(n_new+n_map):
					if ii >= n_new:
						process_list[ii%n_new].join()
						if process_list[ii%n_new].exitcode != 0:
							raise RuntimeError("A thred faild with exit code"
								+ str(process_list[ii%n_new].exitcode))
					if ii < n_map:
						jkmap_fname = \
							hr_str+'jk'+str(ii)+mid[0]+end+'.npy'
						jknmap_fname = \
							hr_str+'jk'+str(ii)+mid[1]+end+'.npy'
						process_list[ii%n_new] = mp.Process(
							target=self.process_map, args=(list_abandon[ii],
								imap, nmap, jkmap_fname, jknmap_fname))
						process_list[ii%n_new].start()




	def process_map(self, abandon, imap, nmap, jkmap_fname, jknmap_fname):
		params = self.params
		out_root = params['output_root']
		jkmap = imap.copy()
		njkmap = nmap.copy()

		begin0 = abandon[0]
		stop0  = abandon[1]
		begin1 = abandon[2]
		stop1  = abandon[3]
		begin2 = abandon[4]
		stop2  = abandon[5]

		for i in range(begin0, stop0):
			for j in range(begin1, stop1):
				for k in range(begin2, stop2):
					jkmap[i][j][k] = 0.
					njkmap[i][j][k] = 0.

		box, nbox = self.fill(jkmap, njkmap)

		algebra.save(out_root + 'fftbox_' + jkmap_fname, box)
		algebra.save(out_root + 'fftbox_' + jknmap_fname, nbox)
		
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

		print "JeckknifeMapPrepare: Filling the FFT BOX"
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

		#return imap, nmap
		return box, nbox

	def fq2r(self, freq, freq0=1.4e9 , c_H0 = 2.99e3, Omegam=0.27, Omegal=0.73):
		"""change the freq to distence"""
		zz =  freq0/freq - 1.
		for i in range(0, zz.shape[0]):
			zz[i] = c_H0*self.funcdl(zz[i], Omegam, Omegal)
		return zz

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

	def discrete(self, array):
		"""discrete the data pixel into small size"""
		newarray = sp.zeros(self.params['discrete']*(array.shape[0]-1)+1)
		for i in range(0, array.shape[0]-1):
			delta = (array[i+1]-array[i])/float(self.params['discrete'])
			for j in range(0, self.params['discrete']):
				newarray[i*self.params['discrete']+j] = array[i] + j*delta
		newarray[-1] = array[-1]
		return newarray

if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		JackKnifeErrorMap(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		JackKnifeErrorMap().execute()
