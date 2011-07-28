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


params_init = {
	'processes' : 1,
	'plot' : False,
	'input_root' : '../newmaps/',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : './map/',

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
			if imap.axes != ('freq', 'ra', 'dec') :
				raise ce.DataError('AXES ERROR!')

			imap_fname = in_root + hr_str + mid[1] + end + '.npy'
			nmap = algebra.load(imap_fname)
			nmap = algebra.make_vect(nmap)
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
			jkmap = algebra.info_array(sp.zeros(imap.shape))
			jkmap.info = dict(imap.info)
			jkmap = algebra.make_vect(jkmap)

			njkmap = algebra.info_array(sp.zeros(nmap.shape))
			njkmap.info = dict(nmap.info)
			njkmap = algebra.make_vect(njkmap)

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

						jkmap = imap.copy()
						njkmap = nmap.copy()

						for i in range(begin0, stop0):
							for j in range(begin1, stop1):
								for k in range(begin2, stop2):
									jkmap[i][j][k] = 0.
									njkmap[i][j][k] = 0.

						jkmap_fname = \
							out_root+hr_str+'jk'+str(num)+mid[0]+end+'.npy'
						algebra.save(jkmap_fname, jkmap)
						jkmap_fname = \
							out_root+hr_str+'jk'+str(num)+mid[1]+end+'.npy'
						algebra.save(jkmap_fname, njkmap)

						num = num + 1

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

if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		JackKnifeErrorMap(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		JackKnifeErrorMap().execute()
