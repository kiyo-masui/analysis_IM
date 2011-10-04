#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import matplotlib.pyplot as plt
import mkpower
from mpi4py import MPI


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
	'plot' : False,
	'input_root' : '../../../jkmap/map/',
	'jknumber' : 9,
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : '../../../powerresult/',

	'boxshape' : (128,128,128),
	'boxunit' : 15., # in unit Mpc
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-64*15,64*15),
	'Zrange' : (0.,128*15),
	'kbinNum' : 20,

	'FKPweight' : False,
}
prefix = 'jk_'

class JackKnifeError(object):
	"""Calculate The JackKnifeErrot for Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=1):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback
	
	def execute(self, nprocesses=1):
		
		comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
		size = comm.Get_size()

		params = self.params
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		# Make parent directory and write parameter file.
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='pk_')
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']

		FKPweight = params['FKPweight']
		n_processes = params['processes']

		#### Process ####
		n_new = n_processes - 1
		n_map = params['jknumber']

		kbn = params['kbinNum']
		PK = np.zeros(shape=(n_map,kbn))

		self.q = mp.JoinableQueue()

		if n_new <= 0:
			for ii in range(n_map):
				self.process_map(ii, PK[ii])
		elif n_new > 32:
			raise ValueError('Process limit is 32')
		else:
			process_list = range(n_new)
			for ii in xrange(n_map + n_new):
				if ii >= n_new:
					PK[ii-n_new] = self.q.get()
					process_list[ii%n_new].join()
					if process_list[ii%n_new].exitcode != 0:
						raise RuntimeError("A thred faild with exit code"
							+ str(process_list[ii%n_new].exitcode))
				if ii < n_map:
					process_list[ii%n_new] = mp.Process(
						target=self.process_map, 
						args=(ii, ii%n_new))
					process_list[ii%n_new].start()

		if FKPweight:
			sp.save(params['output_root']+\
				'PKjk_fkp_' + resultf, PK)
		else:
			sp.save(params['output_root']+'PKjk_' + resultf, PK)
		#PKmean = sp.load(params['input_root'] + 'PK.npy')
		PKmean = PK.mean(axis=0)
		PK[:] = (PK[:]-PKmean)**2
		PKvar = np.sum(PK, axis=0)
		PKvar = PKvar*(params['jknumber']-1)/params['jknumber']
		PKvar = np.sqrt(PKvar)
		print PKvar
		if FKPweight:
			sp.save(params['output_root']+\
				'PKvar_fkp_' + resultf, PKvar)
		else:
			sp.save(params['output_root']+\
				'PKvar_' + resultf, PKvar)




	def process_map(self, jknum, rank):
		params = self.params
		mid = params['mid']
		params['mid'] = ('jk'+str(jknum)+mid[0], 'jk'+str(jknum)+mid[1])

		kiyopy.utils.mkparents(params['output_root'])
		inifile = params['output_root']+ 'rank' + str(rank) +'params.ini'
		parse_ini.write_params(params, inifile ,prefix='pk_')
		PK = mkpower.PowerSpectrumMaker(
			inifile, feedback=self.feedback).execute()
		self.q.put_nowait(PK)

#		jkbin = int(params['jknumber']/size)
#		jk = np.array(range(jkbin))
#		jk = (rank*jkbin)+jk
#
#
#		kbn = params['kbinNum']
#		PK = np.zeros(shape=(jk.shape[0],kbn))
#		num = 0
#		for i in jk:
#			params['mid'] = ('jk'+str(i)+mid[0], 'jk'+str(i)+mid[1])
#			#params['mid'] = ('jk'+str(i)+mid[0], mid[1])
#			#print params['mid'][0]
#			kiyopy.utils.mkparents(params['output_root'])
#			inifile = params['output_root']+ 'rank' + str(rank) +'params.ini'
#			parse_ini.write_params(params, inifile ,prefix='pk_')
#			PK[num] = mkpower.PowerSpectrumMaker(
#				inifile, feedback=self.feedback).execute()
#			num = num + 1
#		
#		if rank !=0 :
#			comm.send(PK, dest=0, tag=11)
#		if rank ==0:
#			print 'Calculate the error!'
#			for i in range(1,size):
#				print 'Receive ' + str(i)
#				PK = np.append(PK,comm.recv(source=i, tag=11),axis=0)
#
#			if FKPweight:
#				sp.save(params['output_root']+\
#					'PKjk_fkp_' + resultf, PK)
#			else:
#				sp.save(params['output_root']+'PKjk_' + resultf, PK)
#			#PKmean = sp.load(params['input_root'] + 'PK.npy')
#			PKmean = PK.mean(axis=0)
#			PK[:] = (PK[:]-PKmean)**2
#			PKvar = np.sum(PK, axis=0)
#			PKvar = PKvar*(params['jknumber']-1)/params['jknumber']
#			PKvar = np.sqrt(PKvar)
#			print PKvar
#			if FKPweight:
#				sp.save(params['output_root']+\
#					'PKvar_fkp_' + resultf, PKvar)
#			else:
#				sp.save(params['output_root']+\
#					'PKvar_' + resultf, PKvar)

if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		JackKnifeError(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		JackKnifeError().execute()
