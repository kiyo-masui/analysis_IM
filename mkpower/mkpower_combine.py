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
#from mpi4py import MPI


pi = 3.1415926
deg2rad = pi/180.

params_init = {
	'processes' : 1,
	'plot' : False,
	'input_root' : '../../../jkmap/map/',
	'hrlist' : (),
	'ltlist' : (),
	'resultf' : '',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : '../../../powerresult/',

	'cldir' : '',
	'cllist' : (),

	'boxshape' : (128,128,128),
	'boxunit' : 15., # in unit Mpc
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-64*15,64*15),
	'Zrange' : (0.,128*15),
	'kbinNum' : 20,
	'kmin' : -1,
	'kmax' : -1,

	'FKPweight' : False,
	'FKPpk' : 1.e-3,
	'jkerror' : False,
	'sme' : True,
}
prefix = 'pkc_'


class PowerSpectrumMaker(mkpower.PowerSpectrumMaker):
	"""Calculate The Power Spectrum"""
	B = []
	Bk = []
	q = mp.JoinableQueue()

	def __init__(self, parameter_file_or_dict=None, feedback=1):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, 
			params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback
	
	def execute(self, nprocesses=1):
		
		#comm = MPI.COMM_WORLD
		#rank = comm.Get_rank()
		#size = comm.Get_size()

		params = self.params
		resultf = params['resultf']

		# Make parent directory and write parameter file.
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']

		FKPweight = params['FKPweight']
		n_processes = params['processes']

		if params['cllist']!=():
			# Read in the bias calibration data
			#self.B = sp.load(params['cldir']+'b_each_bias.npy')
			self.B = sp.load(params['cldir']+'b_bias.npy')
			self.Bk = sp.load(params['cldir']+'k_bias.npy')
			#print self.B


		#### Process ####
		n_new = n_processes - 1
		n_map = len(params['hrlist'])

		kbn = params['kbinNum']
		kmin = params['kmin']
		kmax = params['kmax']
		PK = np.zeros(shape=(n_map,kbn))

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
			if params['jkerror']:
				sp.save(params['output_root']+resultf+'_p_each_jk_fkp', PK)
			else:
				sp.save(params['output_root']+resultf+'_p_each_fkp', PK)
		else:
			if params['jkerror']:
				sp.save(params['output_root']+resultf+'_p_each_jk', PK)
			else:
				sp.save(params['output_root']+resultf+'_p_each', PK)
		#PKmean = sp.load(params['input_root'] + 'PK.npy')
		print PK
		PKmean = PK.mean(axis=0)
		PKvar = PK.std(axis=0)
		if params['jkerror']:
			print '\tJackKnife Error:'
			PKvar = PKvar*sqrt((n_map-1))
		elif params['sme']:
			print '\tUsing the Standard Deviation of Sample Mean as error'
			PKvar = PKvar/sqrt(n_map)
		print PKvar
			

		#print 
		#print '===power before cal==='
		#print PKmean
		#print PKvar
		if params['cllist']!=():
			#print 
			#print '=== cal==='
			#print self.B
			PKmean = PKmean*self.B
			PKvar  = PKvar*self.B
		#print 
		#print '===power after cal==='
		#print PKmean
		#print PKvar
		print 
		print '===power after cal==='
		print PKmean
		print PKvar

		kunit = 2.*pi/(params['boxunit'])
		if (kmin==-1) or (kmax==-1):
			k = np.logspace(
				log10(1./params['boxshape'][0]), log10(sqrt(3)), num=kbn+1)
		else:
			kmin = kmin/kunit
			kmax = kmax/kunit
			k = np.logspace(log10(kmin), log10(kmax), num=kbn+1)
		k = k*2.*pi/params['boxunit']
		k = k[:-1]
		sp.save(params['output_root']+resultf+'_k_combined', k)
		if FKPweight:
			if params['jkerror']:
				sp.save(params['output_root']+resultf+'_p_var_jk_fkp', PKvar)
			else:
				sp.save(params['output_root']+resultf+'_p_var_combined_fkp', PKvar)
				sp.save(params['output_root']+resultf+'_p_combined_fkp', PKmean)
		else:
			if params['jkerror']:
				sp.save(params['output_root']+resultf+'_p_var_jk', PKvar)
			else:
				sp.save(params['output_root']+resultf+'_p_var_combined', PKvar)
				sp.save(params['output_root']+resultf+'_p_combined', PKmean)
	
#	def findbias(self, hr, last, cllist, B):
#
#		def findidx(aa):
#			idx=('A', 'B', 'C', 'D')
#			for s in idx:
#				if aa.find(s)!=-1:
#					return s
#			return 'NO'
#
#		bias = np.ones(B.shape[1])
#
#		for i in range(len(hr)):
#			idx = findidx(hr[i]) + findidx(last[i])
#			print idx
#			if idx=='NONO':
#				print 'no bias found'
#				continue
#			clidx = cllist.index(idx)
#			bias = bias*np.sqrt(B[clidx])
#			#print bias
#		return bias


	def process_map(self, mapnum, rank):
		params = self.params
		params['hr'] = (params['hrlist'][mapnum][0],params['hrlist'][mapnum][1])
		params['last'] =(params['ltlist'][mapnum][0],params['ltlist'][mapnum][1])

		PK, k, PK2, k2 = self.GetPower()
		print PK

#		if params['cllist']!=():
#			bias = self.findbias(params['hr'], 
#				params['last'], params['cllist'], self.B)
#
#		if params['cllist']!=():
#			PK = PK*bias

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
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
