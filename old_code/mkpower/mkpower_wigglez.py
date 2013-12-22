#! /usr/bin/env python

import scipy as sp
import numpy as np
#from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp

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

import mkpower_combine

import functions

params_init = {
	'processes' : 1,
	'plot' : True,
	'saveweight' : False,
	'input_root' : '../newmaps/',
	'hrlist' : (),
	'ltlist' : (),
	'resultf' : '',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : './',

	'cldir' : '',
	'cllist' : (),

	'boxshape' : (60,12,6),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),

	'kbinNum' : 20,
	'kmin' : -1,
	'kmax' : -1,

	'FKPweight' : False,
	'FKPpk' : 1.e-3,
	'jkerror' : False,
	'OmegaHI' : 1.e-3,
	'Omegam' : 0.23,
	'OmegaL' : 0.74,
	'z' : 1,
	'sme' : True,
}

pi = 3.1415926
deg2rad = pi/180.
prefix = 'wgz_'

#q = mp.JoinableQueue()

class PowerSpectrumMaker(mkpower_combine.PowerSpectrumMaker):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=1):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict,
			params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def process_map(self, mapnum, rank):
		"""
			We rewrite this function to add short noise term
		"""
		params = self.params
		params['hr'] = (params['hrlist'][mapnum][0],params['hrlist'][mapnum][1])
		params['last'] =(params['ltlist'][mapnum][0],params['ltlist'][mapnum][1])

		if params['cllist']!=():
			bias = self.findbias(params['hr'], 
				params['last'], params['cllist'], self.B)

		PK, k, PK2, k2 = self.GetPower()

		resultf = functions.getresultf(params)
		shortnoise_fname = \
			params['input_root'] + resultf + '_shortnoise_p_combined.npy'
		shortnoise = np.load(shortnoise_fname)
		PK = PK - shortnoise
		#print PK
		if params['cllist']!=():
			PK = PK*bias
		self.q.put_nowait(PK)
	
#	def GetDelta(self, box, nbox, mbox):
#		params = self.params
#		V = params['boxunit']**3
#
#		alpha = box.flatten().sum()/mbox.flatten().sum()
#		#print alpha
#		box = box - nbox
#		box[nbox!=0] = box[nbox!=0]/nbox[nbox!=0]
#
#		#print nbox.flatten().max()
#		#print nbox.flatten().min()
#		#print box.flatten().max()
#		#print box.flatten().min()
#		#print box[100][20]
#		#print box.flatten().max()
#
#		self.shortnoise = \
#			(1.+alpha)*(nbox).flatten().sum()/(nbox**2/V).flatten().sum()
#		print self.shortnoise
#
#		# get the fkp weight
#		#nbox = nbox/(1.+nbox*1000)
#
#		return box, nbox


if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
