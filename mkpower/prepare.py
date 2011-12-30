#! /usr/bin/env python
"""
	This module used to convert the maps into fftbox
	
	----  Dec 11, 2011. Y.-CH. LI  ----
	If the intensity map is already overdensity, the parameter 'mid' 
	should have two elements:
		middle part for intensity map name, and 
		middle part for inverse noise map name.
	The inverse noise map will be the weight for power spectrum estimation.

	If the intensity map is not the overdensity, the parameter 'mid'
	should have three elements:
		middle part for intensity map name, 
		middle part for inverse noise map name(or selection function), and
		middle part for mock map name.
		
	
"""
import os
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

# public functions are defined in this module
import functions

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

	'cutlist': [],
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
		#parse_ini.write_params(params, 
		#	params['output_root']+'params.ini',prefix='pk_')
		hr = params['hr']
		mid = params['mid']
		last = params['last']
		pol_str = params['polarizations'][0]
		n_processes = params['processes']
		out_root = params['output_root']

		# make directions for fftbox saving
		if os.path.exists(out_root+'fftbox/')==False:
			os.mkdir(out_root+'fftbox/')
		
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

					# prepare the mock maps for optical 
					if len(mid)==3:
						mock_fname = hr[ii] + mid[2] + end + '.npy'
						process_list[ii%n_new] = mp.Process(
							target=self.process_map, 
							args=(imap_fname, nmap_fname, mock_fname))
					else:
						process_list[ii%n_new] = mp.Process(
							target=self.process_map, 
							args=(imap_fname, nmap_fname))
					process_list[ii%n_new].start()
		return 0

	def process_map(self, imap_fname, nmap_fname, mock_fname=None):
		params = self.params
		out_root = params['output_root']
		in_root = params['input_root']
		
		if mock_fname != None:
			imap, nmap, mmap = functions.getmap(
				in_root+imap_fname, in_root+nmap_fname, in_root+mock_fname)
			box , nbox, mbox = functions.fill(params, imap, nmap, mmap)
			#print mmap.flatten().mean()
			#print nmap.flatten().max()
			#print box.flatten().max()

			## get the overdencity
			#alpha = box.flatten().sum()/mbox.flatten().sum()
			#print alpha
			#box = box - alpha*mbox
			#box[nbox!=0] = box[nbox!=0]/nbox[nbox!=0]

			## get the fkp weight
			#nbox = nbox/(1.+nbox*1000)

			pkrm_mfname = out_root + 'fftbox/' + 'fftbox_' +  mock_fname
			np.save(pkrm_mfname, mbox)
			#algebra.save(pkrm_nfname, mbox)
		else:
			#box, nbox = self.fill(imap, nmap)
			imap, nmap = functions.getmap(in_root+imap_fname, in_root+nmap_fname)

			# subtract the mean value of the imaps
			#print "--The mean value for imap is:",imap.flatten().mean(),"--"
			imap = imap - imap.flatten().mean()

			# cut off some bad frequencies. 
			if len(params['cutlist'])!=0:
				print '\t:Bad frequencies cutting off'
				nmap[params['cutlist']]=0

			box , nbox = functions.fill(params, imap, nmap)

			#nbox = nbox/(1.+3.e3*nbox)
			#nbox = nbox/(1.+7.e-5*nbox)

		pkrm_fname = out_root + 'fftbox/' + 'fftbox_' + imap_fname
		np.save(pkrm_fname, box)
		#algebra.save(pkrm_fname, box)

		pkrm_nfname = out_root + 'fftbox/' + 'fftbox_' +  nmap_fname
		np.save(pkrm_nfname, nbox)
		#algebra.save(pkrm_nfname, nbox)

#		print mmap.shape

#		print 'Removing Peak... Map:' + hr_str[:-1]
#		self.pkrm(imap,nmap,
#			out_root+'pkrm'+hr_str+'dirty_map_'+pol_str+'.png', threshold=2.5)

		#print imap.flatten().max()
		#print imap.flatten().min()
		#print box.flatten().max()
		#print box.flatten().min()


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


if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		Prepare(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		Prepare().execute()
