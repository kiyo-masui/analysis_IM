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
from scipy.interpolate import interp1d

import functions

pi = 3.1415926
deg2rad = pi/180.

params_init = {
	# the process number used in this module, now we just support processes=1
	'processes' : 1,

	# indicated the location of the data to be plot
	'input_root' : '../newmaps/',
	# the dirction where the figure will be saved in
	'output_root' : './',

	'resultf': '',

	'jknumber' : 100.,
	# set to import the FKP weighted power spectrum
	'FKPweight' : False,

	# set the parameters for the theoretical power spectrum from CAMB
	'OmegaHI' : 1.e-3,
	'Omegam' : 0.24,
	'OmegaL' : 0.76,
	'z' : 0.7,

	# set the lower limit of the y-axis in the figure
	'ymin' : 1.e3,
	'ymax' : 1.e6,
	# set the upper and lower limit of the x-axis in the figure
	'kmin' : -1,
	'kmax' : -1,
	# set the unit of the power spectrum
	'PKunit' : 'mk',
	# set to plot the cross power spectrum. 
	# Then the unit will be changed to square root of the auto power spectrum
	'cross' : False,

	'optical' : False,
	# set to calibrate with transfer function
	'biascal' : False,
	# set the location of the transfer function file
	'biascaldir_k' : '', 
	'biascaldir_b' : '',
	# set to plot Delta instead of Pk
	'delta_plot' : True
}
prefix = 'pt_'

class PowerSpectrumPlot(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict,
			params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

	def execute(self, nprocesses=1):
		params = self.params

		resultf = functions.getresultf(params)

		FKPweight = params['FKPweight']
		kmin = params['kmin']
		kmax = params['kmax']

		# import the power spectrum result
		try:
			k  = sp.load(params['input_root']+resultf+'_k.npy')
			PK = sp.load(params['input_root']+resultf+'_p.npy')
		except IOError:
			try:
				k  = sp.load(params['input_root']+resultf+'_k_combined.npy')
				PK = sp.load(params['input_root']+resultf+'_p_combined.npy')
			except IOError:
				print 'Error: NO power to be plotted!!'

		dk = sqrt(k[1]/k[0])
		PK[PK<0]=0

		non0 = PK.nonzero()

		#print k
		PK = PK.take(non0)[0]
		k = k.take(non0)[0]
		#print k

		## show each power spectrum
		#self.powershow_all()
		#return 0


		# import the jk error
		#PKjk = sp.load(params['input_root'] + 'PKjk_'+resultf+'.npy')
		JKplot = True
		try:
			PKvar = sp.load(params['input_root']+resultf+'_p_var_jk.npy')
			PKvar = PKvar
			PKerr = np.ndarray(shape=(2,len(non0[0])))
			PKerr[0] = PKvar.take(non0)[0]
			PKerr[1] = PKvar.take(non0)[0]

			for i in range(len(PKerr[0])):
				if PKerr[0][i] >=PK[i]:
					PKerr[0][i] = PK[i]-1.e-10
		except IOError:
			print '\t::No Jeckknife Error!!'
			JKplot = False

		# import the average error
		Averageplot = True
		try:
			PKvar = sp.load(params['input_root']+resultf+'_p_var_combined.npy')
			PKvar = PKvar #*2
			PKerr_average = np.ndarray(shape=(2,len(non0[0])))
			PKerr_average[0] = PKvar.take(non0)[0]
			PKerr_average[1] = PKvar.take(non0)[0]

			for i in range(len(PKerr_average[0])):
				if PKerr_average[0][i] >=PK[i]:
					PKerr_average[0][i] = PK[i]-1.e-15
		except IOError:
			print '\t::No Average Error!!'
			Averageplot = False

		# import the fkp result
		if FKPweight:
			kfkp 	  = sp.load(params['input_root']+resultf+'_k_fkp.npy')
			PKfkp   = sp.load(params['input_root']+resultf+'_p_fkp.npy')
			PKjkfkp = sp.load(params['input_root']+resultf+'_p_jk_fkp.npy')
			PKvarfkp= sp.load(params['input_root']+resultf+'_p_var_jk_fkp.npy')
			non0fkp = PKfkp.nonzero()
			kfkp = kfkp.take(non0fkp)[0]
			PKfkp = PKfkp.take(non0fkp)[0]
			PKerrfkp = np.ndarray(shape=(2,len(non0fkp[0])))
			PKerrfkp[0] = PKvarfkp.take(non0fkp)[0]
			PKerrfkp[1] = PKvarfkp.take(non0fkp)[0]
			for i in range(len(PKerrfkp[0])):
				if PKerrfkp[0][i]>= PKfkp[i]:
					PKerrfkp[0][i] = PKfkp[i]-1.e-10

		# import the theoretical power spectrum 
		PKcamb = functions.getpower_th(params)

		klow = k
		k = k*dk
		kup = k*dk
		kerr = np.ndarray(shape=(2, len(k)))
		kerr[0] = k-klow
		kerr[1] = kup-k
		print k
		
		# plot
		plt.figure(figsize=(10,5))
		
		# change variance into power spectrum
		if params['delta_plot']==False:
			PK = PK/(k*k*k/2./pi/pi)
			if JKplot:
				PKerr = PKerr/(k*k*k/2./pi/pi)
			if Averageplot:
				PKerr_average = PKerr_average/(k*k*k/2./pi/pi)
			if FKPweight:
				PKfkp = PKfkp/(kfkp*kfkp*kfkp/2./pi/pi)
				PKerrfkp = PKerrfkp/(kfkp*kfkp*kfkp/2./pi/pi)
			PKcamb[1] = PKcamb[1]/(PKcamb[0]*PKcamb[0]*PKcamb[0]/2./pi/pi)
		
		plt.subplot('111')
		plt.plot(PKcamb[0], PKcamb[1], 'g-', linewidth=2, 
			label='Theoretical Power Spectrum')
		#plt.plot(PKnonl_k, PKnonl, 'k-', linewidth=2,
		#	label='Theoretical Convolved wiht Window')
		if JKplot:
			#plt.errorbar(k, PK, PKerr, fmt='o', c='r', 
			plt.errorbar(k, PK, PKerr, kerr, fmt='o', c='r', 
				label='WiggleZ Auto Power Spectrum',
				capsize=4.5, elinewidth=1)
		elif Averageplot:
			if params['cross']:
				#plt.errorbar(k, PK, PKerr_average, fmt='o', c='r',
				plt.errorbar(k, PK, PKerr_average, kerr, fmt='o', c='r',
					#label='Noise Inv Weight', 
					label='21cm and WiggleZ Cross Power Spectrum',
					capsize=4.5, elinewidth=1)
			elif params['optical']:
				#plt.errorbar(k, PK, PKerr_average, fmt='o', c='r',
				plt.errorbar(k, PK, PKerr_average, kerr, fmt='o', c='r',
					#label='Noise Inv Weight', 
					label='WiggleZ Auto Power Spectrum',
					capsize=4.5, elinewidth=1)
			else:
				#plt.errorbar(k, PK, PKerr_average, fmt='o', c='r',
				plt.errorbar(k, PK, PKerr_average, kerr, fmt='o', c='r',
					#label='Noise Inv Weight', 
					label='21cm Auto Power Spectrum',
					capsize=4.5, elinewidth=1)
				
		else:
			plt.scatter(k, PK, c='b', label='Noise Inv Weight')

		if FKPweight:
			plt.errorbar(kfkp, PKfkp, PKerrfkp, fmt='o', c='r',
				label='FPK Weight', capsize=4)

		plt.loglog()
		#plt.semilogy()
		ymin = params['ymin']
		ymax = params['ymax']
		if params['delta_plot']==False:
			ymin = ymin*1.e3
			ymax = ymax*1.e3
		plt.ylim(ymin=ymin, ymax=ymax)	
		xmin = 0.01
		xmax = 0.90
		if (kmin!=-1):
			xmin = kmin-0.01*kmin
		if (kmax!=-1):
			xmax = kmax+0.01*kmax
		plt.xlim(xmin=xmin, xmax=xmax)
		plt.xlabel('$k (hMpc^{-1})$')
		if params['cross']:
			if params['delta_plot']:
				ylabel = '$\Delta_k^2 (%(PKunit)s)$' %params
			else:
				ylabel = '$P_k (%(PKunit)s(h^{-1}Mpc)^3)$' %params
		elif params['optical']:
			if params['delta_plot']:
				ylabel = '$\Delta_k^2$' 
			else:
				ylabel = '$P_k (h^{-1}Mpc)^3$' %params
		else:
			if params['delta_plot']:
				ylabel = '$\Delta_k^2 (%(PKunit)s^{2})$' %params
			else:
				ylabel = '$P_k (%(PKunit)s^{2}(h^{-1}Mpc)^3)$' %params
		plt.ylabel(ylabel)
		#plt.xticks(np.arange(kmin, kmax, (kmax-kmin)/5.), 
		#	np.arange(kmin, kmax, (kmax-kmin)/5.))
		plt.tick_params(length=6, width=1.)
		plt.tick_params(which='minor', length=3, width=1.)
		plt.legend(loc=0, scatterpoints=1)

		plt.savefig(params['output_root']+resultf+'_p_err.eps', format='eps')
		plt.savefig(params['output_root']+resultf+'_p_err.png', format='png')

		plt.show()

	def powershow_all(self):
		params = self.params

		resultf = functions.getresultf(params)

		if params['resultf']=='':
			k   = sp.load(params['input_root']+resultf+'_k.npy')
		else:
			k   = sp.load(params['input_root']+resultf+'_k_combined.npy')
			
		PKeach = sp.load(params['input_root']+resultf+'_p_each.npy')

		for ii in range(PKeach.shape[0]):
			plt.scatter(k, PKeach[ii])
			#plt.plot(k, PKeach[ii])
		plt.loglog()
		#plt.semilogx()
		ymin = params['ymin']
		plt.ylim(ymin=ymin)	
		#plt.xlim(xmin=0.01, xmax=0.9)
		plt.xlim(xmin=k.min()-0.1*k.min(), xmax=k.max()+0.1*k.max())
		#plt.show()

		return 0



if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
