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

params_init = {
	'processes' : 1,
	'plot' : True,
	'saveweight' : False,
	'input_root' : '../newmaps/',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : ('dirty_map_',),
	'polarizations' : ('I',),
	'last' : (),
	'output_root' : './',

	'boxshape' : (60,12,6),
	'boxunit' : 15., # in unit Mpc h-1
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-6*15,6*15),
	'Zrange' : (0.,6*15),

	'kbinNum' : 20,
	'kmin' : None,
	'kmax' : None,

	'FKPweight' : False,
	'FKPpk' : 1.e-3,
	'OmegaHI' : 1.e-3,
	'Omegam' : 0.23,
	'OmegaL' : 0.74,
	'z' : 1,
}

pi = 3.1415926
deg2rad = pi/180.
prefix = 'pk_'

class PowerSpectrumMaker(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict,
			params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):

		params = self.params
		out_root = params['output_root']
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]
		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		fkpp = params['FKPpk']
		FKPweight = params['FKPweight']


		fftbox = self.GetRadioFFTbox()


		kbn = params['kbinNum']
		kmin = params['kmin']
		kmax = params['kmax']

		kunit = 2.*pi/(params['boxunit'])
		PK = np.zeros(kbn)
		if (kmin==None) or (kmax==None):
			k = np.logspace(
				log10(1./params['boxshape'][0]), log10(sqrt(3)), num=kbn+1)
		else:
			kmin = kmin/kunit
			kmax = kmax/kunit
			k = np.logspace(log10(kmin), log10(kmax), num=kbn+1)
		PK2 = np.zeros(shape=(10, 10))
		k2 = np.zeros(shape=(2, 10))
		MakePower.Make(fftbox, PK, k, PK2, k2)
		k = k*kunit
		k2 = k2*kunit
		boxN = params['boxshape'][0]*params['boxshape'][1]*params['boxshape'][2]
		PK = PK
		PK2 = PK2
		#PK = PK/(V*boxN)
		#PK2 = PK2/V/boxN

#		OmegaHI = params['OmegaHI']
#		Omegam = params['Omegam']
#		OmegaL = params['OmegaL']
#		z = params['z']
#		a3 = (1+z)**(-3)
#		Tb = (OmegaHI) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
#			* ((1.+z)/2.5)**0.5
#		#if params['PKunit']=='mK':
#		#	Tb = Tb/1.e-3
#
#		xx = ((1.0/Omegam)-1.0)/(1.0+z)**3
#		num = 1.0 + 1.175*xx + 0.3046*xx**2 + 0.005335*xx**3
#		den = 1.0 + 1.875*xx + 1.021 *xx**2 + 0.1530  *xx**3
#
#		G = (1.0 + xx)**0.5/(1.0+z)*num/den
#		#print G**2*(1+z)**(-2)
#		PKrand = sp.load('/home/lixiating/Work/GBT/crandomresult/PK_combined_auto-100randmaps.npy')
#		PKrand = PKrand*Tb**2
#
#		PK = PKrand


		non0 = PK.nonzero()
		#print PK.min(), PK.max()
		#print PK2
		#print PK
		#print k[non0]
		#return 0
		if params['saveweight']:
			if FKPweight:
				sp.save(out_root+'PK_fkp_'+resultf, PK)
				sp.save(out_root+'PK2_fkp_'+resultf, PK2)
				sp.save(out_root+'k_fkp_'+resultf, k)
			else:
				sp.save(out_root+'PK_'+resultf, PK)
				sp.save(out_root+'PK2_'+resultf, PK2)
				sp.save(out_root+'k_'+resultf, k)

		power_th = np.load('../maps/test/power_yuebin.npy')
		power_th[0] = power_th[0]/0.705
		power_th[1] = power_th[1]*0.705**3

		#power_0 = np.load(out_root+'PK_test_yuebin_2-test_yuebin_2_0.npy')
		#k_0 = np.load(out_root+'k_test_yuebin_2-test_yuebin_2_0.npy')
		#non0_0 = power_0.nonzero()
		#power_0 = power_0.take(non0_0)
		#k_0 = k_0.take(non0_0)

		if self.plot==True:
			print PK
			#print power_0
			plt.figure(figsize=(8,8))
			#print k
			#print PK
			plt.subplot('211')
			plt.scatter(k.take(non0), PK.take(non0))
			#plt.scatter(k_0, power_0, c='w')
			plt.plot(power_th[0], power_th[1], c='r')
			plt.loglog()
			#plt.semilogx()
			plt.ylim(ymin=1.e-1)	
			plt.xlim(xmin=k.min(), xmax=k.max())
			plt.title('Power Spectrum')
			plt.xlabel('$k$')
			plt.ylabel('$P(k) (Kelvin^{2}(h^{-1}Mpc)^3)$')

			PK = PK*k*k*k/2./pi/pi
			power_th[1] = power_th[1]*power_th[0]*power_th[0]*power_th[0]/2/pi/pi

			#print PK
			plt.subplot('212')
			plt.scatter(k.take(non0), PK.take(non0))
			plt.plot(power_th[0], power_th[1], c='r')
			plt.loglog()
			#plt.semilogx()
			plt.ylim(ymin=1.e-4)	
			plt.xlim(xmin=k.min(), xmax=k.max())
			plt.xlabel('$k (h Mpc^{-1})$')
			plt.ylabel('$\Delta^2 (Kelvin^{2})$')
			#plt.show()
			if FKPweight:
				plt.savefig(out_root+'power_fkp_'+resultf+'.eps', format='eps')
			else:
				plt.savefig(out_root+'power_'+resultf+'.eps', format='eps')

	#		#PK2 = np.log10(PK2)
	#		plt.figure(figsize=(6,6))
	#		extent = (k2[0][0], k2[0][-1], k2[1][0], k2[1][-1])
	#		plt.imshow(PK2, origin='lower', extent = extent, interpolation='nearest')
	#		plt.xlabel('$k vertical (h Mpc^{-1})$')
	#		plt.ylabel('$k parallel (h Mpc^{-1})$')
	#		cb = plt.colorbar()
	#		cb.set_label('$lg(P^{2D}_{k_pk_v}) (Kelvin^2(h^{-1}Mpc)^3)$')
	#		plt.loglog()
	#		if FKPweight:
	#			plt.savefig(out_root+'power2_fkp_'+resultf+'.eps', format='eps')
	#		else:
	#			plt.savefig(out_root+'power2_'+resultf+'.eps', format='eps')

			plt.show()
			#print 'Finished @_@ '
		return PK

	def GetRadioFFTbox(self):
		params = self.params
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, 
			params['output_root']+'params.ini',prefix='pk_' )
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']
		all_out_fname_list = []
		all_in_fname_list = []
		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		fkpp = params['FKPpk']
		FKPweight = params['FKPweight']
		
		#### Process ####
		pol_str = params['polarizations'][0]
		hr_str = params['hr'][0]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][0]
		box_fname = in_root + 'fftbox_' + hr_str + mid[0] + end + '.npy'
		box = np.load(box_fname)

		nbox_fname = in_root + 'fftbox_' + hr_str + mid[1] + end + '.npy'
		nbox = np.load(nbox_fname)

		#Using map in different day 
		hr_str = params['hr'][1]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][1]
		box_fname = in_root + 'fftbox_' + hr_str + mid[0] + end + '.npy'
		box2 = np.load(box_fname)

		nbox_fname = in_root + 'fftbox_' + hr_str + mid[1] + end + '.npy'
		nbox2 = np.load(nbox_fname)

		normal = (nbox**2).flatten().sum()
		normal2 = (nbox2**2).flatten().sum()
		normal = sqrt(normal)*sqrt(normal2)
		#print normal
		box = box*nbox
		box2 = box2*nbox2
		#box = nbox*nbox
		#box2 = nbox2*nbox2

		V = params['boxunit']**3

		print "PowerMaker: FFTing "
		inputa = np.zeros(params['boxshape'], dtype=complex)
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		inputa.imag = 0.
		inputa.real = box
		FFTW.execute(fft)

		#print outputa[10][10]

		inputb = np.zeros(params['boxshape'], dtype=complex)
		outputb = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputb,outputb, direction='forward', flags=['measure'])
		inputb.imag = 0.
		inputb.real = box2
		FFTW.execute(fft)
		
		#print outputb[10][10]

		fftbox = (outputa.__mul__(outputb.conjugate())).real

		#fftbox = (outputa*(outputa.conjugate())).real
		fftbox = fftbox*V/normal #/2./pi/pi/pi
		#boxN = params['boxshape'][0]*params['boxshape'][1]*params['boxshape'][2]
		#fftbox = fftbox*V/boxN #/2./pi/pi/pi

		return fftbox

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
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
