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

	'OmegaHI' : 1.e-3,
	'Omegam' : 0.23,
	'OmegaL' : 0.74,
}
prefix = 'fpk_'

class FKPPowerSpectrumMaker(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params
		resultf = params['hr'][0]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][0]
		resultf = resultf + '-' + params['hr'][1]
		if len(params['last']) != 0:
			resultf = resultf + params['last'][1]

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='pk_' )
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']
		all_out_fname_list = []
		all_in_fname_list = []
		
		#### Process ####
		pol_str = params['polarizations'][0]
		hr_str = params['hr'][0]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][0]
		imap_fname = in_root + hr_str + mid[0] + end + '.npy'
		imap = algebra.load(imap_fname)
		imap = algebra.make_vect(imap)
		if imap.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		nmap_fname = in_root + hr_str + mid[1] + end + '.npy'
		try:
			nmap = algebra.load(nmap_fname)
			nmap = algebra.make_vect(nmap)

			bad = nmap<1.e-5*nmap.flatten().max()
			nmap[bad] = 0.
			#non0 = nmap.nonzero()
			#imap[non0] = imap[non0]/nmap[non0]
		except IOError:
			print 'NO Noise File :: Set Noise to One'
			nmap = algebra.info_array(sp.ones(imap.shape))
			nmap.axes = imap.axes
			nmap = algebra.make_vect(nmap)



		##  Using map in different day 
		hr_str = params['hr'][1]
		end = pol_str
		if len(params['last']) != 0:
			end = end + params['last'][1]
		imap_fname = in_root + hr_str + mid[0] + end + '.npy'
		imap2 = algebra.load(imap_fname)
		imap2 = algebra.make_vect(imap2)
		if imap2.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		nmap_fname = in_root + hr_str + mid[1] + end + '.npy'
		try:
			nmap2 = algebra.load(nmap_fname)
			nmap2 = algebra.make_vect(nmap2)
	
			bad = nmap2<1.e-5*nmap2.flatten().max()
			nmap2[bad] = 0.
			#non0 = nmap2.nonzero()
			#imap2[non0] = imap2[non0]/nmap2[non0]
		except IOError:
			print 'NO Noise File :: Set Noise to One'
			nmap2 = algebra.info_array(sp.ones(imap2.shape))
			nmap2.axes = imap2.axes
			nmap2 = algebra.make_vect(nmap2)

		## FKP ##
		fkp = algebra.info_array(sp.ones(imap.shape))
		fkp.axes = imap.axes
		fkp = algebra.make_vect(fkp)

		fkp2 = algebra.info_array(sp.ones(imap.shape))
		fkp2.axes = imap.axes
		fkp2 = algebra.make_vect(fkp2)

		OmegaHI = params['OmegaHI']
		Omegam = params['Omegam']
		OmegaL = params['OmegaL']
		freq = imap.get_axis('freq')
		z = 1.4e9/freq -1.
		a3 = (1+z)**(-3)
		Tb = 0.3e-3 * (OmegaHI/1.e-3) * ((Omegam + a3*OmegaL)/0.29)**(-0.5)\
			* ((1.+z)/2.5)**0.5

		for i in range(fkp.shape[0]):
			fkp[i] = fkp[i]*Tb[i]/(1.+Tb[i]*1.e4)
			fkp2[i] = fkp2[i]*Tb[i]/(1.+Tb[i]*1.e4)
		
			

		#print imap2.flatten().min()
		#print dmap2.flatten().min()
		#print nmap2.flatten().sum()

		mapshape = np.array(imap.shape)
		#print imap.shape

		r  = self.discrete(self.fq2r(imap.get_axis('freq')))
		ra = self.discrete(imap.get_axis('ra'))*deg2rad
		de = self.discrete(imap.get_axis('dec'))*deg2rad
		ra0= ra[int(ra.shape[0]/2)]
		ra = ra - ra0
		dr = r.ptp()/r.shape[0]
		dra= ra.ptp()/ra.shape[0]
		dde= de.ptp()/de.shape[0]
		disc_n = params['discrete']

		#print r.min(), r.max()
		#print self.xyz(ra.min(), de.min(), r.min())
		#print self.xyz(ra.max(), de.min(), r.min())
		#print self.xyz(ra.min(), de.max(), r.min())
		#print self.xyz(ra.max(), de.max(), r.min())
		#print self.xyz(ra.min(), de.min(), r.max())
		#print self.xyz(ra.max(), de.min(), r.max())
		#print self.xyz(ra.min(), de.max(), r.max())
		#print self.xyz(ra.max(), de.max(), r.max())

		#return 0

		mapinf = [dr, dra, dde, disc_n]
		mapinf = np.array(mapinf)

		box = algebra.info_array(sp.zeros(params['boxshape']))
		box.axes = ('x','y','z')
		box = algebra.make_vect(box)
		boxshape = np.array(box.shape)

		box2 = algebra.info_array(sp.zeros(params['boxshape']))
		box2.axes = ('x','y','z')
		box2 = algebra.make_vect(box2)
		
		xrange0 = params['Xrange'][0]
		yrange0 = params['Yrange'][0]
		zrange0 = params['Zrange'][0]
		boxunit = params['boxunit']
		shapex = params['boxshape'][2]
		shapera = ra.shape[0]

		boxinf = [xrange0, yrange0, zrange0, boxunit]
		boxinf = np.array(boxinf)

		print "PowerMaker: Filling the BOX"
		MakePower.Filling(imap, imap2, box, box2, r, ra, de, boxinf, mapinf)
		#MakePower.Filling(dmap, dmap2, box, box2, r, ra, de, boxinf, mapinf)

		# normalize 
		nbox = algebra.info_array(sp.zeros(params['boxshape']))
		nbox.axes = ('x','y','z')
		nbox = algebra.make_vect(nbox)
		nbox2 = algebra.info_array(sp.zeros(params['boxshape']))
		nbox2.axes = ('x','y','z')
		nbox2 = algebra.make_vect(nbox2)

		#non0 = nmap.nonzero()
		#nmap[non0] = np.sqrt(1./nmap[non0])
		#non0 = nmap2.nonzero()
		#nmap2[non0] = np.sqrt(1./nmap2[non0])

		MakePower.Filling(nmap, nmap2, nbox, nbox2, r, ra, de, boxinf, mapinf)
		if params['saveweight']:
			sp.save(out_root+'Weight_'+resultf, nbox)
			sp.save(out_root+'Weight2_'+resultf, nbox2)
			print '\t::Weight Saved'

		## FKP box ##
		fkpbox = algebra.info_array(sp.zeros(params['boxshape']))
		fkpbox.axes = ('x','y','z')
		fkpbox = algebra.make_vect(fkpbox)
		fkpbox2 = algebra.info_array(sp.zeros(params['boxshape']))
		fkpbox2.axes = ('x','y','z')
		fkpbox2 = algebra.make_vect(fkpbox2)
		MakePower.Filling(fkp, fkp2, fkpbox, fkpbox2, r, ra, de, boxinf, mapinf)
		fkpnormal = (fkpbox**2).flatten().sum()

		#nbox = np.ones(params['boxshape'])
		#nbox2 = np.ones(params['boxshape'])

		#bad = nbox<1.e-5*nbox.flatten().max()
		#nbox[bad] = 0.
		#non0 = nbox.nonzero()
		#nbox[non0] = 1./nbox[non0]
		normal = (nbox**2).flatten().sum()
		#bad = nbox2<1.e-5*nbox2.flatten().max()
		#nbox2[bad] = 0.
		#non0 = nbox2.nonzero()
		#nbox2[non0] = 1./nbox2[non0]
		normal2 = (nbox2**2).flatten().sum()
		normal = sqrt(normal)*sqrt(normal2)
		#print normal

		#box = box*nbox
		#box2 = box2*nbox2

		box = box*fkpbox
		box2 = box2*fkpbox2

		print "PowerMaker: FFTing "
		inputa = np.zeros(params['boxshape'], dtype=complex)
		inputa.real = box.copy()
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		FFTW.execute(fft)
		box = outputa.real**2 + outputa.imag**2
		inputa = np.zeros(params['boxshape'], dtype=complex)
		inputa.real = box2.copy()
		outputa = np.zeros(params['boxshape'], dtype=complex)
		fft = FFTW.Plan(inputa,outputa, direction='forward', flags=['measure'])
		FFTW.execute(fft)
		box2 = outputa.real**2 + outputa.imag**2
		fftbox = box.__pow__(0.5)*box2.__pow__(0.5)

		V = params['boxunit']**3

		#fftbox = fftbox*V*V/normal #/2./pi/pi/pi

		fftbox = fftbox*V*V/fkpnormal #/2./pi/pi/pi
      
		#boxN = params['boxshape'][0]*params['boxshape'][1]*params['boxshape'][2]
		#fftbox = fftbox*V*V/boxN #/2./pi/pi/pi

		PK = np.zeros(40)
		k = np.zeros(40)
		PK2 = np.zeros(shape=(10, 10))
		k2 = np.zeros(shape=(2, 10))
		MakePower.Make(fftbox, PK, k, PK2, k2)
		kunit = 2.*pi/(boxunit)
		k = k*kunit
		k2 = k2*kunit
		PK = PK/V
		PK2 = PK2/V
		#PK = PK/(V*boxN)
		#PK2 = PK2/V/boxN

		non0 = PK.nonzero()
		#print PK2
		#print PK
		#print k[non0]
		#return 0
		sp.save(out_root+'fkpPK_'+resultf, PK)
		sp.save(out_root+'fkpPK2_'+resultf, PK2)
		sp.save(out_root+'fkpk_'+resultf, k)

		if self.plot==True:
			print PK[non0]
			plt.figure(figsize=(8,8))
			#print k
			#print PK
			plt.subplot('211')
			plt.scatter(k.take(non0), PK.take(non0))
			plt.loglog()
			plt.ylim(ymin=1.e-8)	
			plt.xlim(xmin=k.min(), xmax=k.max())
			plt.title('Power Spectrum')
			plt.xlabel('$k$')
			plt.ylabel('$P(k) (Kelvin^{2}(h^{-1}Mpc)^3)$')

			PK = PK*k*k*k/2./pi/pi
			#print PK
			plt.subplot('212')
			plt.scatter(k.take(non0), PK.take(non0))
			plt.loglog()
			plt.ylim(ymin=1.e-16)	
			plt.xlim(xmin=k.min(), xmax=k.max())
			plt.xlabel('$k (h Mpc^{-1})$')
			plt.ylabel('$\Delta^2 (Kelvin^{2})$')
			#plt.show()
			plt.savefig(out_root+'fkppower_'+resultf+'.eps', format='eps')

			PK2 = np.log10(PK2)
			plt.figure(figsize=(6,6))
			extent = (k2[0][0], k2[0][-1], k2[1][0], k2[1][-1])
			plt.imshow(PK2, origin='lower', extent = extent, interpolation='nearest')
			plt.xlabel('$k vertical (h Mpc^{-1})$')
			plt.ylabel('$k parallel (h Mpc^{-1})$')
			cb = plt.colorbar()
			cb.set_label('$lg(P^{2D}_{k_pk_v}) (Kelvin^2(h^{-1}Mpc)^3)$')
			plt.loglog()
			plt.savefig(out_root+'fkppower2_'+resultf+'.eps', format='eps')

			plt.show()
			#print 'Finished @_@ '
		return PK

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
		FKPPowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		FKPPowerSpectrumMaker().execute()
