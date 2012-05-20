#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
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
	'input_root' : '../newmaps/',
	'hr' : ('15hr_40-41-43_','15hr_42_',),
	'mid' : 'dirty_map_'
	'polarizations' : ('I',),
	'output_root' : './',

	'boxshape' : (128,128,128),
	'boxunit' : 15., # in unit Mpc
	'discrete' : 3,
	'Xrange' : (1400,),
	'Yrange' : (-64*15,64*15),
	'Zrange' : (0.,128*15),
}
prefix = 'pk_'

class PowerSpectrumMaker(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params

		# Make parent directory and write parameter file.
		kiyopy.utils.mkparents(params['output_root'])
		parse_ini.write_params(params, params['output_root']+'params.ini',prefix='pk_')
		in_root = params['input_root']
		out_root = params['output_root']
		mid = params['mid']
		all_out_fname_list = []
		all_in_fname_list = []
		
		#### Process ####
		pol_str = params['polarizations'][0]
		hr_str = params['hr'][0]
		imap_fname = in_root + hr_str + mid + pol_str + '.npy'
		imap = algebra.load(imap_fname)
		imap = algebra.make_vect(imap)
		if imap.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		nmap_fname = in_root + hr_str + 'noise_inv_diag_' + pol_str + '.npy'
		nmap = algebra.load(nmap_fname)
		nmap = algebra.make_vect(nmap)

		#noise normal
		normal = (nmap**2).sum()

		#Using map in different day 
		hr_str = params['hr'][1]
		imap_fname2 = in_root + hr_str + mid + pol_str + '.npy'
		imap2 = algebra.load(imap_fname2)
		imap2 = algebra.make_vect(imap2)
		if imap2.axes != ('freq', 'ra', 'dec') :
			raise ce.DataError('AXES ERROR!')

		nmap_fname = in_root + hr_str + 'noise_inv_diag_' + pol_str + '.npy'
		nmap2 = algebra.load(nmap_fname)
		nmap2 = algebra.make_vect(nmap2)

		#noise normal
		normal2 = (nmap2**2).sum()

		normal = sqrt(normal)*sqrt(normal2)

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
		#imap = imap.swapaxes(1,2)  # change the ra and dec
		#print imap.shape

		mapinf = [dr, dra, dde, disc_n]
		mapinf = np.array(mapinf)
		#print mapinf

		#print r	
		#print ra
		#print de

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
		V = params['boxunit']**3

		boxinf = [xrange0, yrange0, zrange0, boxunit]
		boxinf = np.array(boxinf)

		print "Filling the BOX"
		MakePower.Filling(imap, imap2, box, box2, r, ra, de, boxinf, mapinf)

		print "FFTing "
		fftbox = fftn(box)
		fftbox = fftbox.real**2 + fftbox.imag**2
		fftbox2= fftn(box2)
		fftbox2 = fftbox2.real**2 + fftbox2.imag**2
		fftbox = fftbox2
		#fftbox = fftbox.__pow__(0.5)*fftbox2.__pow__(0.5)

		PK = np.zeros(40)
		k = np.zeros(40)
		PK2 = np.zeros(shape=(10, 10))
		k2 = np.zeros(shape=(2, 10))
		MakePower.Make(fftbox, PK, k, PK2, k2)
		kunit = 2.*pi/(boxshape[0]*boxunit)
		k = k*kunit
		k2 = k2*kunit
		PK = PK*V*params['boxshape'][0]**3/normal
		PK2 = PK2*V*params['boxshape'][0]**3/normal

		sp.save(out_root+'PK', PK)
		sp.save(out_root+'PK2', PK2)

		non0 = PK.nonzero()

		if self.plot==True:
			plt.figure(figsize=(8,8))
			#print k
			#print PK
			plt.subplot('211')
			plt.scatter(k.take(non0), PK.take(non0))
			plt.loglog()
			plt.ylim(ymin=1.e1)	
			plt.xlim(xmin=k.min())
			plt.title('Power Spectrum')
			plt.xlabel('$k$')
			plt.ylabel('$P(k) (Kelvin^{2}(h^{-1}Mpc)^3)$')

			PK = PK*V*params['boxshape'][0]**3/1.2e12*k*k*k/2./pi/pi
			#print PK
			plt.subplot('212')
			plt.scatter(k.take(non0), PK.take(non0))
			plt.loglog()
			plt.ylim(ymin=1.e-9)	
			plt.xlim(xmin=k.min())
			plt.xlabel('$k (h Mpc^{-1})$')
			plt.ylabel('$\Delta^2 (Kelvin^{2})$')
			#plt.show()
			plt.savefig(out_root+'power.eps', format='eps')

			PK2 = np.log10(PK2)
			plt.figure(figsize=(6,6))
			extent = (k2[0][0], k2[0][-1], k2[1][0], k2[1][-1])
			plt.imshow(PK2, origin='lower', extent = extent, interpolation='nearest')
			plt.xlabel('$k vertical (h Mpc^{-1})$')
			plt.ylabel('$k parallel (h Mpc^{-1})$')
			cb = plt.colorbar()
			cb.set_label('$lg(P^{2D}_{k_pk_v}) (Kelvin^2(h^{-1}Mpc)^3)$')
			plt.loglog()
			plt.savefig(out_root+'power2.eps', format='eps')

			#plt.show()
			print 'Finished @_@ '
		return PK

#
#		
#		#print r.min(), r.max(), r.max()-r.min()
#		#print ra.min(), ra.max(), ra.max()-ra.min()
#		#print de.min(), de.max(), de.max()-de.min()


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

	def fgrm(self, imap, fname):
		freq = imap.get_axis('freq')/1.e6
		imap = imap.swapaxes(0,2)
		fg = algebra.info_array(sp.zeros(imap.shape))
		fg.axes = ('ra','dec','freq')
		fg = algebra.make_vect(fg)

		if self.plot==True:
			plt.figure(figsize=(8,12))
			#plot the map without foreground removeing
			plt.subplot(311)
			#plt.title('Map without foreground remove')
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) With Foreground')
			for i in range(0,imap.shape[0]):
				#for j in range(1, 2):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap[i][j])

		for i in range(0,imap.shape[0]):
			#for j in range(1, 2):
			for j in range(0, imap.shape[1]):
				y = imap[i][j]
				y = np.asmatrix(y).T
				X = np.ones(shape=(3,imap.shape[2]))
				X[1] = np.log10(freq)
				X[2] = np.square(np.log10(freq))
				XT = np.asmatrix(X)
				X = XT.T
				a = ((XT*X).I*XT)*y
				yy = np.asarray((X*a).T)
				y = np.asarray((y-X*a).T)
				imap[i][j] = y
				fg[i][j] = yy

		if self.plot==True:
			plt.subplot(312)
			#plt.title('Foreground')
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) Foreground')
			for i in range(0,imap.shape[0]):
				#for j in range(1, 2):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, fg[i][j])
	
			plt.subplot(313)
			#plt.title('Map with foreground remove')
			plt.xlabel('Frequece (MHz)')
			plt.ylabel('$\Delta$ T(Kelvin) Without Foreground')
			for i in range(0,imap.shape[0]):
				#for j in range(1, 2):
				for j in range(0, imap.shape[1]):
					plt.plot(freq, imap[i][j])
			plt.savefig(fname, format='png')
			#plt.show()


if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		PowerSpectrumMaker(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		PowerSpectrumMaker().execute()
