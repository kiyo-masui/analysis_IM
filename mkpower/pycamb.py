#! /usr/bin/env python

import scipy as sp
import numpy as np
import scipy.linalg as linalg
import os
from math import *
from sys import *
import matplotlib.pyplot as plt

import pycamb

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce

pi = 3.1415926

params_init = {
	'processes' : 1,
	'plot' : 'T',
	#'hr' : ('15hr_40-41-43_','15hr_42_',),
	#'mid' : 'dirty_map_',
	#'polarizations' : ('I',),
	'output_root' : './cambini/',

	'get_scalar_cls' : 'F',
	'get_vector_cls' : 'F',
	'get_tensor_cls' : 'F',
	'get_transfer'   : 'T',

	'do_lensing'     : 'F',

	#nonliner: 0, linear;  1 nonlinear matter power spectrum; 2, CMB lensing
	'do_nonlinear' : 1,

	'l_max_scalar'      : 2000,
	'k_eta_max_scalar'  : 4000,
	'l_max_tensor'      : 1500,
	'k_eta_max_tensor'  : 3000,
	'use_physical'   : 'T',
	'ombh2'          : 0.0226,
	'omch2'          : 0.112,
	'omnuh2'         : 0,
	'omk'            : 0,
	'hubble'         : 70.2,
	'w'              : -1,
	'cs2_lam'        : 1,
	'temp_cmb'           : 2.725,
	'helium_fraction'    : 0.24,
	'massless_neutrinos' : 3.04,
	'massive_neutrinos'  : 0,
	'nu_mass_eigenstates' : 1,
	'nu_mass_degeneracies' : 0,  
	'nu_mass_fractions' : 1,
	'initial_power_num'         : 1,
	'pivot_scalar'              : 0.05,
	'pivot_tensor'              : 0.05,
	'scalar_amp(1)'             : 2.1e-9,
	'scalar_spectral_index(1)'  : 0.96,
	'scalar_nrun(1)'            : 0,
	'tensor_spectral_index(1)'  : 0,
	'initial_ratio(1)'          : 1,
	'reionization'         : 'T',
	're_use_optical_depth' : 'T',
	're_optical_depth'     : 0.09,
	're_redshift'          : 11,
	're_delta_redshift'   : 1.5,
	're_ionization_frac'   : -1,
	'RECFAST_fudge' : 1.14,
	'RECFAST_fudge_He' : 0.86,
	'RECFAST_Heswitch' : 6,
	'RECFAST_Hswitch'  : 'T',
	'initial_condition'   : 1,
	'initial_vector' : '-1 0 0 0 0',
	'vector_mode' : 0,
	'COBE_normalize' : 'F',
	'CMB_outputscale' : 7.4311e12,
	'transfer_high_precision' : 'F',
	'transfer_kmax'           : 2,
	'transfer_k_per_logint'   : 0,
	'transfer_num_redshifts'  : 1,
	'transfer_interp_matterpower' : 'T',
	'transfer_redshift'		  : 0,
	'transfer_redshift(1)'    : 0,
	'transfer_filename(1)'    : 'Transfer_out.dat',
	'transfer_matterpower(1)' : 'matterpower.dat',
	'scalar_output_file' : 'scalCls.dat',
	'vector_output_file' : 'vecCls.dat',
	'tensor_output_file' : 'tensCls.dat',
	'total_output_file'  : 'totCls.dat',
	'lensed_output_file' : 'lensedCls.dat',
	'lensed_total_output_file'  :'lensedtotCls.dat',
	'lens_potential_output_file' : 'lenspotentialCls.dat',
	'FITS_filename'      : 'scalCls.fits',
	'do_lensing_bispectrum' : 'F',
	'do_primordial_bispectrum' : 'F',
	'bispectrum_nfields' : 2,
	'bispectrum_slice_base_L' : 0,
	'bispectrum_ndelta':3,
	'bispectrum_delta(1)':0,
	'bispectrum_delta(2)':2,
	'bispectrum_delta(3)':4,
	'bispectrum_do_fisher': 'F',
	'bispectrum_fisher_noise':0,
	'bispectrum_fisher_noise_pol':0,
	'bispectrum_fisher_fwhm_arcmin':7,
	'bispectrum_full_output_file': '',
	'feedback_level' : 1,
	'lensing_method' : 1,
	'accurate_BB': 'F',
	'massive_nu_approx' : 3,
	'accurate_polarization'   : 'T',
	'accurate_reionization'   : 'T',
	'do_tensor_neutrinos'     : 'F',
	'do_late_rad_truncation'   : 'T',
	'number_of_threads'       : 0,
	'accuracy_boost'          : 1,
	'l_accuracy_boost'        : 1,
	'l_sample_boost'          : 1


}
prefix = 'pcb_'

class CAMB(object):
	"""Calculate Power Spectrum"""

	def __init__(self, parameter_file_or_dict=None, feedback=2):
		# Read in the parameters.
		self.params = parse_ini.parse(parameter_file_or_dict, params_init, prefix=prefix, feedback=feedback)

		self.feedback=feedback

		self.plot = bool(self.params['plot'])
	
	def execute(self, nprocesses=1):
		params = self.params
		params['transfer_redshift(1)'] = params['transfer_redshift']


		inifilename= params['output_root'] + 'params.ini'
		inifile = open(inifilename, 'w')
		#parse_ini.write_params(params, inifile ,prefix='')
		try:
			for key in params:
				print  >>inifile, '{0} = {1}'.format(key, params[key])
		finally: inifile.close()

		cambpath = os.getenv('CAMBPATH') + 'camb'

		os.system(cambpath + ' ' + inifilename)

		P = []
		k = []
		fname = params['output_root'] + '_matterpower.dat'
		f = open(fname, 'r')
		data = f.readlines()
		for line in data:
			line = line.split()
			k.append(line[0])
			P.append(line[1])

		f.close()

		PK = np.ndarray(shape=(2,len(k)))
		PK[0] = k
		PK[1] = P
		sp.save(params['output_root']+'PKcamb', PK)


		#plt.figure(figsize=(8,4))
		#plt.subplot('111')
		#plt.plot(k, P, c='r')
		#plt.plot(PK[0], PK[1], c='r')
		##plt.scatter(k, P, c='r')
		#plt.loglog()
		#plt.ylim(ymin=1.e1)
		#plt.xlim(xmin=1.e-2, xmax=1.)
		#plt.show()

	#def funcdl(self, z, omegam, omegal):



if __name__ == '__main__':
	import sys
	if len(sys.argv)==2 :
		CAMB(str(sys.argv[1])).execute()
	elif len(sys.argv)>2 :
		print 'Maximun one argument, a parameter file name.'
	else :
		CAMB().execute()
