#! /usr/bin/env python
'''
    This module is used to calculate the bias(Transfer Function)
    caused by the mode subtraction. 

    The bias is defined as the ratio of the power spectrum 
    between the simulation maps with and without mode subtraction.

'''

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
from scipy.optimize import leastsq

import functions

params_init = {
    'processes' : 1,
    'plot' : True,
    'input_root' : '../newmaps/',
    'resultf' : '',
    'resultf0' : '',
    'output_root' : './',

    'simmap_root' : './',

    'kbinNum' : 20,
    'kmin' : None,
    'kmax' : None,

    'FKPweight' : False,

    'OmegaHI' : 1.e-3,
    'Omegam' : 0.23,
    'OmegaL' : 0.74,
    'z' : 1,
    'PKunit' : 'K',
    'cross' : False,
}

pi = 3.1415926
deg2rad = pi/180.
prefix = 'bc_'

class BiasCalibrate(object):
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
        in_root = params['input_root']
        resultf = functions.getresultf(params)
        resultf0= params['resultf0']

        FKPweight = params['FKPweight']

        # Read the Power Spectrum Result
        # clean_{map+sim}(map+sim)
        k1 = sp.load(in_root + resultf + '_k_combined.npy')
        if FKPweight:
            pk1 = sp.load(in_root + resultf + '_p_fkp_combined.npy')
            pk12= sp.load(in_root + resultf + '_p2_fkp_combined.npy')
        else:
            pk1 = sp.load(in_root + resultf + '_p_combined.npy')
            pk12= sp.load(in_root + resultf + '_p2_combined.npy')

        #non0 = pk1.nonzero()
        #pk1 = pk1.take(non0)[0]
        #k1 = k1.take(non0)[0]
        # clean_{map}(map)
        k0 = sp.load(in_root + resultf0 + '_k.npy')
        if FKPweight:
            pk0 = sp.load(in_root + resultf0 + '_p_fkp.npy')
            pk02= sp.load(in_root + resultf0 + '_p2_fkp.npy')
        else:
            pk0 = sp.load(in_root + resultf0 + '_p.npy')
            pk02= sp.load(in_root + resultf0 + '_p2.npy')

        #pk0 = pk0.take(non0)[0]
        #k0 = k0.take(non0)[0]

        # Read the Pwer Spectrum without mode subtraction
        simmap_root = params['simmap_root']
        k  = sp.load(simmap_root + 'simmaps_k.npy')
        pk = sp.load(simmap_root + 'simmaps_p.npy')
        pk2= sp.load(simmap_root + 'simmaps_p2.npy')

        # Test if k and k0 are match
        #print simmap_root 
        #print 
        #print k0
        #print
        #print k
        #print
        print k.shape
        print k0.shape
        print k1.shape
        if (k-k0).any() or (k-k1).any():
            print "k and k0 are not match!!"
            return 0

        dpk = pk1-pk0
        dpk[dpk<=0] = np.inf
        dpk = pk/dpk

        dpk2= (pk12-pk02)
        dpk2[dpk2<=0] = np.inf
        dpk2= pk2/dpk2

        B = interp1d(k, dpk, kind='cubic')

        ki = np.logspace(log10(k.min()+0.001), log10(k.max()-0.001), num=300)
        bias = B(ki)

        if params['cross']:
            dpk = np.sqrt(dpk)
            dpk2 = np.sqrt(dpk2)

        sp.save(out_root + 'k_bias', k)
        sp.save(out_root + 'b_bias', dpk)
        sp.save(out_root + 'b2_bias', dpk2)


        if self.plot==True:
            plt.figure(figsize=(8,4))
            plt.subplot('111')
            plt.scatter(k, dpk, s=30, c='w', marker='s')
            plt.plot(ki, bias)
            #plt.scatter(k, dpk)
            plt.loglog()
            #plt.semilogx()
            plt.ylim(ymin=dpk.min())
            plt.xlim(xmin=k.min(), xmax=k.max())
            plt.title('Power Spectrum Bias')
            plt.xlabel('$k$')
            plt.ylabel('$dP(k) (mk^{2}(h^{-1}Mpc)^3)$')

            plt.savefig(out_root+'b_bias.eps', format='eps')
            #plt.show()

#   def getbias(self, pke, begin, end, pk0):
#       print pke[begin:end]
#       pkA = pke[begin:end].mean(axis=0)
#       pke[begin:end] = (pke[begin:end][:]-pkA)**2
#       pkvarA = np.sum(pke[begin:end], axis=0)
#       pkvarA = pkvarA/(end-begin)
#       pkvarA = np.sqrt(pkvarA)
#       print pk0
#       print pkA
#       print pkvarA
#       pkerrA = np.ndarray(shape=(2, len(pkA)))
#       pkerrA[0] = pkvarA
#       pkerrA[1] = pkvarA
#       dpkA = pk0/pkA
#       print dpkA
#       return pkA, pkerrA, dpkA

if __name__ == '__main__':
    import sys
    if len(sys.argv)==2 :
        BiasCalibrate(str(sys.argv[1])).execute()
    elif len(sys.argv)>2 :
        print 'Maximun one argument, a parameter file name.'
    else :
        BiasCalibrate().execute()
