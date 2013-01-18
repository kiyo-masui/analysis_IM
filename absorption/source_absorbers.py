"""
"""
import os
from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT
import utils.misc as utils

params_init = {
	       "input_root" : "./",
	       "file_middles" : ("GBTdata",),
               "input_end" : ".raw.acs.fits",
               "output_root" : "absorbers",
	       "scans" : (),
               "IFs" : (),
               }
prefix = 'abs_'

class SourceAbsorb(object):
    """ Stacks a set of source data to measure the absorption trough"""
  
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                          prefix=prefix)
        self.feedback = feedback

    def execute(self,nprocesses=1) :
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params,params['output_root']+'params.ini',
                               prefix=prefix)
        output_root = params['output_root']
        output_end = params['output_end']
        

        for file_middle in params['file_middles']:
            input_fname = (params['input_root']+file_middle+
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname)
	    n_IFs = len(Reader.IF_set)
            n_scans = len(Reader.scan_set)
            OnBlocks = Reader.read(0,n_scans,2),0,force_tuple=True)
            OffBlocks = Reader.read(1,n_scans,2),0,force_tuple=True)
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)

            
            m =0
            
