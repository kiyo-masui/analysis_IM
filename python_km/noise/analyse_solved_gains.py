"""Analyse the gains solved for and stored by time_stream/subtract_map.py.
"""

import cPickle
import time

import scipy as sp
import numpy.ma as ma

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce

# TODO: See if gain correlates with LST, checks if map might not have constant
# gain.

# Parameters prefixed with 'sg_' when read from file.
params_init = {
              # IO:
              'input_root' : './',
              # The unique part of every fname
              'file_middles' : ("testfile_GBTfits",),
              'input_end' : "_gain.pickle"
              }
prefix = 'sg_'


class SolvedGains(object) :
    
    def __init__(self, parameter_file_or_dict=None) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix)

    def execute(self) :
        
        params = self.params
        n_files = len(params['file_middles'])
        scan_number = 0

        for file_middle in params['file_middles'] :
            file_name = (params['input_root'] + file_middle + 
                         params['input_end'])
            f = open(file_name, 'r')
            scan_list = cPickle.load(f)
            n_scans_file = len(scan_list)
            for scan in scan_list :
                if scan_number == 0 :
                    n_scans = n_scans_file*n_files
                    gain = ma.empty((n_scans,) + scan['gain'].shape)
                    time_arr = sp.empty(n_scans, dtype=int)
                gain[scan_number, ...] = scan['gain']
                to = time.strptime(str(scan['time']), "%Y_%m_%d_%H:%M:%S")
                time_arr[scan_number] = to.tm_sec + 60*(to.tm_min + 
                    60*(to.tm_hour + 24*(to.tm_yday + 365*(to.tm_year-2000))))
                scan_number += 1
        self.time = time_arr[:scan_number]
        self.gain = gain[:scan_number,...]


        



