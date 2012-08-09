"""Plotting code for fits files in order to look at output data"""
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np
import pylab as pl

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT
from core import utils


# Define a dictionary with keys the names of parameters to be read from
# file and values the defaults.
params_init = {
               # Input and output.
               "input_root" : "./",
               "file_middles" : ("GBTdata",),
               "input_end" : ".raw.acs.fits",
               "output_root" : "calibrate",
               "output_end" : '',
               # Select data to process.
               "scans" : (),
               "IFs" : (),
               }
prefix = 'plot_'


class Plot_Data(object) :
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
# Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                          prefix=prefix)
        self.feedback = feedback

    def execute(self, nprocesses=1) :

        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        file_name = params['file_middles'][0].split('/')[1]

        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
# Read in the data, and loop over data blocks.
            Reader = core.fitsGBT.Reader(input_fname)
            n_IFs = len(Reader.IF_set) # Should be 1 given that we've stitched windows for the spectrometer or by def in guppi
            n_scans = len(Reader.scan_set) #Should be 4 for the spectrometer, 2 for guppi
#            print n_scans
            OnBlocks = Reader.read(range(0,n_scans,2),0,force_tuple=True)
            OffBlocks = Reader.read(range(1,n_scans,2),0,force_tuple=True)

#force_tuple=True makes the ouput of Reader.read a tuple even if thre is only one Block to return.
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)

# Setting labels for indices for later
            on_ind = 0
            off_ind = 1
            I_ind = 0
            XX_ind = 0
            YY_ind = 3
          
            for Data in Blocks:
                freq_len = Data.dims[3]
                Data.calc_freq()
                freq_val = Data.freq
                freq_val = freq_val/1000000

                S_med_calon_src = sp.zeros(freq_len)
                S_med_caloff_src = sp.zeros(freq_len)
                S_med_calon = sp.zeros(freq_len)
                S_med_caloff = sp.zeros(freq_len)

                for Data in OnBlocks:
                    S_med_caloff_src = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                    S_med_calon_src = ma.median(Data.data[:,I_ind,on_ind,:],axis=0)
#                    S_med_caloff_src = 0.5*(ma.median(Data.data[:,XX_ind,off_ind,:],axis=0)+ma.median(Data.data[:,YY_ind,off_ind,:],axis=0))
#                    S_med_calon_src = 0.5*(ma.median(Data.data[:,XX_ind,on_ind,:],axis=0)+ma.median(Data.data[:,YY_ind,off_ind,:],axis=0)) 
                    
                for Data in OffBlocks:
                    S_med_caloff = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                    S_med_calon = ma.median(Data.data[:,I_ind,on_ind,:],axis=0)
#                    S_med_caloff = 0.5*(ma.median(Data.data[:,XX_ind,off_ind,:],axis=0)+ma.median(Data.data[:,YY_ind,off_ind,:],axis=0))
#                    S_med_calon = 0.5*(ma.median(Data.data[:,XX_ind,on_ind,:],axis=0)+ma.median(Data.data[:,YY_ind,off_ind,:],axis=0))

                data = 0.5*(S_med_calon_src+S_med_caloff_src-S_med_calon-S_med_caloff)
            
            Isrc_3C286 = []
            Iact_3C286 = []
            freq_trunc = []
            for f in range(0,freq_len):
                if 838.00<freq_val[f]<841.00:
                    Isrc_3C286.append(19.74748409*pow((750.0/freq_val[f]),0.49899785)*(2.28315426-0.000484307905*freq_val[f]))
                    Iact_3C286.append(data[f])
                    freq_trunc.append(freq_val[f])

            pl.plot(freq_trunc,Isrc_3C286,label='Isrc')
            pl.scatter(freq_trunc,Iact_3C286,label='Iact')
            pl.legend()
            pl.xlabel("Frequency (MHz)")
            pl.ylabel("Temperature (K)")
            title0 = '3C286_absorber_test.png'
            pl.savefig(title0)
            pl.clf()

if __name__ == "__main__":
    import sys
    Plot_Data(str(sys.argv[1])).execute()

     
