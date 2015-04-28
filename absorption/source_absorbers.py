"""
"""
import os
from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np
import pylab 

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT
import utils.misc as utils

params_init = {
	       "input_root" : "./",
	       "file_middles" : ("GBTdata",),
               "input_end" : ".raw.acs.fits",
               "cal" : "N",
               "plot_label" : "",
               "output_root" : "absorbers",
               "output_end" : ".txt",
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
       
        self.file_num = len(params['file_middles'])
        print self.file_num

        sample_file = params['file_middles'][0]
        input_sample = (params['input_root']+sample_file+params['input_end'])
        Reader = core.fitsGBT.Reader(input_sample)
        off_sample = Reader.read(1,0,force_tuple=True)
        for Data in off_sample:
            freq_len = Data.dims[3]        
        
        On_tot = []
        Off_tot = [] 
        RAs_on = []
        DECs_on = []
        RAs_off = []
        DECs_off = []
                        
        for file_middle in params['file_middles']:
            input_fname = (params['input_root']+file_middle+
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname)
	    n_IFs = len(Reader.IF_set)
            n_scans = len(Reader.scan_set)
            if n_scans == 2:
	        OnBlocks = Reader.read(range(0,n_scans,2),0,force_tuple=True)
        	OffBlocks = Reader.read(range(1,n_scans,2),0,force_tuple=True)
	        Blocks = Reader.read(params['scans'], params['IFs'],
        	                     force_tuple=True)
            for Data in OnBlocks:
                freq_len = Data.dims[3]
                time_len = Data.dims[0] 
                Data.calc_freq() 
                freq_val = Data.freq
                freq_val = freq_val/1000000 
                Full_date = Data.field['DATE-OBS']
                print Full_date[0],Full_date[-1],freq_val[1]-freq_val[0]
#                p#rint Full_date
                Az = Data.field['CRVAL2'] 
#                print Az 
                El = Data.field['CRVAL3'] 
#                print El 
#                print len(Az),len(El),len(Full_date)
                RA_data = sp.zeros(len(Az)) 
                DEC_data = sp.zeros(len(Az))
                for i in range(0,len(Az)): 
                    RA_data[i], DEC_data[i] = utils.azel2radecGBT(Az[i],El[i],Full_date[i])
#                    print RA_data[i], DEC_data[i]
#                print ma.mean(RA_data)
                off_data = Data.data[:,0,1,:]
                off_med = ma.median(off_data,axis=0)# This is cal off, not src off
                On_tot.append(off_med)
                RAs_on.append(ma.mean(RA_data))
                DECs_on.append(ma.mean(DEC_data))
 
	    for Data in OffBlocks:
	        freq_len = Data.dims[3]
                time_len = Data.dims[0]
                Data.calc_freq()
                freq_val = Data.freq
                freq_val = freq_val/1000000
                Full_date = Data.field['DATE-OBS']
#                p#rint Full_date
                Az = Data.field['CRVAL2']
#                print Az
                El = Data.field['CRVAL3']
#                print El
#                print len(Az),len(El),len(Full_date)
                RA_data = sp.zeros(len(Az))
                DEC_data = sp.zeros(len(Az))
                for i in range(0,len(Az)):
                    RA_data[i], DEC_data[i] = utils.azel2radecGBT(Az[i],El[i],Full_date[i])
#                    print RA_data[i], DEC_data[i]
                
                off_data = Data.data[:,0,1,:]
                off_med = ma.median(off_data,axis=0)# This is cal off, not src off
                Off_tot.append(off_med)
                RAs_off.append(ma.mean(RA_data))
                DECs_off.append(ma.mean(DEC_data))
                
#        print np.shape(Data_tot)        
#        print RAs
#        print DECs
#        Data_tot = np.array(Data_tot)

        offlarge = []
        offsmall = []
        for i in range(0,len(Off_tot)):
            if RAs_off[i]>203.0:
                offlarge.append(Off_tot[i])
            else:
                offsmall.append(Off_tot[i])
        
        On_med = ma.median(On_tot,axis=0)
        Off_med_lrg = ma.median(offlarge,axis=0)
        Off_med_sml = ma.median(offsmall,axis=0)
        
        min =0
        max = len(freq_val)
        for i in range(0,len(freq_val)):
            if freq_val[i]>=835.0:
                min = i
            if freq_val[i]>=845.0:
                max = i
        
        plot_label = params['plot_label']
        calibrated = params['cal']
#        print freq_val[min],freq_val[max]
        pylab.scatter(freq_val[max:min],On_med[max:min],c='b',edgecolor='b',s=3,label='On Src')
        pylab.vlines(839.4,0.,80.,color='r')
        if calibrated == 'N':
            pylab.ylim(27,32)
            pylab.ylabel('Temperature (Noise Cal Units)')
        elif calibrated == 'Y':
            pylab.ylim(53,57)
            pylab.ylabel('Temperature (K)')
        pylab.xlabel('Frequency (MHz)')
        pylab.savefig('absorb_stack_3C286_'+plot_label)
        pylab.clf()

        pylab.scatter(freq_val[max:min],Off_med_lrg[max:min],c='g',edgecolor='g',s=3,label='Off Above')
        pylab.scatter(freq_val[max:min],Off_med_sml[max:min],c='r',edgecolor='r',s=3,label='Off Below')
        if calibrated=='N':
           pylab.ylim(8,12)
           pylab.ylabel('Temperature (Noise Cal Units)')
        elif calibrated=='Y':
            pylab.ylim(15,25)
            pylab.ylabel('Temperature (K)')
        pylab.xlabel('Frequency (MHz)')
        pylab.savefig('absorb_stack_off_src_'+plot_label)
        pylab.clf()

       
#        np.savetxt('3C286_absorb_data.txt',stack_data,delimiter=' ')


         
if __name__ == "__main__":
    import sys
    SourceAbsorb(str(sys.argv[1])).execute()  
