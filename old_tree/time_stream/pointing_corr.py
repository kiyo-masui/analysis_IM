"""
Module that performs pointing corrections using a map and time stream data. 
Based upon the code written by Aravind. 
Original code found in cal/pointing.py 
Must be run in IQUV basis.
"""

import pylab as pl
import numpy as np
from numpy import *
from scipy import optimize

import kiyopy.custom_exceptions as ce
import base_single
import map.tools

from core import fits_map
from core import fitsGBT
import utils.misc as utils

#  base_single.BaseSingle is a base class that knows how to read an input file
# loop over files, scans and IFs.  It does all the input and output.  The only
# thing it doesn't know how to do is the science, which is what will be added
# here.
class Pointing(base_single.BaseSingle) :
    """Pipeline module that corrects for pointing errors."""
    # Here we define a bunch of stuff that BaseSingle needs to know to do its
    # thing.
    # prefix is a few letters that are added to all parameter names that are
    # read from a configuration file.
    prefix = 'pt_'
    # These are the parameters that should be read from file.  These are in
    # addition to the ones defined at the top of base_single.py.
    params_init = {
                   'map_dir' : 'default_fname', 
                   'map_file' : 'default_fname',
                    }
    
    # The base single initialization method does a bunch of stuff, but we want
    # to add one thing.  We want to read a mueler matrix from file.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
       	base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)

   
    # This function tells BaseSingle what science to do.  Data is a
    # core.data_block.DataBlock object.  It holds all the data for a single
    # scan and a single IF.  BaseSingle knows how to loop over all of these.
    # More on DataBlock objects in the calibrate function below.
    def action(self, Data) :
        # Main Action
        map_npy_file = self.params['map_dir']+self.params['map_file']+'.npy'

        pointing_corr(Data, map_npy_file)
       	Data.add_history('Corrected for Pointing Errors')
      	return Data

def read_map_data(map_npy_file):

   """  Reads the median map.
        Need map file as .npy and also a .npy.meta file.
   """
   in_file = open(map_npy_file+'.meta','r')
   info = in_file.readline()
   in_file.close()
   info1 = safe_eval(info)
   data_map = load(map_npy_file)
   dec_map = info1['dec_delta']*(arange(data_map.shape[2]) - data_map.shape[2]/2) + info1['dec_centre']
   freq_map = info1['freq_delta']*(arange(data_map.shape[0]) - data_map.shape[0]/2) + info1['freq_centre']
   ra_map = info1['ra_delta']*(arange(data_map.shape[1]) - data_map.shape[1]/2) + info1['ra_centre']

   return freq_map,dec_map,ra_map,data_map

def get_val(ra_val,dec_val,ra_map,dec_map,data,fnum,scan_point):

   """    Returns Temp from the map closest to (RA,DEC) in the scan.
   """
   if (ra_val > ra_map[0]) or (ra_val < ra_map[len(ra_map)-1]):
      return scan_point   # Outside the map dimensions, so return the scan value rather than the map value.
   #
   #   Alert the user here that a few scan points are too close to the map boundary 
   #   to compute a derivative. 
   #
   if (dec_val < dec_map[0]) or (dec_val > dec_map[len(dec_map)-1]):
      return scan_point   # Outside the map dimensions.

   for i in range(0,len(ra_map)):
      if ra_map[i] < ra_val:
         break
   for j in range(0,len(dec_map)):
      if dec_map[j] > dec_val:
         break

   # Interpolate in RA
   avg_minus = data[fnum][i-1][j-1] + ( ((data[fnum][i][j-1]-data[fnum][i-1][j-1])/(ra_map[i]-ra_map[i-1]))*(ra_val-ra_map[i-1]) )
   avg_plus  = data[fnum][i-1][j] + ( ((data[fnum][i][j]-data[fnum][i-1][j])/(ra_map[i]-ra_map[i-1]))*(ra_val-ra_map[i-1]) )
   # Interpolate in DEC
   avg = avg_minus + ( ((avg_plus-avg_minus)/(dec_map[j]-dec_map[j-1])) * (dec_val-dec_map[j-1]) )


   return avg

def get_pointing_error(fnum,scan_ra,scan_dec,map_ra,map_dec,map_data,scan_data):

   """ Computes the pointing errors, and returns the corrected scan.
       dRA and dDEC are arbitrary, but shouldn't be too small or too large.
       fnum can be 0,1,2,3 for 4 frequency bands 875,825,775,725 MHz.
   """
   dRA = 0.1
   dDEC = 0.1

   aver = mean(scan_data)
   scan_data = [val - aver for val in scan_data]
   temp_plus = [get_val(scan_ra[i]+dRA,scan_dec[i],map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp_plus = [temp_plus[i] - mean(temp_plus) for i in range(0,len(temp_plus)) if temp_plus[i] != -1.]
   temp = [get_val(scan_ra[i],scan_dec[i],map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp = [temp[i] - mean(temp) for i in range(0,len(temp))]
   prime_RA = [(temp_plus[i]-temp[i])/dRA for i in range(0,len(scan_ra))]

   temp_plus_temp = [get_val(scan_ra[i],scan_dec[i]+dDEC,map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp_plus = [val for val in temp_plus_temp if val != -1.]
   temp_plus = [temp_plus[i] - mean(temp_plus) for i in range(0,len(temp_plus))]
   prime_DEC = [(temp_plus[i]-temp[i])/dDEC for i in range(0,len(scan_ra))]

   deltaT = [scan_data[i]-temp[i] for i in range(0,len(scan_ra))]

   A1 = mean([deltaT[i]*prime_RA[i] for i in range(0,len(scan_ra))])
   A2 = mean([prime_RA[i]*prime_RA[i] for i in range(0,len(scan_ra))])
   A3 = mean([prime_RA[i]*prime_DEC[i] for i in range(0,len(scan_ra))])
   B1 = mean([deltaT[i]*prime_DEC[i] for i in range(0,len(scan_ra))])
   B2 = A3
   B3 = mean([prime_DEC[i]*prime_DEC[i] for i in range(0,len(scan_ra))])

   DeltaRA = (A1*B3 - B1*A3) / (A2*B3 - B2*A3)
   DeltaDEC = (A1*B2 - B1*A2) / (A3*B2 - B3*A2)
   new_scan = [aver + scan_data[i] - (prime_RA[i]*DeltaRA + prime_DEC[i]*DeltaDEC) for i in range(0,len(scan_ra))]
   return DeltaRA, DeltaDEC, new_scan

def pointing_corr(Data, map_npy_file) :
    """Corrects the RA and DEC using a map."""
        
    # Data is a DataBlock object.  It holds everything you need to know about
    # the data in a single scan and IF.  You should get to know them very well.
    # Data.data is a numpy masked array (see numpy documentation) and holds the
    # acctual data.  It is a 4 dimensional array.  The demensions are (in
    # order): (time, pol, cal, freq).  Each dimension can be any length which
    # you can figure out by looking at Data.dims = sp.shape(Data.data).
    # Data.field is a python dictionary that holds all the other data that you
    # might care about from the origional fits file.  For instance,
    # Data.field['CAL'] is an array with length dims[2].  It normally has
    # values ['T', 'F']. Data.field['CRVAL4'] tells you about the polarization
    # axis of Data.data.  By SDfits convension each polarization is represented
    # by an integer: 1=I, 2=Q, 3=U, 4=V, -5=XX, -6=YY, -7=XY, -8=YX.
    freq_map,dec_map,ra_map,data_map = read_map_data(map_npy_file)

    Data.calc_pointing()
    Data.calc_freq()

    for fnum in range(0,len(Data.freq)):
        scan_val = [Data.data[i,0,0,fnum] for i in range(0,len(Data.ra))]
        deltaRA,deltaDEC,new_scan = get_pointing_error(fnum,Data.ra,Data.dec,ra_map,dec_map,data_map,scan_val)
        for i in range(0,len(Data.ra)):
            Data.data[i,0,0,fnum] = new_scan[i]

# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Pointing(str(sys.argv[1])).execute()
