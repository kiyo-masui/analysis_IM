"""Module that performs polarization and flux calibration."""

import scipy as sp
import numpy.ma as ma
import pylab as pl
import numpy as np

import kiyopy.custom_exceptions as ce
import base_single
import map.tools
from core import fits_map

# base_single.BaseSingle is a base class that knows how to read an input file
# loop over files, scans and IFs.  It does all the input and output.  The only
# thing it doesn't know how to do is the science, which is what will be added
# here.
class Calibrate(base_single.BaseSingle) :
    """Pipeline module that corrects for polarization leakage and flux calibrates."""
    # Here we define a bunch of stuff that BaseSingle needs to know to do its
    # thing.
    # prefix is a few letters that are added to all parameter names that are
    # read from a configuration file.
    prefix = 'tc_'
    # These are the parameters that should be read from file.  These are in
    # addition to the ones defined at the top of base_single.py.
    params_init = {'mueler_file' : 'default_fname' }
    
    # The base single initialization method does a bunch of stuff, but we want
    # to add one thing.  We want to read a mueler matrix from file.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
       	base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
#        print self.params
        # Read in the mueler matrix file.
#        i = self.file_ind
#        file_middle = self.params['file_middles'][i]
#        sess_num = file_middle.split('_')[0]
#        sess_num = int(sess_num)
#        print sess_num
#        mueler_file_name = self.params['mueler_file']+str(sess_num)+'_mueller_matrix_from_params.txt'
#        self.mueler
#        mueler_file_name = self.params['mueler_file']+'41_mueller_matrix_from_params.txt'
#        self.mueler = mueller(mueler_file_name)
    
    # This function tells BaseSingle what science to do.  Data is a
    # core.data_block.DataBlock object.  It holds all the data for a single
    # scan and a single IF.  BaseSingle knows how to loop over all of these.
    # More on DataBlock objects in the calibrate function below.
    def action(self, Data) :
#        Data.calc_freq()
#        frequency = Data.freq/1000000
#        pl.plot(frequency,Data.data[0,0,0,:],label='I-init')
#        pl.plot(frequency,Data.data[0,1,0,:],label='Q-init')
#        pl.plot(frequency,Data.data[0,2,0,:],label='U-init')
#        pl.plot(frequency,Data.data[0,3,0,:],label='V-init')

        # Main Action
        i = self.file_ind
        file_middle = self.params['file_middles'][i]
        sess_num = file_middle.split('_')[0]
        sess_num = int(sess_num)
        print sess_num
        mueler_file_name = self.params['mueler_file']+str(sess_num)+'_mueller_matrix_from_params.txt'
        self.mueler = mueller(mueler_file_name)
        calibrate_pol(Data, self.mueler)
       	Data.add_history('Flux calibrated and Corrected for polarization leakage.', 
               	         ('Mueller matrix file: ' + self.params['mueler_file'],))
        
#        pl.plot(frequency,Data.data[0,0,0,:],label='I-mod')
#        pl.plot(frequency,Data.data[0,1,0,:],label='Q-mod')
#        pl.plot(frequency,Data.data[0,2,0,:],label='U-mod')
#        pl.plot(frequency,Data.data[0,3,0,:],label='V-mod')
#        pl.legend()
#        pl.ylim(-20,130)
#        pl.xlabel("Frequency (MHz)")
#        pl.ylabel("Sample Data")
#        title0 = str(Data.field['SCAN'])+'_caloff_pol_params_'
#        pl.savefig(title0+'Comparison_Test_for_3C286.png')
#        pl.clf()

       	return Data

# m_total is the final Mueller matrix that needs to multiply by the stokes parameters.

# Note: I have it set up so that the third index represents the frequency 
# bin where 0 is the lowest frequency bin and 7 is the highest, and the 
# first and second indices represent the mueller matrix for each frquency bin. 

def mueller(mueler_file_name) :
    mp = np.loadtxt(mueler_file_name)
#     print mp
#This is a file with first index being freq, second index being matrix element:
# 0 = Freq, 1 = mII, 2 = mIQ, ... 16 = mVV

# Generates Mueller Matrix for use
    freq_limit = len(mp[:,0])
#     print freq_limit
    m_total = sp.zeros((4,4,freq_limit),float)
    for i in range(0,freq_limit):
        m_total[0,0,i] = mp[i,1]
        m_total[0,1,i] = mp[i,2]
        m_total[0,2,i] = mp[i,3]
        m_total[0,3,i] = mp[i,4]
        m_total[1,0,i] = mp[i,5]
        m_total[1,1,i] = mp[i,6]
        m_total[1,2,i] = mp[i,7]
        m_total[1,3,i] = mp[i,8]
        m_total[2,0,i] = mp[i,9]
        m_total[2,1,i] = mp[i,10]
        m_total[2,2,i] = mp[i,11]
        m_total[2,3,i] = mp[i,12]
        m_total[3,0,i] = mp[i,13]
        m_total[3,1,i] = mp[i,14]
        m_total[3,2,i] = mp[i,15]
        m_total[3,3,i] = mp[i,16]
        M_total = sp.mat(m_total[:,:,i])
#        M_total = M_total.I
#        print M_total
        for j in range(0,4):
            for k in range(0,4):
                m_total[j,k,i] = M_total[j,k]

    return m_total
def calibrate_pol(Data, m_total) :
    """Subtracts a Map out of Data."""
        
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

    # Also this depends on having the polarizations rotated correctly to IQUV.
    # Some code to do this has been hacked together in the rotate_pol module,
    # but I don't trust it yet.

    # Some dimension checks.
    # We expect 4 polarizations.
    if not Data.dims[1] == 4 :
       	raise ce.DataError('Require 4 polarizations.')
    # We expect polarizations to be in order IQUV.
    if (Data.field['CRVAL4'][0] != 1 or Data.field['CRVAL4'][1] != 2 or
        Data.field['CRVAL4'][2] != 3 or Data.field['CRVAL4'][3] != 4) :
       	raise ce.DataError('Expected the polarization basis to be IQUV.')

    # A useful function that might need:
    Data.calc_freq()
    # Now data has an atribute Data.freq which is an array that gives the
    # frequency along the last axis.
          
    # Data.field['CRVAL1'] is center frequency in Hz. 
    # Data.data 4 dim array 2nd index polarization, 4th index frequency. 

    # Need to get parallactic angle:
    Data.calc_PA()
    # This gives an array (Data.PA) of PA values of length = time dim.

    for time_index in range(0,Data.dims[0]):

    #Generate a sky matrix for this time index:
        m_sky = sp.zeros((4,4))
        m_sky[0,0] = 1
        m_sky[1,1] = ma.cos(2*Data.PA[time_index]*sp.pi/180)
        m_sky[1,2] = ma.sin(2*Data.PA[time_index]*sp.pi/180)
        m_sky[2,1] = -ma.sin(2*Data.PA[time_index]*sp.pi/180)
        m_sky[2,2] = ma.cos(2*Data.PA[time_index]*sp.pi/180)
        m_sky[3,3] = 1

        M_sky = sp.mat(m_sky)
        M_sky = M_sky.I
#        print M_sky

        for cal_index in range(0,Data.dims[2]):
        # Determines the Mueller Matrix to use   
            for freq in range(0,Data.dims[3]):

     # Tells which mueller matrix to use. 
               freq_limit = len(m_total[0,0,:])
               frequency = int(Data.freq[freq]/1000)
               bin = int((900000-frequency)*freq_limit/200000)
#               if freq_limit == 200:
#                   bin = 900-frequency
#Not setup to work with spectrometer data.
#               elif freq_limit == 260:
#                   bin = 929-frequency
#               print bin
    # Converts files into matrix format 
               STOKES = Data.data[time_index,:,cal_index,freq]       
#               print STOKES
               MUELLER = sp.mat(m_total[:,:,bin])
#               print MUELLER

    # Next there is a matrix multiplication that will generate 
    # a new set of stokes values.
               stokesmod = np.dot(MUELLER,STOKES)
               stokesmod = np.dot(M_sky,stokesmod)
#               print stokesmod
               for i in range(0,Data.dims[1]):
                    Data.data[time_index,i,cal_index,freq] = stokesmod[i]	

    # At this point the polarization values should be adjusted. 


# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
