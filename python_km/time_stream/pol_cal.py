"""Module that performs polarization calibration."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single
import map.tools
from core import fits_map

# base_single.BaseSingle is a base class that knows how to read an input file
# loop over files, scans and IFs.  It does all the input and output.  The only
# thing it doesn't know how to do is the science, which is what will be added
# here.
class Calibrate(base_single.BaseSingle) :
    """Pipeline module that corrects for polarization leakage."""
    
    # Here we define a bunch of stuff that BaseSingle needs to know to do its
    # thing.
    # prefix is a few letters that are added to all parameter names that are
    # read from a configureation file.
    prefix = 'pc_'
    # These are the parameters that should be read from file.  These are in
    # addition to the ones defined at the top of base_single.py.
    params_init = {'mueler_file' : 'default_fname' }

    # The base single initialization method does a bunch of stuff, but we want
    # to add one thing.  We want to read a mueler matrix from file.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
       	base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
        # Read in the mueler matrix file.
       	mueler_file_name = self.params['mueler_file']
       	self.mueler = mueller()
    
    # This function tells BaseSingle what science to do.  Data is a
    # core.data_block.DataBlock object.  It holds all the data for a single
    # scan and a single IF.  BaseSingle knows how to loop over all of these.
    # More on DataBlock objects in the calibrate function below.
    def action(self, Data) :
       	calibrate_pol(Data, self.mueler)
       	Data.add_history('Corrected for polarization leakage.', 
               	         ('Mueler matrix file: ' + self.params['mueler_file'],))

       	return Data

# m_total is the final Mueller matrix that needs to be multiplied by the stokes parameters to adjust them.

#Note: I have it set up so that the third index represents the frequency bin where 0 is the lowest frequency bin and 7 is the highest, and the first and second indices represent the mueller matrix for each frquency bin. 

def mueller() :
# Mueller Matrix parameters from polcal calculations
    deltaG = [-0.185,-0.113,-0.213,-0.295,-0.394,-0.012,0.546,0.680]
    psi = [-21.9,179.9,-175.9,178.9,171.2,152.5,159.9,122.5]
    alpha = [-1.2,3.5,2.9,-1.3,10.0,-6.2,-2.2,-7.5]
    epsilon = [0.008,0.013,0.015,0.019,0.024,0.016,0.016,0.017]
    phi = [-164.9,-4.0,-7.9,0.5,6.1,23.2,17.6,55.4]
    CFR = [694,724,754,784,814,844,874,904]
    delf = sp.array([(CFR[i]-1485.8) for i in range(0,8)])
    psi = sp.array([(psi[i]-0.0086*delf[i]) for i in range(0,8)])

    mp = sp.array([deltaG,psi,alpha,epsilon,phi,CFR])

# Amplifier Mueller Matrix
    m_a = sp.zeros((4,4,8), float)
    for i in range(0,8):
        m_a[0,0,i] = 1
        m_a[1,0,i] = 0.5*mp[0,i]
        m_a[0,1,i] = 0.5*mp[0,i]
        m_a[1,1,i] = 1
        m_a[2,2,i] = ma.cos(mp[1,i]*sp.pi/180)
        m_a[3,2,i] = -ma.sin(mp[1,i]*sp.pi/180)
        m_a[2,3,i] = -m_a[3,2,i]
        m_a[3,3,i] = m_a[2,2,i]

# Feed Mueller Matrix
    m_f = sp.zeros((4,4,8), float)
    for i in range(0,8):
        m_f[0,0,i] = 1
        m_f[1,1,i] = ma.cos(mp[2,i]*sp.pi/180)*ma.cos(mp[2,i]*sp.pi/180)-ma.sin(mp[2,i]*sp.pi/180)*ma.sin(mp[2,i]*sp.pi/180)
        m_f[3,1,i] = 2*ma.cos(mp[2,i]*sp.pi/180)*ma.sin(mp[2,i]*sp.pi/180)
        m_f[2,2,i] = ma.cos(mp[2,i]*sp.pi/180)*ma.cos(mp[2,i]*sp.pi/180)+ma.sin(mp[2,i]*sp.pi/180)*ma.sin(mp[2,i]*sp.pi/180)
        m_f[1,3,i] = -2*ma.cos(mp[2,i]*sp.pi/180)*ma.sin(mp[2,i]*sp.pi/180)
        m_f[3,3,i] = ma.cos(mp[2,i]*sp.pi/180)*ma.cos(mp[2,i]*sp.pi/180)-ma.sin(mp[2,i]*sp.pi/180)

# Feed Imperfections Mueller Matrix
    m_ilfr = sp.zeros((4,4,8), float)
    for i in range(0,8):
        m_ilfr[0,0,i] = 1
        m_ilfr[2,0,i] = 2*mp[3,i]*ma.cos(mp[4,i]*sp.pi/180)
        m_ilfr[3,0,i] = 2*mp[3,i]*ma.sin(mp[4,i]*sp.pi/180)
        m_ilfr[1,1,i] = 1
        m_ilfr[0,2,i] = m_ilfr[2,0,i]
        m_ilfr[2,2,i] = 1
        m_ilfr[0,3,i] = m_ilfr[3,0,i]
        m_ilfr[3,3,i] = 1

# Rotates Matrix into Astronomical frame
    m_astron = sp.zeros((4,4,8), float)
    for i in range(0,8):
        m_astron[0,0,i] = 1
        m_astron[3,3,i] = 1
        m_astron[1,1,i] = -1
        m_astron[2,2,i] = -1

# Generates Inverse Mueller Matrix for use
    m_total = sp.zeros((4,4,8), float)
    for i in range(0,8):
        M_a = sp.mat(m_a[:,:,i])
        M_f = sp.mat(m_f[:,:,i])
        M_ilfr = sp.mat(m_ilfr[:,:,i])
        M_tot = M_a*M_f*M_ilfr
        M_total = M_tot.I
        M_astron = sp.mat(m_astron[:,:,i])
        M_total = M_astron*M_total
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
    # These arrays should be of the form stokes[i,4] where i is the number 
    # of data points in the particular frequency bin. Data.field['CRVAL1'] 
    # is center frequency in Hz. 
    # Data.data 4 dim array 2nd column polarization, 4th column frequency. 
   
    for time_index in range(0,Data.dims[0]):
        for cal_index in range(0,Data.dims[2]):

            CenterFrequency = int(Data.field['CRVAL1']/1000000)
            
            if CenterFrequency in range(690,700):
                bin = 0
            elif CenterFrequency in range(720,730):
                bin = 1
            elif CenterFrequency in range(750,760):
                bin = 2
            elif CenterFrequency in range(780,790):
                bin = 3
            elif CenterFrequency in range(810,820):
                bin = 4
            elif CenterFrequency in range(840,850):
                bin = 5
            elif CenterFrequency in range(870,880):
                bin = 6
            elif CenterFrequency in range(900,910):
                bin = 7
            else :
                raise ce.DataError('The center frequency does not match expected')
    # Need to figure out how to list the polarization values into the needed
    # matrix so that I can perform the conversion 
            STOKES = sp.mat(Data.data[time_index,:,cal_index,:])
            MUELLER = sp.mat(m_total[:,:,bin])
    # Next there is a matrix multiplication that will generate 
    # a new set of stokes values.
       	    stokesmod = MUELLER*STOKES
    # Need to check that the matrix multiplication worked, should 
    # now have a matrix with dimensions [nfreq,4] where nfreq is the 
    # number of frequencies in that bin (aka same dim as original stokes)
            for i in range(0,Data.dims[1]):
                for j in range(0,Data.dims[3]):
                    Data.data[time_index,i,cal_index,j] = stokesmod[i,j]	

# At this point the polarization values should be adjusted. 
# Now want to plot the polarizations I, Q, U, V as a function of Frequency.
# So I will want 4 plots for each cal_index, time_index. Not sure how to read
# the Data.data file so that I can generate those plots. Will need to also 
# generate plots for the pre-mueller data. 

# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
