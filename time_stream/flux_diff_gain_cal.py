"""Module that performs flux and differential gain calibration. Needs to be run in XY basis"""

import scipy as sp
import numpy.ma as ma
import pylab as pl
import numpy as np
from scipy import optimize

import kiyopy.custom_exceptions as ce
import base_single
import map.tools

from core import fits_map
import utils.misc as utils
#from core import utils
import time
import ephem

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
    prefix = 'fgc_'
    # These are the parameters that should be read from file.  These are in
    # addition to the ones defined at the top of base_single.py.
    params_init = {
                   'mueler_file' : 'default_fname', 
                   'RM_file' : 'default_fname',
                   'R_to_sky' : True,
                   'DP_correct' : False,  
                   'RM_correct' : False,
                   'Flux_special' : False,
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
        i = self.file_ind
        file_middle = self.params['file_middles'][i]
        separate = file_middle.split('/')[1]
        sess_num = separate.split('_')[0]
        sess_num = int(sess_num)
        if (sess_num<10):
            sess_num = '0'+str(sess_num)
        project = file_middle.split('/')[0]

#        fg_file_name = self.params['mueler_file']+project+'/'+str(sess_num)+'_diff_gain_calc.txt'
        if self.params['RM_correct']==True:
            fg_file_name = self.params['mueler_file']+project+'/'+str(sess_num)+'_diff_gain_calc_RM.txt'

#Alternative for average flux/differential gain calibration.
        if self.params['Flux_special']==True:
            fg_file_name = self.params['mueler_file']

        self.flux_diff = flux_dg(fg_file_name)
        RM_dir = self.params['RM_file']
        calibrate_pol(Data, self.flux_diff,RM_dir,self.params['R_to_sky'],self.params['DP_correct'],self.params['RM_correct'])
       	Data.add_history('Flux calibrated and Corrected for leakage.', 
               	         ('Gain file: ' + self.params['mueler_file'],))

       	return Data

def flux_dg(fg_file_name) :
    mp = np.loadtxt(fg_file_name)
#This is a file with first index being freq, second index being XX_gain, third index being YY_gain

    freq_limit = len(mp[:,0])

    m_total = sp.zeros((2,freq_limit),float)
    for i in range(0,freq_limit):
        m_total[0,i] = mp[i,1]
        m_total[1,i] = mp[i,2]
    return m_total

def calibrate_pol(Data, m_total,RM_dir,R_to_sky,DP_correct,RM_correct) :
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
    if (Data.field['CRVAL4'][0] != -5 or Data.field['CRVAL4'][1] != -7 or
        Data.field['CRVAL4'][2] != -8 or Data.field['CRVAL4'][3] != -6) :
       	raise ce.DataError('Expected the polarization basis to be XY.')

    # A useful function that might need:
    Data.calc_freq()
    # Now data has an atribute Data.freq which is an array that gives the
    # frequency along the last axis.
          
    # Data.field['CRVAL1'] is center frequency in Hz. 
    # Data.data 4 dim array 2nd index polarization, 4th index frequency. 

    # Need to get parallactic angle:
#    Data.calc_PA()
    # This gives an array (Data.PA) of PA values of length = time dim.
# But gives the wrong answer. New code in misc.utils has to be built in for each time step.


#This segment of the code is for Rotation Measure Component
# Since the RM Tables have half hour time divisions and scans are shorter, we can do 1 selection.

    if RM_correct==True:
        Comp_Time = 0.0
        Full_date = Data.field['DATE-OBS'][Data.dims[0]/2]
        Date = Full_date.split('T')[0]
        Year = Date.split('-')[0]
        Month = Date.split('-')[1]
        Day = Date.split('-')[2]
        Full_time = Full_date.split('T')[1]
        Hour = Full_time.split(':')[0]
        Min = Full_time.split(':')[1]
        Sec = Full_time.split(':')[2]
        if int(Min)<=15:
            Comp_Time = float(Hour)+0.0
        elif int(Min)<=45:
            Comp_Time = float(Hour)+0.5
        else :
            Comp_Time = float(Hour)+1.0
    #Victor's tables have time in format Hour (xx.xx), Az (deg), El (deg), RM
    # Angle phi = RM*(wavelength)^2 where phi is in radians and wavelength is in meters

        RM_file_name = RM_dir + Year + Month + Day + '_RM.txt'
        RM_data = np.loadtxt(RM_file_name)
        RA_RM = sp.zeros(len(RM_data[:,0]))
        DEC_RM = sp.zeros(len(RM_data[:,0]))
        for i in range(0,len(RM_data[:,0])):
            RM_Hr = int(RM_data[i,0])
            if RM_data[i,0]%1 == 0 :
                RM_Min = '00'
                minutes = 0.0
            else:
                RM_Min = '30'
                minutes = 0.5
            Test = float(RM_Hr)+minutes
            if str(Comp_Time) == str(Test):
                UT_RM = Year+'-'+Month+'-'+Day+'T'+str(RM_Hr)+':'+RM_Min+':00.00'
                EL_RM = RM_data[i,2]
                AZ_RM = RM_data[i,1]
                RA_RM[i], DEC_RM[i] = utils.elaz2radecGBT(EL_RM,AZ_RM,UT_RM)
    #Now have tables of RA/DEC to compare to actual RA/DEC
        RM = 0


#This segment of the code is for Differential Phase Correction Generation

# Can determine the differential phase prior to the loop:
    if DP_correct==True:
# Set up a table of data to examine (calon-caloff to get Tcal)
        Tcal = ma.mean(Data.data[:,:,0,:],axis=0)-ma.mean(Data.data[:,:,1,:],axis=0)

# This segment of the code is for solving differential phase for each scan.
        fitfunc = lambda p,x: sp.cos(p[0]*x+p[1])
        errfunc = lambda p,x,y: fitfunc(p,x)-y
        freqs = sp.zeros(Data.dims[3])
        U_data = sp.zeros(Data.dims[3])
        V_data = sp.zeros(Data.dims[3])
        R_data = sp.zeros(Data.dims[3])
        for i in range(0,Data.dims[3]):
            freqs[i] = Data.freq[Data.dims[3]-i-1]/1000000
            U_data[i] = Tcal[1,Data.dims[3]-i-1]
            V_data[i] = Tcal[2,Data.dims[3]-i-1]
            R_data[i] = U_data[i]/sp.sqrt(U_data[i]**2+V_data[i]**2)

        mask_num = 0
        mask_num2 = -1
        for j in range(0,Data.dims[3]):
            if int(freqs[j])==710:
                mask_num = j
            if int(freqs[j])==740:
                mask_num2 = j
        print mask_num,mask_num2,Data.dims[3]
        Datain = R_data[mask_num:mask_num2]
        fin = freqs[mask_num:mask_num2]
        bad_pts = np.logical_or(np.isnan(Datain),np.isinf(Datain))
        good_ind = np.where(np.logical_not(bad_pts))
        Datain = Datain[good_ind]
        fin = fin[good_ind]
        R0 = [0.18,1.0]
        if len(good_ind[0])>1:
            R,success = optimize.leastsq(errfunc,R0[:],args=(fin,Datain),maxfev=10000)
            R[1] = R[1]%(2*sp.pi)
            print R
        else:
            R=[0.0,0.0]
            print "Not able to resolve a noise cal phase, setting phase to zero."
  

# This starts the actual data processing for the given scan
    RM_set = []     
    for time_index in range(0,Data.dims[0]):
        AZ = Data.field['CRVAL2'][time_index]
        EL = Data.field['CRVAL3'][time_index]
	UT = Data.field['DATE-OBS'][time_index]
	PA = utils.azel2pGBT(AZ,EL,UT)
#Extra data needed for RM Correction
#	if RM_correct==True:
#	    RA = AZ
#	    DEC = EL

# Extra data needed for Rotation Measure Correction
        if RM_correct==True:
#            RA = Data.field['CRVAL2'][time_index]
#            DEC = Data.field['CRVAL3'][time_index]
            RA, DEC = utils.azel2radecGBT(AZ,EL,UT)
#            print RA, DEC
            RM = 0
            dtheta0 = sp.pi
            RA_comp = RA
            DEC_comp = DEC
            for i in range(0,len(RA_RM)):
                dtheta = ma.arccos(ma.sin(DEC*sp.pi/180.)*ma.sin(DEC_RM[i]*sp.pi/180.)+ma.cos(DEC*sp.pi/180)*ma.cos(DEC_RM[i]*sp.pi/180.)*ma.cos(abs(RA-RA_RM[i])*sp.pi/180))
                if dtheta<dtheta0:
                    dtheta0 = dtheta
                    RM = RM_data[i,3]
                    RA_comp = RA_RM[i]
                    DEC_comp = DEC_RM[i]
            
            RM_set.append(RM)
#            if RM==0:
#                print 'Scan RM match not found'
#                print 'Scan RA and DEC:',RA, DEC
#            else:
#                print 'Scan RA and DEC:',RA,DEC
#                print 'RMtable RA and DEC:',RA_comp,DEC_comp
#                print 'RM Correction',RM
#                print '-------------------------------'


#            valid = []
#            for i in range(0,len(RA_RM)):
#                if RA_RM[i] != 0:
#                    if abs(RA-RA_RM[i])<10.0:
#                        if abs(DEC-DEC_RM[i])<10.0:
#                            RM = RM_data[i,3]     
#                            valid.append(i)
#            RA_M = 10.0
#            DEC_M = 10.0
#            for j in range(0,len(valid)):
#                if abs(RA-RA_RM[valid[j]])<RA_M:
#                    if abs(DEC-DEC_RM[valid[j]])<DEC_M:
#                        RM = RM_data[valid[j],3]  
#            print RM
                         
 
#Generate a sky matrix for this time index (assumes a XY basis):
        m_sky = sp.zeros((4,4))
        m_sky[0,0] = 0.5*(1+ma.cos(2*PA*sp.pi/180))
        m_sky[0,1] = -ma.sin(2*PA*sp.pi/180)
        m_sky[0,3] = 0.5*(1-ma.cos(2*PA*sp.pi/180))
        m_sky[1,0] = 0.5*ma.sin(2*PA*sp.pi/180)
        m_sky[1,1] = ma.cos(2*PA*sp.pi/180)
        m_sky[1,3] = -0.5*ma.sin(2*PA*sp.pi/180)
        m_sky[2,2] = 1
        m_sky[3,0] = 0.5*(1-ma.cos(2*PA*sp.pi/180))
        m_sky[3,1] = ma.sin(2*PA*sp.pi/180)
        m_sky[3,3] = 0.5*(1+ma.cos(2*PA*sp.pi/180))

        M_sky = sp.mat(m_sky)
        M_sky = M_sky.I

#        for cal_index in range(0,Data.dims[2]):
        for cal_index in range(0,2):
# Determines the Gains to use   
            for freq in range(0,Data.dims[3]):
# Tells which mueller matrix to use. 
               freq_limit = len(m_total[0,:])
               frequency = int(Data.freq[freq]/1000)
               bin = int((900000-frequency)*freq_limit/200000)

# Generate a sky matrix for this time index:
# With faraday rotation  sky matrix now frequency dependent
               if RM_correct==True:
                   wavelength = 300000.0/float(frequency) # should be in meters given that frequency is in kHz
                   Phi = RM*wavelength*wavelength
                   m_sky = sp.zeros((4,4)) 
                   m_sky[0,0] = 0.5*(1+ma.cos(2*PA*sp.pi/180+Phi))
                   m_sky[0,1] = -ma.sin(2*PA*sp.pi/180+Phi) 
                   m_sky[0,3] = 0.5*(1-ma.cos(2*PA*sp.pi/180+Phi))
                   m_sky[1,0] = 0.5*ma.sin(2*PA*sp.pi/180+Phi) 
                   m_sky[1,1] = ma.cos(2*PA*sp.pi/180+Phi) 
                   m_sky[1,3] = -0.5*ma.sin(2*PA*sp.pi/180+Phi)
                   m_sky[2,2] = 1
                   m_sky[3,0] = 0.5*(1-ma.cos(2*PA*sp.pi/180+Phi))
                   m_sky[3,1] = ma.sin(2*PA*sp.pi/180+Phi)
                   m_sky[3,3] = 0.5*(1+ma.cos(2*PA*sp.pi/180+Phi))
  
                   M_sky = sp.mat(m_sky)
                   M_sky = M_sky.I 

# Converts files into vector format 
               XY_params = Data.data[time_index,:,cal_index,freq]       
# Next there is a matrix multiplication that will generate 
# a new set of xy values. (Differential gain correction)
               XY_params[0] = XY_params[0]*m_total[0,bin]
               XY_params[3] = XY_params[3]*m_total[1,bin]
               XY_params[1] = XY_params[1]*sp.sqrt(m_total[0,bin]*m_total[1,bin])
               XY_params[2] = XY_params[2]*sp.sqrt(m_total[0,bin]*m_total[1,bin])

# Add in correction for differential phase
               if DP_correct==True:
                   XY_params[1] = XY_params[1]*sp.cos(R[0]*frequency/1000+R[1])-XY_params[2]*sp.sin(R[0]*frequency/1000+R[1])
                   XY_params[2] = XY_params[1]*sp.sin(R[0]*frequency/1000+R[1])+XY_params[2]*sp.cos(R[0]*frequency/1000+R[1])

#Rotate to sky coordinates (and RM correct if set)
               if R_to_sky==True:
                   XY_params = np.dot(M_sky,XY_params)

#Write corrected data to the new file. 
               for i in range(0,Data.dims[1]):
                    Data.data[time_index,i,cal_index,freq] = XY_params[i]	

    if RM_correct==True:
        print 'Number of time_steps in this scan is: ', len(RM_set)
        print 'Average RM correction for this scan is: ', ma.mean(RM_set)


# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
