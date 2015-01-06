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
		   'OD_correct' : False,
		   'DP_linear': False,
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
#        Data.calc_freq()
#        frequency = Data.freq/1000000
#        print Data.data[0,0,1,:]
#        pl.plot(frequency,ma.median(Data.data[:,0,0,:],axis=0),label='XX-init')
#        pl.plot(frequency,Data.data[0,1,0,:]-Data.data[0,1,1,:],label='Q-init')
#        pl.plot(frequency,Data.data[0,2,0,:]-Data.data[0,2,1,:],label='U-init')
#        pl.plot(frequency,ma.median(Data.data[:,3,0,:],axis=0),label='YY-init')

        # Main Action
        i = self.file_ind
        file_middle = self.params['file_middles'][i]
        separate = file_middle.split('/')[1]
        sess_num = separate.split('_')[0]
        sess_num = int(sess_num)
        if (sess_num<10):
            sess_num = '0'+str(sess_num)
        project = file_middle.split('/')[0]
#        print sess_num

        fg_file_name = self.params['mueler_file']+project+'/'+str(sess_num)+'_new_fcal.txt'
#        fg_file_name = self.params['mueler_file']+project+'/'+str(sess_num)+'_diff_gain_calc.txt'
#        if self.params['RM_correct']==True:
#            fg_file_name = self.params['mueler_file']+project+'/'+str(sess_num)+'_diff_gain_calc_RM.txt'

#Alternative for average flux/differential gain calibration.
        if self.params['Flux_special']==True:
            fg_file_name = self.params['mueler_file']+'all_one.txt'
#        fg_file_name = self.params['mueler_file']+'1hr_fdg_calc_avg.txt'
        self.flux_diff = flux_dg(fg_file_name)
        RM_dir = self.params['RM_file']
	Mpfn=self.params['mueler_file']+project+'/'+'Mueller_parameters.txt'
	self.J_params=np.loadtxt(Mpfn)
	Mp4fn=self.params['mueler_file']+project+'/'+'theta_0.txt'
	Mp56fn=self.params['mueler_file']+project+'/'+'eps_diff.txt'
	J4=np.loadtxt(Mp4fn)
	J56=np.loadtxt(Mp56fn)
	J4=np.array(J4)
	nfreq=J4.shape[0]-1
	#print int(sess_num)
	dist=np.array(np.where(J4[0]<int(sess_num))).shape[1]-1
	theta_0=J4[1:nfreq+1,dist]
	eps_d=J56[1:nfreq+1,dist]
	eps_sum=self.J_params[:,4]+self.J_params[:,5]
	#print self.J_params[:,3].shape, theta_0.shape
	self.J_params[:,3]=theta_0
	self.J_params[:,4]=eps_sum-0.5*eps_d
	self.J_params[:,5]=eps_sum+0.5*eps_d
	
	DPfn=self.params['mueler_file']+project+'/'+'calibrated_UV_111019.txt'
	self.ncal_phase=np.loadtxt(DPfn)/180.0*np.pi
	
	#print self.params['OD_correct']
	#print self.J_params
        calibrate_pol(Data, self.flux_diff,RM_dir,self.params['R_to_sky'],self.params['DP_correct'],self.params['RM_correct'],self.params['OD_correct'],self.params['DP_linear'],self.J_params,self.ncal_phase)
       	Data.add_history('Flux calibrated and Corrected for leakage.', 
               	         ('Gain file: ' + self.params['mueler_file'],))
        
#        pl.plot(frequency,ma.median(Data.data[:,0,0,:],axis=0),label='XX-mod')
#        pl.plot(frequency,Data.data[0,1,0,:]-Data.data[0,1,1,:],label='Q-mod')
#        pl.plot(frequency,Data.data[0,2,0,:]-Data.data[0,2,1,:],label='U-mod')
#        pl.plot(frequency,ma.median(Data.data[:,3,0,:],axis=0),label='YY-mod')
#        pl.plot(frequency,19.74748409*pow((750.0/frequency),0.49899785)*(2.28315426-0.000484307905*frequency),label='3C286-Isrc')
#        pl.plot(frequency,25.15445092*pow((750.0/frequency),0.75578842)*(2.28315426-0.000484307905*frequency),label='3C48-Isrc')
#        pl.legend()
#        pl.ylim(-20,130)
#        pl.xlabel("Frequency (MHz)")
#        pl.ylabel("Sample Data")
#        title0 = str(sess_num)+'_'+str(Data.field['SCAN'])+'_caloff_flux_diff_gain_'
#        pl.savefig(title0+'Comparison_Test.png')
#        pl.clf()

       	return Data

# m_total is the final Mueller matrix that needs to multiply by the stokes parameters.

# Note: I have it set up so that the third index represents the frequency 
# bin where 0 is the lowest frequency bin and 7 is the highest, and the 
# first and second indices represent the mueller matrix for each frquency bin. 

def flux_dg(fg_file_name) :
    mp = np.loadtxt(fg_file_name)
#This is a file with first index being freq, second index being XX_gain, third index being YY_gain

    freq_limit = len(mp[:,0])

    m_total = sp.zeros((2,freq_limit),float)
    for i in range(0,freq_limit):
        m_total[0,i] = mp[i,1]
        m_total[1,i] = mp[i,2]
    return m_total

def calibrate_pol(Data, m_total,RM_dir,R_to_sky,DP_correct,RM_correct,OD_correct,DP_linear,J_params,ncal_phase):
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
    #Data.calc_PA()
    # This gives an array (Data.PA) of PA values of length = time dim.
#   print Data.dims[0]



# This segment of the code is for Rotation Measure Component

    # Since the RM Tables have half hour time divisions and scans are shorter, we can do 1 selection.

    if RM_correct==True:
        Comp_Time = 0.0
        Full_date = Data.field['DATE-OBS'][Data.dims[0]/2]
        #print Full_date
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
        if Comp_Time==24:
            Comp_Time=23.5
        #print Comp_Time
    #Victor's tables have time in format Hour (xx.xx), Az (deg), El (deg), RM
    # Angle phi = RM*(wavelength)^2 where phi is in radians and wavelength is in meters

        RM_file_name = RM_dir + Year + Month + Day + '_RM.txt'
        #print RM_file_name
        RM_data = np.loadtxt(RM_file_name)
        AZ_RM = sp.zeros(len(RM_data[:,0]))
        EL_RM = sp.zeros(len(RM_data[:,0]))
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
                EL_RM[i] = RM_data[i,2]
                AZ_RM[i] = RM_data[i,1]
#                EL_RM = RM_data[i,2]
#                AZ_RM = RM_data[i,1]
#                RA_RM[i], DEC_RM[i] = utils.elaz2radecGBT(EL_RM,AZ_RM,UT_RM)
    #Now have tables of RA/DEC to compare to actual RA/DEC
        RM = 0




#This segment of the code is for Differential Phase Correction Generation

# Can determine the differential phase prior to the loop:
    Tcal = ma.mean(Data.data[:,:,0,:],axis=0)-ma.mean(Data.data[:,:,1,:],axis=0)
    #print (np.angle(Tcal[1,242:203:-1]+1j*Tcal[2,242:203:-1])/np.pi*180)%360
    Tcals=Tcal
    Tcals[0,:]=(Tcal[0,:]+Tcal[3,:])*0.5
    Tcals[3,:]=(-Tcal[0,:]+Tcal[3,:])*0.5
    #print Tcals[0,150:160]
    #print (Tcal[0,:]+Tcal[3,:])*0.5
    #print Tcal.shape
    if DP_correct==True:
# Set up a table of data to examine (calon-caloff to get Tcal)
        
#    Tcal = ma.mean(Data.data[:,:,0,:]-Data.data[:,:,1,:],axis=0)



# This version was if we arbitrarily set 4 possible phases and found closest match. There seems to be
# enough variability in the phase within the four categories that this doesn't quite work.

    # Randomly pick frequency bin near one of the zero crossings to compare U's
#    U_test = Tcal[1,230]/sp.sqrt(Tcal[1,230]**2+Tcal[2,230]**2)
#    print Tcal[:,191]
#    print Tcal[1,:]
#    U_test = Tcal[1,191]
#    print U_test
#    chi_sq =sp.zeros(4)
#    dp_dat = sp.zeros((4,2))
#    dp_dat[0] = [0.1354,2.341]
#    dp_dat[1] = [0.0723, 2.4575]
#    dp_dat[1] = [0.0730,2.611] #calculated specifically for sess 81
#    dp_dat[2] = [0.1029,0.045]
#    dp_dat[3] = [0,0]
#    dp_dat[3] = [0.1669,5.609] # giving problems because closer for sess 81 at given freq
#    min = 10
#    val = 5
#    for i in range(0,4):
#       chi_sq[i] = U_test-sp.cos(dp_dat[i,0]*Data.freq[230]/1000000+dp_dat[i,1])/sp.sqrt(Tcal[1,230]**2+Tcal[2,230]**2)
#        print sp.cos(dp_dat[i,0]*Data.freq[191]/1000000+dp_dat[i,1])
#        chi_sq[i] = U_test-sp.cos(dp_dat[i,0]*Data.freq[191]/1000000+dp_dat[i,1])
#        if abs(chi_sq[i]) < min:
#            min = abs(chi_sq[i])
#            val = i
# val tells which of the correction functions to use.    
#    print chi_sq
#    print val
#    print Data.freq[191]



# Alternate code for solving differential phase for each scan.
        if DP_linear==True:
            fitfunc = lambda p,x: sp.cos(p[0]*x+p[1])+1j*sp.sin(p[0]*x+p[1])
            errfunc = lambda p,x,y: np.abs(fitfunc(p,x)-y)
            freqs = sp.zeros(Data.dims[3])
            U_data = sp.zeros(Data.dims[3])
            V_data = sp.zeros(Data.dims[3])
            R_data = sp.zeros(Data.dims[3])+1j*sp.zeros(Data.dims[3])
            for i in range(0,Data.dims[3]):
                freqs[i] = Data.freq[Data.dims[3]-i-1]/1000000
                U_data[i] = Tcal[1,Data.dims[3]-i-1]
                V_data[i] = Tcal[2,Data.dims[3]-i-1]
                R_data[i] = (U_data[i]+1j*V_data[i])/sp.sqrt(U_data[i]**2+V_data[i]**2)
#    print np.any(np.isnan(R_data))
#    print np.any(np.isinf(R_data))
#    print freqs    

            for j in range(0,Data.dims[3]):
                if int(freqs[j])==750:
                    mask_num = j
                if int(freqs[j])==850:
                    mask_num2 = j
            print mask_num,mask_num2
            Datain = R_data[mask_num:mask_num2]
	    #print Datain
            fin = freqs[mask_num:mask_num2]
            bad_pts = np.logical_or(np.isnan(Datain),np.isinf(Datain))
            good_ind = np.where(np.logical_not(bad_pts))
            Datain = Datain[good_ind]
            fin = fin[good_ind]
	    ngi=Datain.shape[0]
	    hngi=int(ngi/2.0)
            R0=[0,0]
	    R0[0]=np.angle(Datain[5]/Datain[0])/(fin[5]-fin[0])
	    R0[1]=np.angle(Datain[hngi])-fin[hngi]*R0[0]
	    #print R0
            if len(good_ind[0])>1:
                #print good_ind[0]
                R,success = optimize.leastsq(errfunc,R0[:],args=(fin,Datain),maxfev=10000)
	        #print Datain
	        #print fin
	        #print R
                R[1] = R[1]%(2*sp.pi)
                print R
	        print np.mean(errfunc(R,fin,Datain)**2.0)
	        #print ((R[0]*fin+R[1])/np.pi*180)%360
            else:
                R=[0.0,0.0]
                print "Not able to resolve a noise cal phase, setting phase to zero."
        #masn1=255-mask_num
        #masn2=255-mask_num2
        #VV=np.mean(Data.data[:,2,0,masn2:masn1]-Data.data[:,2,1,masn2:masn1],axis=0)
        #UU=np.mean(Data.data[:,1,0,masn2:masn1]-Data.data[:,1,1,masn2:masn1],axis=0)
        #print 'ang_b4',np.angle(UU+1j*VV)
        else:
            ang_change=np.angle(Tcal[1,:]+1j*Tcal[2,:])+ncal_phase

    if OD_correct==1:
        invm=np.zeros((Data.dims[3],4,4))
        for freq in range(0,Data.dims[3]):
            G=J_params[freq][0]
            gamma=J_params[freq][1]
            psi=J_params[freq][2]
            theta_0=0
            theta_1=J_params[freq][3]
            eps_0=J_params[freq][4]
            eps_1=J_params[freq][5]
            
            #theta_1=0.01
            #eps_0=0
            #eps_1=0
            
            jones=parameters_to_jones(G,gamma,psi,theta_0,theta_1,eps_0,eps_1)
            invm[freq]=jones_to_mueller(jones)

# This starts the actual data processing for the given scan    
    for time_index in range(0,Data.dims[0]):
        #print time_index
        AZ = Data.field['CRVAL2'][time_index]
        EL = Data.field['CRVAL3'][time_index]
        UT = Data.field['DATE-OBS'][time_index]
        PA = utils.azel2pGBT(AZ,EL,UT)
        if time_index==100:
            print 'PA=',PA,Data.field['DATE-OBS'][time_index]

# Extra data needed for Rotation Measure Correction
        if RM_correct==True:
#            RA = Data.field['CRVAL2'][time_index]
#            DEC = Data.field['CRVAL3'][time_index]
#        print RA
#        print DEC
            AZ=np.mod(AZ,360.0)
            RM = 0
            valid = []
            for i in range(0,len(AZ_RM)):
                if AZ_RM[i] != 0:
                    if abs(AZ-AZ_RM[i])<10.0:
                        if abs(EL-EL_RM[i])<10.0:
                            RM = RM_data[i,3]     
                            valid.append(i)
            ABS_M = 15.0
            #print valid
            for j in range(0,len(valid)):
            #print len(AZ_RM)
            #print(AZ_RM[30000])
            #for j in range(0,len(AZ_RM)):
                azd=(AZ-AZ_RM[valid[j]])*np.cos(EL/180.0*np.pi)
                eld=EL-EL_RM[valid[j]]
                #azd=(AZ-AZ_RM[j])*np.cos(EL/180.0*np.pi)
                #eld=EL-EL_RM[j]
                abs_m=np.abs(azd+1j*eld)
                if abs_m<ABS_M:
                    RM = RM_data[valid[j],3]
                    #RM = RM_data[j,3]
                    ABS_M=abs_m
                    if time_index==100:
                        print AZ,EL,RM,AZ_RM[valid[j]],EL_RM[valid[j]],ABS_M				
         
            if time_index==100:
                print RM,Data.field['DATE-OBS'][time_index],AZ,EL
 
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
        #M_sky = M_sky.I
        #if time_index==100:
        #    print 2*PA,ma.cos(2*PA*sp.pi/180),ma.sin(2*PA*sp.pi/180)
        #    print M_sky

#        for cal_index in range(0,Data.dims[2]):
        for cal_index in range(0,2):
        # Determines the Gains to use   
            for freq in range(0,Data.dims[3]):
     # Tells which mueller matrix to use. 
               freq_limit = len(m_total[0,:])
               frequency = int(Data.freq[freq]/1000)
#               print frequency
               bin = int((900000-frequency)*freq_limit/200000)
#               print bin
#               if freq_limit == 200:
#                   bin = 900-frequency
#Not setup to work with spectrometer data.
#               elif freq_limit == 260:
#                   bin = 929-frequency
#               print bin

    #Generate a sky matrix for this time index:
    #With faraday rotation  sky matrix now frequency dependent
               if RM_correct==True:
                   wavelength = 300000.0/float(frequency) # should be in meters given that frequency is in kHz
#               print wavelength
                   Phi = RM*wavelength*wavelength
                   #if time_index==100 and freq==170:
                       #print RM,wavelength,Phi
                       #print 2*PA,(2*PA*sp.pi/180-2*Phi)/sp.pi*180
#               print Phi
                   m_sky = sp.zeros((4,4)) 
                   m_sky[0,0] = 0.5*(1+ma.cos(2*PA*sp.pi/180-2*Phi))
                   m_sky[0,1] = -ma.sin(2*PA*sp.pi/180-2*Phi) 
                   m_sky[0,3] = 0.5*(1-ma.cos(2*PA*sp.pi/180-2*Phi))
                   m_sky[1,0] = 0.5*ma.sin(2*PA*sp.pi/180-2*Phi) 
                   m_sky[1,1] = ma.cos(2*PA*sp.pi/180-2*Phi) 
                   m_sky[1,3] = -0.5*ma.sin(2*PA*sp.pi/180-2*Phi)
                   m_sky[2,2] = 1
                   m_sky[3,0] = 0.5*(1-ma.cos(2*PA*sp.pi/180-2*Phi))
                   m_sky[3,1] = ma.sin(2*PA*sp.pi/180-2*Phi)
                   m_sky[3,3] = 0.5*(1+ma.cos(2*PA*sp.pi/180-2*Phi))
  
                   M_sky = sp.mat(m_sky)
                   #M_sky = M_sky.I 
#               print M_sky 

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
		   U=XY_params[1]
		   V=XY_params[2]
		   if DP_linear==True:
		       r_ang=(R[0]*frequency/1000+R[1])
		   else:
		       r_ang=ang_change[freq]
		   XY_params[1] = U*sp.cos(r_ang)+V*sp.sin(r_ang)
		   XY_params[2] = -U*sp.sin(r_ang)+V*sp.cos(r_ang)
		   #if time_index==10 and cal_index==0 and freq%5==0:
                       #print frequency, (r_ang/np.pi*180)%360, (np.angle(Tcal[1,freq]+1j*Tcal[2,freq])/np.pi*180)%360,(np.angle((U+1j*V)/(XY_params[1]+1j*XY_params[2]))/np.pi*180)%360
		       #print r_ang%(2.0*np.pi),sp.cos(r_ang),sp.sin(r_ang),U+1j*V,XY_params[1]+1j*XY_params[2]
	       
	       if OD_correct==True:
		   inv_mueller=invm[freq]
                   XY_params = np.dot(inv_mueller,XY_params)
		    #if freq==10:
			#print G,gamma,psi,theta_0,theta_1,eps_0,eps_1
			#print mueller

    #Rotate to sky coordinates (and RM correct if set)

               if R_to_sky==True:
                   XY_params = np.dot(M_sky,XY_params)

    #Write corrected data to the new file. 

               for i in range(0,Data.dims[1]):
                    Data.data[time_index,i,cal_index,freq] = XY_params[i]
    #VV=np.mean(Data.data[:,2,0,masn2:masn1]-Data.data[:,2,1,masn2:masn1],axis=0)
    #UU=np.mean(Data.data[:,1,0,masn2:masn1]-Data.data[:,1,1,masn2:masn1],axis=0)
    #print 'ang_after',np.angle(UU+1j*VV)
    Tcal2 = ma.mean(Data.data[:,:,0,:],axis=0)-ma.mean(Data.data[:,:,1,:],axis=0)
    ang_diff=(np.angle((Tcal2[1,0:256:5]+1j*Tcal2[2,0:256:5])/(Tcal[1,0:256:5]+1j*Tcal[2,0:256:5]))/np.pi*180)%360
    #print 'b4',np.angle(Tcal[1,0:256:5]+1j*Tcal[2,0:256:5])/np.pi*180
    #print 'af',np.angle(Tcal2[1,0:256:5]+1j*Tcal2[2,0:256:5])/np.pi*180
    #print ang_diff



def parameters_to_jones(G,gamma,psi,theta_0,theta_1,eps_0,eps_1):
	inv=np.linalg.inv
	cos=np.cos
	sin=np.sin
	exp=np.exp
	dot=np.dot
	s0a=np.array([[cos(eps_0),1j*sin(eps_0)],[1j*sin(eps_0),cos(eps_0)]])
	s1a=np.array([[cos(eps_1),1j*sin(eps_1)],[1j*sin(eps_1),cos(eps_1)]])
	s0b=np.array([[cos(theta_0),sin(theta_0)],[-sin(theta_0),cos(theta_0)]])
	s1b=np.array([[cos(theta_1),sin(theta_1)],[-sin(theta_1),cos(theta_1)]])
	s0=dot(s0a,s0b)
	s1=dot(s1a,s1b)
	cm=np.array([s0[0],s1[1]])
	ext=gamma+1j*psi
	rs=np.array([[exp(ext),0],[0,exp(-ext)]])
	jones=G*dot(rs,cm)
	return jones

def jones_to_mueller(jones):
	inv=np.linalg.inv
	dot=np.dot
	conj=np.conjugate
	cjones=conj(jones)
	mup=np.kron(jones,cjones)
	liya=np.array([[0,0,0,2],[0,1,1,0],[0,1j,-1j,0],[2,0,0,0]])
	iliya=inv(liya)
	mueller=dot(liya,dot(mup,iliya))
	inv_mueller=np.real(inv(mueller))
	return inv_mueller
	

# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
