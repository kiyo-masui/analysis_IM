"""
Procedure to calculate the Mueller parameters for each frequency from on-off scans of a calibrator such as 3C286."""
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT

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
               "Guppi_test" : False,
               }
prefix = 'mg_'


class MuellerGen(object) :
    """Calculates mueller matrix parameters for each frequency using on-off scans at three separate parallactic angles. 
    """
   
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
# Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                          prefix=prefix)
        self.feedback = feedback

    def peval(self,p):
        dG = p[0]
#        al = 90*sp.pi/180 # Note that this parameter is set based on observed results with chi = 0.
        al = ((p[1]+180)%360-180)*sp.pi/180
        ps = ((p[2]+180)%360-180)*sp.pi/180
        ph = ((p[3]+180)%360-180)*sp.pi/180
        ep = p[4]
        Q = p[5]
        U = p[6]
#        ch = 0
        ch = ((p[7]+45)%180-45)*sp.pi/180
        theta = self.theta
        t = self.function
        for i in range(0,len(t)-1,3):
            t[i] = 0.5*dG + Q*(sp.cos(2.*al)*sp.cos(2.*theta[i])-sp.sin(2*al)*sp.cos(ch)*sp.sin(2*theta[i]))+U*(sp.cos(2*al)*sp.cos(2*theta[i])+sp.sin(2*al)*sp.cos(ch)*sp.sin(2*theta[i]))
            t[i+1] = 2*ep*sp.cos(ps+ph)+Q*(-sp.sin(2*al)*sp.cos(ps+ch)*sp.cos(2.*theta[i])-(sp.cos(al)*sp.cos(al)*sp.cos(ps)-sp.sin(al)*sp.sin(al)*sp.cos(ps+2*ch))*sp.sin(2*theta[i]))+U*(-sp.sin(2*al)*sp.cos(ps+ch)*sp.cos(2*theta[i])+(sp.cos(al)*sp.cos(al)*sp.cos(ps)-sp.sin(al)*sp.sin(al)*sp.cos(ps+2*ch))*sp.sin(2*theta[i]))
            t[i+2] = 2*ep*sp.sin(ps+ph)+Q*(-sp.sin(2*al)*sp.sin(ps+ch)*sp.cos(2*theta[i])-(sp.sin(ps)*sp.cos(al)*sp.cos(al)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2*ch))*sp.sin(2*theta[i]))+U*(-sp.sin(2*al)*sp.sin(ps+ch)*sp.cos(2*theta[i])+(sp.sin(ps)*sp.cos(al)*sp.cos(al)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2*ch))*sp.sin(2*theta[i]))
        t[len(t)-1] = U/Q
        return t 

    def residuals(self, p, d, errors, f):
        err = (d[:,f]-self.peval(p))/errors
        return err
    
    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        guppi_result = params['Guppi_test']

        self.file_num = len(params['file_middles']) # getting a variable for number of calibrator files being used

# Need to remove count for calibrator files that are not the right size.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle + 
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname)
            n_scans = len(Reader.scan_set)
            if guppi_result == True : 
                if n_scans != 2 :
                    self.file_num -=1
            elif guppi_result == False :
                if n_scans != 4 :
                    self.file_num -=1

# Need to know the general frequency binning (going to assume that it's 200 for guppi, 260 for spectrometer, aka 1 MHz binning)
        if guppi_result == True :
            freq_num = 200 
        if guppi_result == False :
            freq_num = 260    
            self.file_num *= 2 #because there are two sets of scans per file for spectrometer, need to double the number of values


        self.function = sp.zeros(3*self.file_num+1) #setting a variable of the right size for peval to use
        self.theta = sp.zeros(3*self.file_num+1) #setting a variable for parallactic angle in radians
        d = sp.zeros((3*self.file_num+1,freq_num)) #setting a variable for measured values
        
# Loop over files to process.      
        k=0
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
# Read in the data, and loop over data blocks.
            Reader = core.fitsGBT.Reader(input_fname)
            n_IFs = len(Reader.IF_set) # Should be 1 given that we've stitched windows for the spectrometer or by def in guppi
            n_scans = len(Reader.scan_set) #Should be 4 for the spectrometer, 2 for guppi
            OnBlocks = Reader.read(range(0,n_scans,2),0,force_tuple=True)
            OffBlocks = Reader.read(range(1,n_scans,2),0,force_tuple=True)
#force_tuple=True makes the ouput of Reader.read a tuple even if thre is only one Block to return.
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)

# Setting labels for indices for later
            on_ind = 0
            off_ind = 1
            I_ind = 0
            Q_ind = 1
            U_ind = 2
            V_ind = 3
            
#Calculating Parallactic angles for the cal file
            PA = sp.zeros(n_scans)
            m = 0           
            for Data in Blocks:
                freq_len = Data.dims[3]
                Data.calc_freq()
                freq_val = Data.freq
                freq_val = freq_val/1000000       
                Data.calc_PA()
                PA[m] = ma.mean(Data.PA)  
                m+=1
            
            if guppi_result == False :
                if n_scans == 4 : # To make sure that we don't use incomplete data
                    self.theta[k] = 0.5*(PA[0]+PA[1]) # the average between the on and off values
                    self.theta[k+1] =0.5*(PA[0]+PA[1])
                    self.theta[k+2] = 0.5*(PA[0]+PA[1])
                    self.theta[k+3] = 0.5*(PA[2]+PA[3]) 
                    self.theta[k+4] = 0.5*(PA[2]+PA[3])
                    self.theta[k+5] = 0.5*(PA[2]+PA[3])

                    S_med_on = sp.zeros((2,freq_len,4))
                    S_med = sp.zeros((2,freq_len,4))  
            
                    i=0  
                    for Data in OnBlocks :
                        S_med_on[i,:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med_on[i,:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med_on[i,:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
                        S_med_on[i,:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                        i+=1
             
                    j=0
                    for Data in OffBlocks :
                        S_med[j,:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med[j,:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med[j,:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
                        S_med[j,:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                        j+=1

                    I_onoff_1 = S_med_on[0,:,0]-S_med[0,:,0] # for first set of on and off scans
                    I_onoff_2 = S_med_on[1,:,0]-S_med[1,:,0] # for second set of on and off scans

# Setting the measured stokes values for each file (to be used in the least squares fit)
                    d[k,:] = (S_med_on[0,:,1]-S_med[0,:,1])/I_onoff_1
                    d[k+1,:] = (S_med_on[0,:,2]-S_med[0,:,2])/I_onoff_1
                    d[k+2,:] = (S_med_on[0,:,3]-S_med[0,:,3])/I_onoff_1
                    d[k+3,:] = (S_med_on[1,:,1]-S_med[1,:,1])/I_onoff_2 
                    d[k+4,:] = (S_med_on[1,:,2]-S_med[1,:,2])/I_onoff_2
                    d[k+5,:] = (S_med_on[1,:,3]-S_med[1,:,3])/I_onoff_2
                    k+=6            

            if guppi_result == True : #This is the same as above only there is a single set of on and off scans in this case.
                if n_scans == 2 : 
                    self.theta[k] = ma.mean(PA)
                    self.theta[k+1] = ma.mean(PA)
                    self.theta[k+2] = ma.mean(PA)

                    S_med_on = sp.zeros((freq_len,4))
                    S_med = sp.zeros((freq_len,4))

                    for Data in OnBlocks:
                        S_med_on[:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med_on[:,1]  = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med_on[:,2]  = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med_on[:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                    
                    for Data in OffBlocks :
                        S_med[:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med[:,1]  = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med[:,2]  = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med[:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)

                    I_onoff = S_med_on[:,0]-S_med[:,0]
 
                    d[k,:] = (S_med_on[:,1]-S_med[:,1])/I_onoff
                    d[k+1,:] = (S_med_on[:,2]-S_med[:,2])/I_onoff
                    d[k+2,:] = (S_med_on[:,3]-S_med[:,3])/I_onoff
                    k+=3

        d[k,:] = sp.tan(2*33*sp.pi/180) # the 33 degrees here  is specific to 3C286, takes a different value for different calibrators
        #The seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], Qsrc[5], Usrc[6] chi[7] => the parameter vector is p
        p0 = [0.3, 90.0, 170.0, 10.0, 0.016, 0.005, 0.026, 0] # preliminary values based on guesses and heiles generation.
        error = sp.ones(3*self.file_num+1)
        #Note that error can be used to weight the equations if not all set to one.

        p_val_out = sp.zeros((freq_len, 9))
        p_err_out = sp.zeros((freq_len, 9))
     
        for f in range(0,freq_len):   
            plsq = leastsq(self.residuals,p0,args=(d,error,f),full_output=1, maxfev=100000,factor =1)
            pval = plsq[0] # this is the 1-d array of results
            perr = plsq[1] # this is a 2d array representing the estimated covariance of the results.

#want to adjust results if angles not in limits
            pval[1]=(pval[1]+180)%360-180
            pval[2]=(pval[2]+180)%360-180
            pval[3]=(pval[3]+180)%360-180
            pval[7]=(pval[7]+180)%360-180

            p_val_out[f,0] = freq_val[f]
            p_val_out[f,1] = pval[0]
            p_val_out[f,2] = pval[1]
            p_val_out[f,3] = pval[2]
            p_val_out[f,4] = pval[3]
            p_val_out[f,5] = pval[4]
            p_val_out[f,6] = pval[5]
            p_val_out[f,7] = pval[6]
            p_val_out[f,8] = pval[7]

# Gets error values, note that this doesn't work if set chi to zero. 
            p_err_out[f,0] = freq_val[f]
            p_err_out[f,1] = perr[0,0]
            p_err_out[f,2] = perr[1,1]
            p_err_out[f,3] = perr[2,2]
            p_err_out[f,4] = perr[3,3]
            p_err_out[f,5] = perr[4,4]
            p_err_out[f,6] = perr[5,5]
            p_err_out[f,7] = perr[6,6]
            p_err_out[f,8] = perr[7,7]

        np.savetxt('mueller_params_calc.txt', p_val_out, delimiter = ' ')
        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

