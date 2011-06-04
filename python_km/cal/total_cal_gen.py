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
prefix = 'tc_'


class MuellerGen(object) :
    """Calculates uncalibrated mueller matrix elements for each frequency using on-off scans at many separate parallactic angles. 
    """
   
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
# Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                          prefix=prefix)
        self.feedback = feedback

    def peval(self,p,f):
        theta = self.theta
        t = self.function
        d = self.d
        for i in range(0,len(t),4):
            t[i] = p[0]*d[i,f]+p[1]*d[i+1,f]+p[2]*d[i+2,f]+p[3]*d[i+3,f]
            t[i+1] =(p[4]*sp.cos(2*theta[i])-p[8]*sp.sin(2*theta[i]))*d[i,f]+(p[5]*sp.cos(2*theta[i])-p[9]*sp.sin(2*theta[i]))*d[i+1,f]+(p[6]*sp.cos(2*theta[i])-p[10]*sp.sin(2*theta[i]))*d[i+2,f]+(p[7]*sp.cos(2*theta[i])-p[11]*sp.sin(2*theta[i]))*d[i+3,f]
            t[i+2] =(p[4]*sp.sin(2*theta[i])+p[8]*sp.cos(2*theta[i]))*d[i,f]+(p[5]*sp.sin(2*theta[i])+p[9]*sp.cos(2*theta[i]))*d[i+1,f]+(p[6]*sp.sin(2*theta[i])+p[10]*sp.cos(2*theta[i]))*d[i+2,f]+(p[7]*sp.sin(2*theta[i])+p[11]*sp.cos(2*theta[i]))*d[i+3,f]
            t[i+3] =p[12]*d[i,f]+p[13]*d[i+1,f]+p[14]*d[i+2,f]+p[15]*d[i+3,f]
        return t 

    def residuals(self, p,errors, f,freq_val):
        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)
        PAsrc = 33.0*sp.pi/180.0
        Psrc = 0.07 # The fraction is about 7% for 3C286 at high frequencies, but not as stable at our frequencies.  
        Qsrc = Isrc*Psrc*sp.cos(2*PAsrc)
        Usrc = Isrc*Psrc*sp.sin(2*PAsrc)
        Vsrc = 0
        source = sp.zeros(4*self.file_num)
        for i in range(0,len(source),4):
            source[i] = Isrc
            source[i+1] = Qsrc
            source[i+2] = Usrc
            source[i+3] = Vsrc
        err = (source-self.peval(p,f))/errors
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


        self.function = sp.zeros(4*self.file_num) #setting a variable of the right size for peval to use
        self.theta = sp.zeros(4*self.file_num) #setting a variable for parallactic angle in radians
        self.d = sp.zeros((4*self.file_num,freq_num)) #setting a variable for measured values
        
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
            
# Going to skip the non guppi version for now.
#            if guppi_result == False :
#                if n_scans == 4 : # To make sure that we don't use incomplete data
#                    self.theta[k] = 0.5*(PA[0]+PA[1]) # the average between the on and off values
#                    self.theta[k+1] =0.5*(PA[0]+PA[1])
#                    self.theta[k+2] = 0.5*(PA[0]+PA[1])
#                    self.theta[k+3] = 0.5*(PA[2]+PA[3]) 
#                    self.theta[k+4] = 0.5*(PA[2]+PA[3])
#                    self.theta[k+5] = 0.5*(PA[2]+PA[3])

#                    S_med_on = sp.zeros((2,freq_len,4))
#                    S_med = sp.zeros((2,freq_len,4))  
            
#                    i=0  
#                    for Data in OnBlocks :
#                        S_med_on[i,:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
#                        S_med_on[i,:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
#                        S_med_on[i,:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
#                        S_med_on[i,:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
#                        i+=1
             
#                    j=0
#                    for Data in OffBlocks :
#                        S_med[j,:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
#                        S_med[j,:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
#                        S_med[j,:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
#                        S_med[j,:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
#                        j+=1

#                    I_onoff_1 = S_med_on[0,:,0]-S_med[0,:,0] # for first set of on and off scans
#                    I_onoff_2 = S_med_on[1,:,0]-S_med[1,:,0] # for second set of on and off scans

# Setting the measured stokes values for each file (to be used in the least squares fit)
#                    d[k,:] = (S_med_on[0,:,1]-S_med[0,:,1])/I_onoff_1
#                    d[k+1,:] = (S_med_on[0,:,2]-S_med[0,:,2])/I_onoff_1
#                    d[k+2,:] = (S_med_on[0,:,3]-S_med[0,:,3])/I_onoff_1
#                    d[k+3,:] = (S_med_on[1,:,1]-S_med[1,:,1])/I_onoff_2 
#                    d[k+4,:] = (S_med_on[1,:,2]-S_med[1,:,2])/I_onoff_2
#                    d[k+5,:] = (S_med_on[1,:,3]-S_med[1,:,3])/I_onoff_2
#                    k+=6            

            if guppi_result == True : #This is the same as above only there is a single set of on and off scans in this case.
                if n_scans == 2 : 
                    self.theta[k] = ma.mean(PA)
                    self.theta[k+1] = ma.mean(PA)
                    self.theta[k+2] = ma.mean(PA)
                    self.theta[k+3] = ma.mean(PA)

                    S_med_calon_src = sp.zeros((freq_len,4))
                    S_med_caloff_src = sp.zeros((freq_len,4))
                    S_med_calon = sp.zeros((freq_len,4))
                    S_med_caloff = sp.zeros((freq_len,4))

                    for Data in OnBlocks:
                        S_med_caloff_src[:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)

                        S_med_calon_src[:,0] = ma.median(Data.data[:,I_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,1] = ma.median(Data.data[:,Q_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,2] = ma.median(Data.data[:,U_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,3] = ma.median(Data.data[:,V_ind,on_ind,:],axis=0)
                    
                    for Data in OffBlocks:
                        S_med_caloff[:,0] = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                        S_med_caloff[:,1] = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                        S_med_caloff[:,2] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med_caloff[:,3] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
 
                        S_med_calon[:,0] = ma.median(Data.data[:,I_ind,on_ind,:],axis=0)
                        S_med_calon[:,1] = ma.median(Data.data[:,Q_ind,on_ind,:],axis=0)
                        S_med_calon[:,2] = ma.median(Data.data[:,U_ind,on_ind,:],axis=0)
                        S_med_calon[:,3] = ma.median(Data.data[:,V_ind,on_ind,:],axis=0)
 
                     
                    self.d[k,:] = 0.5*(S_med_calon_src[:,0]+S_med_caloff_src[:,0]-S_med_calon[:,0]-S_med_caloff[:,0])
                    self.d[k+1,:] = 0.5*(S_med_calon_src[:,1]+S_med_caloff_src[:,1]-S_med_calon[:,1]-S_med_caloff[:,1])
                    self.d[k+2,:] = 0.5*(S_med_calon_src[:,2]+S_med_caloff_src[:,2]-S_med_calon[:,2]-S_med_caloff[:,2])
                    self.d[k+3,:] = 0.5*(S_med_calon_src[:,3]+S_med_caloff_src[:,3]-S_med_calon[:,3]-S_med_caloff[:,3])
                    k+=4
        for a in range(0,4*self.file_num):
            for b in range(0,freq_num):
#                print self.d[a,b]
                if self.d[a,b] > 1000 :
                   self.d[a,b] = 1000

        #The seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], Qsrc[5], Usrc[6] chi[7] => the parameter vector is p
        p0 = [1.0, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 1.0] # guessed preliminary values
        error = sp.ones(4*self.file_num)
        #Note that error can be used to weight the equations if not all set to one.

        p_val_out = sp.zeros((freq_len, 17))
 #       p_err_out = sp.zeros((freq_len, 17))
     
        for f in range(0,freq_len):   
            plsq = leastsq(self.residuals,p0,args=(error,f,freq_val),full_output=0, maxfev=5000)
            pval = plsq[0] # this is the 1-d array of results
#            perr = plsq[1] # this is a 2d array representing the estimated covariance of the results. - Not working properly. 

            Mueller = sp.mat([[pval[0],pval[1],pval[2],pval[3]],[pval[4],pval[5],pval[6],pval[7]],[pval[8],pval[9],pval[10],pval[11]],[pval[12],pval[13],pval[14],pval[15]]])
#            Mueller = Mueller.I  
#            print Mueller

            p_val_out[f,0] = freq_val[f]
            p_val_out[f,1] = Mueller[0,0]
            p_val_out[f,2] = Mueller[0,1]
            p_val_out[f,3] = Mueller[0,2]
            p_val_out[f,4] = Mueller[0,3]
            p_val_out[f,5] = Mueller[1,0]
            p_val_out[f,6] = Mueller[1,1]
            p_val_out[f,7] = Mueller[1,2]
            p_val_out[f,8] = Mueller[1,3]
            p_val_out[f,9] = Mueller[2,0]
            p_val_out[f,10] = Mueller[2,1]
            p_val_out[f,11] = Mueller[2,2]
            p_val_out[f,12] = Mueller[2,3]
            p_val_out[f,13] = Mueller[3,0]
            p_val_out[f,14] = Mueller[3,1]
            p_val_out[f,15] = Mueller[3,2]
            p_val_out[f,16] = Mueller[3,3]

        np.savetxt('flux_mueller_matrix_calc.txt', p_val_out, delimiter = ' ')
#        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

