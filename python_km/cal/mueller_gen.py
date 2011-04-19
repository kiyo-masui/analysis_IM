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
        #Original Version
        dG = p[0]
        al = ((p[1]+180)%360-180)*sp.pi/180
        ps = ((p[2]+180)%360-180)*sp.pi/180
        ph = ((p[3]+180)%360-180)*sp.pi/180
        ep = p[4]
        Q = p[5]
        U = p[6]
        ch = ((p[7]+180)%360-180)*sp.pi/180
        theta = self.theta
        # Alternate values
#        dG = p[0]
#        Q = p[1]
#        U = p[2]
#        ep = p[3]
        t = self.function
        for i in range(0,len(t),3):
            t[i] = 0.5*dG + Q*(ma.cos(2.*al)*ma.cos(2.*theta[i])-ma.sin(2*al)*ma.cos(ch)*ma.sin(2*theta[i]))+U*(ma.cos(2*al)*ma.cos(2*theta[i])+ma.sin(2*al)*ma.cos(ch)*ma.sin(2*theta[i]))
#            t[i] = 0.5*dG+Q*(p[4]*ma.cos(2*theta[i])-p[5]*ma.sin(2*theta[i]))+U*(p[4]*ma.sin(2*theta[i])+p[5]*ma.cos(2*theta[i]))
            t[i+1] = 2*ep*ma.cos(ps+ph)+Q*(-ma.sin(2*al)*ma.cos(ps+ch)*ma.cos(2.*theta[i])-(ma.cos(al)*ma.cos(al)*ma.cos(ps)-ma.sin(al)*ma.sin(al)*ma.cos(ps+2*ch))*ma.sin(2*theta[i]))+U*(-ma.sin(2*al)*ma.cos(ps+ch)*ma.cos(2*theta[i])+(ma.cos(al)*ma.cos(al)*ma.cos(ps)-ma.sin(al)*ma.sin(al)*ma.cos(ps+2*ch))*ma.sin(2*theta[i]))
#            t[i+1] = 2*ep*p[6]+Q*(-p[7]*ma.cos(2*theta[i])-p[8]*ma.sin(2*theta[i]))+U*(-p[7]*ma.sin(2*theta[i])+p[8]*ma.cos(2*theta[i]))
            t[i+2] = 2*ep*ma.sin(ps+ph)+Q*(-ma.sin(2*al)*ma.sin(ps+ch)*ma.cos(2*theta[i])-(ma.sin(ps)*ma.cos(al)*ma.cos(al)-ma.sin(al)*ma.sin(al)*ma.sin(ps+2*ch))*ma.sin(2*theta[i]))+U*(-ma.sin(2*al)*ma.sin(ps+ch)*ma.cos(2*theta[i])+(ma.sin(ps)*ma.cos(al)*ma.cos(al)-ma.sin(al)*ma.sin(al)*ma.sin(ps+2*ch))*ma.sin(2*theta[i]))
#            t[i+2] = 2*ep*p[9]+Q*(-p[10]*ma.cos(2*theta[i])-p[11]*ma.sin(2*theta[i]))+U*(-p[10]*ma.sin(2*theta[i])+p[11]*ma.cos(2*theta[i]))
        return t 

    def residuals(self, p, d, errors, f):
        err = (d[:,f]-self.peval(p))/errors
#        print err
        return err
    
    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        guppi_result = params['Guppi_test']
        # Need to know the general frequency binning (going to assume that it's 200 for guppi, 260 for spectrometer, aka 1 MHz binning)
        if guppi_result == True :
            freq_num = 200
        if guppi_result == False :
            freq_num = 260
       
        self.file_num = len(params['file_middles']) # getting a variable for number of calibrator files being used
        if guppi_result == False :
            self.file_num *= 2 #because there are two sets of scans per file for spectrometer.
        self.function = sp.zeros(3*self.file_num) #setting a variable for peval
        self.theta = sp.zeros(3*self.file_num) #parallactic angle in radians
        d = sp.zeros((3*self.file_num,freq_num)) #measured values
        print self.file_num   
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

            on_ind = 0
            off_ind = 1
            I_ind = 0
            Q_ind = 1
            U_ind = 2
            V_ind = 3
            
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
    
#            print PA

            if guppi_result == False :
                self.theta[k] = 0.5*(PA[0]+PA[1])
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

                I_onoff_1 = S_med_on[0,:,0]-S_med[0,:,0]
                I_onoff_2 = S_med_on[1,:,0]-S_med[1,:,0]

                d[k,:] = (S_med_on[0,:,1]-S_med[0,:,1])/I_onoff_1
                d[k+1,:] = (S_med_on[0,:,2]-S_med[0,:,2])/I_onoff_1
                d[k+2,:] = (S_med_on[0,:,3]-S_med[0,:,3])/I_onoff_1
                d[k+3,:] = (S_med_on[1,:,1]-S_med[1,:,1])/I_onoff_2 
                d[k+4,:] = (S_med_on[1,:,2]-S_med[1,:,2])/I_onoff_2
                d[k+5,:] = (S_med_on[1,:,3]-S_med[1,:,3])/I_onoff_2
                k+=6

            if guppi_result == True :    
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

#        print d
#        print self.theta
#The seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], Qsrc[5], Usrc[6] chi[7] => the parameter vector is p
        p0 = [0.3, -2.0, 170.0, 10.0, 0.016, 0.005, 0.026, 1.0] # preliminary values based on guesses and heiles generation.

        #Alternate p0 if leave chi as variable and  avoid trig functions because of modulo 2pi issue.
#        p0 = [0.3, 0.005,0.026,0.016,0.997,-0.070,-1.0,0.069,-0.983,0,-0.012,0.173]

#        print freq_len

        error = sp.ones(3*self.file_num)
#Note that error can be used to weight the equations if not all set to one.

        p_val_out = sp.zeros((freq_len, 9))
        p_err_out = sp.zeros((freq_len, 9))
        #Alternate Version
#        p_val_out = sp.zeros((freq_len, 17))
#        p_err_out = sp.zeros((freq_len, 17))
      
        for f in range(0,freq_len):   
            plsq = leastsq(self.residuals,p0,args=(d,error,f),full_output=1, maxfev=100000,factor =10)
#            plsq = fmin_slsqp(self.residuals,p0,args=(d,error,f),bounds=[(-1,1),(0,360),(0,360),(0,360),(-1,1),(-1,1),(-1,1),(0,360)],full_output=1,iter=2000)
            pval = plsq[0]
#            print pval
            perr = plsq[1]
            print plsq[3]
# this is a 2d array representing the estimated covariance of the results.
#            print perr
#want to adjust results if angles not between +/- 180 
            while pval[1]>180:
                pval[1] -= 360
            while pval[1]<-180:
                pval[1] += 360 
            while pval[2]>180:
                pval[2] -= 360
            while pval[2]<-180:
                pval[2] += 360
            while pval[3]>180:
                pval[3] -= 360
            while pval[3]<-180:
                pval[3] += 360
            while pval[7]>180:
                pval[7] -= 360
            while pval[7]<-180:
                pval[7] += 360
            p_val_out[f,0] = freq_val[f]
            p_val_out[f,1] = pval[0]
            p_val_out[f,2] = pval[1]
#            print p_val_out[f,2]
            p_val_out[f,3] = pval[2]
#            print p_val_out[f,3]
            p_val_out[f,4] = pval[3]
#            print p_val_out[f,4]
            p_val_out[f,5] = pval[4]
            p_val_out[f,6] = pval[5]
            p_val_out[f,7] = pval[6]
            p_val_out[f,8] = pval[7]
#            print p_val_out[f,8]
            p_err_out[f,0] = freq_val[f]
            p_err_out[f,1] = perr[0,0]
            p_err_out[f,2] = perr[1,1]
            p_err_out[f,3] = perr[2,2]
            p_err_out[f,4] = perr[3,3]
            p_err_out[f,5] = perr[4,4]
            p_err_out[f,6] = perr[5,5]
            p_err_out[f,7] = perr[6,6]
            p_err_out[f,8] = perr[7,7]
#           Alternate version for different set of non-trig parameters
#            p_val_out[f,1] = pval[0]
#            p_err_out[f,1] = perr[0,0]
#            p_val_out[f,4] = pval[3]
#            p_err_out[f,4] = perr[3,3]
#            p_val_out[f,2] = pval[1]
#            p_err_out[f,2] = perr[1,1]
#            p_val_out[f,3] = pval[2]
#            p_err_out[f,3] = perr[2,2]
#            p_val_out[f,5] = pval[4]
#            print pval[4]
#            p_err_out[f,5] = perr[4,4]
#            p_val_out[f,6] = pval[5] 
#            print pval[5]
#            p_err_out[f,6] = perr[5,5] 
#            p_val_out[f,7] = pval[6]
#            print pval[6]
#            p_err_out[f,7] = perr[6,6] 
#            p_val_out[f,8] = pval[7]
#            print pval[7]
#            p_err_out[f,8] = perr[7,7]
#            p_val_out[f,9] = pval[8]
#            p_err_out[f,9] = perr[8,8]
#            p_val_out[f,10] = pval[9]
#            p_err_out[f,10] = perr[9,9]
#            p_val_out[f,11] = pval[10] 
#            p_err_out[f,11] = perr[10,10] 
#            p_val_out[f,12] = pval[11] 
#            p_err_out[f,12] = perr[11,11] 
#            alpha = ma.arccos(pval[4])/2
#            print alpha
#            p_val_out[f,13] = alpha*180/sp.pi
#            alpha_err = ma.arccos(perr[4,4])/2
#            p_err_out[f,13] = alpha_err*180/sp.pi
#            chi = ma.arccos(pval[5]/ma.sin(2*alpha))
#            p_val_out[f,14] = chi*180/sp.pi
#            chi_err = ma.arccos(perr[5,5]/ma.sin(2*alpha_err))
#            p_err_out[f,14] = chi_err*180/sp.pi
#            psi = ma.arccos(pval[7]/ma.sin(2*alpha)) - chi
#            p_val_out[f,15] = psi*180/sp.pi
#            psi_err = ma.arccos(perr[7,7]/ma.sin(2*alpha_err))-chi_err
#            p_err_out[f,15] = psi_err*180/sp.pi
#            phi = ma.arccos(pval[6])-psi
#            p_val_out[f,16] = phi*180/sp.pi
#            phi_err = ma.arccos(perr[6,6])-psi_err
#            p_err_out[f,16] = phi_err*180/sp.pi

#        print perr[-1:]
#        print p_val_out[40]
        np.savetxt('mueller_params_calc.txt', p_val_out, delimiter = ' ')
        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

