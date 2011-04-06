"""
Procedure to calculate the Mueller parameters for each frequency from on-off scans of a calibrator such as 3C286."""
import os

from scipy.optimize import leastsq
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
        al = p[1]*sp.pi/180
        ps = p[2]*sp.pi/180
        ph = p[3]*sp.pi/180
        ep = p[4]
        Q = p[5]
        U = p[6]
        theta = self.theta
        t = sp.zeros(9)
        t[0] = 0.5*dG + Q*ma.cos(2*al)*ma.cos(2*theta[0])+U*ma.cos(2*al)*ma.sin(2*theta[0])
        t[1] = 2*ep*ma.cos(ps+ph)+Q*(ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[0]) - ma.cos(ps)*ma.sin(2*theta[0]))+U*(ma.sin(2*theta[0])*ma.sin(ps) - ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[0]))
        t[2] = 2*ep*ma.sin(ps+ph)+Q*(-ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[0]) - ma.sin(ps)*ma.sin(2*theta[0]))+U*(ma.sin(2*theta[0])*ma.sin(ps) - ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[0]))
        
        t[3] = 0.5*dG + Q*ma.cos(2*al)*ma.cos(2*theta[1])+U*ma.cos(2*al)*ma.sin(2*theta[1])
        t[4] = 2*ep*ma.cos(ps+ph)+Q*(ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[1]) - ma.cos(ps)*ma.sin(2*theta[1]))+U*(ma.sin(2*theta[1])*ma.sin(ps) - ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[1])) 
        t[5] = 2*ep*ma.sin(ps+ph)+Q*(-ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[1]) - ma.sin(ps)*ma.sin(2*theta[1]))+U*(ma.sin(2*theta[1])*ma.sin(ps) - ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[1]))

        t[6] = 0.5*dG + Q*ma.cos(2*al)*ma.cos(2*theta[2])+U*ma.cos(2*al)*ma.sin(2*theta[2])
        t[7] = 2*ep*ma.cos(ps+ph)+Q*(ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[2]) - ma.cos(ps)*ma.sin(2*theta[2]))+U*(ma.sin(2*theta[2])*ma.sin(ps) - ma.sin(2*al)*ma.sin(ps)*ma.cos(2*theta[2])) 
        t[8] = 2*ep*ma.sin(ps+ph)+Q*(-ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[2]) - ma.sin(ps)*ma.sin(2*theta[2]))+U*(ma.sin(2*theta[2])*ma.sin(ps) - ma.sin(2*al)*ma.cos(ps)*ma.cos(2*theta[2]))
        return t 

    def residuals(self, p, d, errors, f):
        err = (d[:,f] - self.peval(p) )/errors
        return err
    
    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)

        first_iteration = True
        self.theta = [295*sp.pi/180,336*sp.pi/180,15*sp.pi/180]
 # Parallactic angle vector if entering manually => Should be in units of fraction of pi.

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
                       
            for Data in Blocks:
                freq_len = Data.dims[3]
                #Calculate Parallactic Angle
                #Data.calc_pointing()
                #RA = Data.ra
                #DEC = Data.dec
                #LST = #Still need to figure this one out.
                #LAT = 38
                #PA[Data] = ma.arctan((ma.sin(LST-RA))/(ma.cos(DEC)*ma.tan(LAT)-ma.sin(DEC)*ma.cos(LST-RA)))
                Data.calc_freq()
                freq_val = Data.freq
                freq_val = freq_val/1000000       
     
            S_med_on = sp.zeros((2,freq_len,4))
            S_med = sp.zeros((2,freq_len,4))  

            i=0  
            for Data in OnBlocks :
                I_med_on = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                Q_med_on = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                U_med_on = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
                V_med_on = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                S_med_on[i,:,0] = I_med_on
                S_med_on[i,:,1] = Q_med_on
                S_med_on[i,:,2] = U_med_on
                S_med_on[i,:,3] = V_med_on
                i+=1


            if n_scans>1 :
                I_on = 0.5*(S_med_on[0,:,0]+S_med_on[1,:,0])
                Q_on = 0.5*(S_med_on[0,:,1]+S_med_on[1,:,1])
                U_on = 0.5*(S_med_on[0,:,2]+S_med_on[1,:,2])
                V_on = 0.5*(S_med_on[0,:,3]+S_med_on[1,:,3])                    
            else :
                I_on = S_med_on[0,:,0]
                Q_on = S_med_on[0,:,1]
                U_on = S_med_on[0,:,2]
                V_on = S_med_on[0,:,3]
 
            j=0
            for Data in OffBlocks :
                I_med = ma.median(Data.data[:,I_ind,off_ind,:],axis=0)
                Q_med = ma.median(Data.data[:,Q_ind,off_ind,:],axis=0)
                U_med = ma.median(Data.data[:,U_ind,off_ind,:],axis=0) 
                V_med = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)

                S_med[j,:,0] = I_med
                S_med[j,:,1] = Q_med
                S_med[j,:,2] = U_med
                S_med[j,:,3] = V_med
                j+=1

            if n_scans>1 :
                I_off = 0.5*(S_med[0,:,0]+S_med[1,:,0])
                Q_off = 0.5*(S_med[0,:,1]+S_med[1,:,1])
                U_off = 0.5*(S_med[0,:,2]+S_med[1,:,2])
                V_off = 0.5*(S_med[0,:,3]+S_med[1,:,3])                    
            else :
                I_off = S_med[0,:,0]
                Q_off = S_med[0,:,1]
                U_off = S_med[0,:,2]
                V_off = S_med[0,:,3]
   
            I_onoff = I_on - I_off
#            print I_onoff
            Q_onoff_ratio = (Q_on - Q_off)/I_onoff
            U_onoff_ratio = (U_on - U_off)/I_onoff
            V_onoff_ratio = (V_on - V_off)/I_onoff

            if k == 0:
                d_3_0_0 = Q_onoff_ratio
                d_3_0_1 = U_onoff_ratio
                d_3_0_2 = V_onoff_ratio
            elif k == 1:
                d_3_1_0 = Q_onoff_ratio
                d_3_1_1 = U_onoff_ratio
                d_3_1_2 = V_onoff_ratio
            elif k == 2:
                d_3_2_0 = Q_onoff_ratio
                d_3_2_1 = U_onoff_ratio
                d_3_2_2 = V_onoff_ratio    
            k+=1
                        
        self.freq_len = freq_len
# The seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], Qsrc[5], Usrc[6] => the parameter vector is p
        p0 = [0.3, -2.0, 170.0, 10.0, 0.016, 0.005, 0.026] # preliminary values based on guesses and heiles generation.
#        print freq_len
        d = sp.zeros((9,freq_len))
        d[0,:] = d_3_0_0
        d[1,:] = d_3_0_1
        d[2,:] = d_3_0_2
        d[3,:] = d_3_1_0
        d[4,:] = d_3_1_1
        d[5,:] = d_3_1_2
        d[6,:] = d_3_2_0
        d[7,:] = d_3_2_1
        d[8,:] = d_3_2_2
#        print d
# t is the array of equations for each frequency that should equal d for each frequency if the parameters are correct.

        error = [1,1,1,1,1,1,1,1,1]
#Note that error can be used to weight the equations if not all set to one.
###############################################################################
#Debugged above this point.
        p_val_out = sp.zeros((freq_len, 8))
        p_err_out = sp.zeros((freq_len, 8))
        for f in range(0,freq_len):   
            plsq = leastsq(self.residuals,p0,args=(d,error,f),full_output=1, maxfev=2000)
            pval = plsq[0]
#            print pval
            perr = plsq[1]
#            print perr
            p_val_out[f,0] = freq_val[f]
            p_val_out[f,1] = pval[0]
            p_val_out[f,2] = pval[1]
            p_val_out[f,3] = pval[2]
            p_val_out[f,4] = pval[3]
            p_val_out[f,5] = pval[4]
            p_val_out[f,6] = pval[5]
            p_val_out[f,7] = pval[6]
            p_err_out[f,0] = freq_val[f]
            p_err_out[f,1] = perr[-1:,0]
            p_err_out[f,2] = perr[-1:,1]
            p_err_out[f,3] = perr[-1:,2]
            p_err_out[f,4] = perr[-1:,3]
            p_err_out[f,5] = perr[-1:,4]
            p_err_out[f,6] = perr[-1:,5]
            p_err_out[f,7] = perr[-1:,6]
#        print perr[-1:]
#        print p_val_out[40]
        np.savetxt('mueller_params_calc.txt', p_val_out, delimiter = ' ')
        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

