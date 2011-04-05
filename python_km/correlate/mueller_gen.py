"""
Procedure to calculate the Mueller parameters for each frequency from on-off scans of a calibrator such as 3C286."""
import os

from scipy.optimize import leastsq
import scipy as sp
import numpy.ma as ma

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

    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)

        first_iteration = True
        # Loop over files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            # Read in the data, and loop over data blocks.
            Reader = core.fitsGBT.Reader(input_fname)
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)

# If I'm understanding this correctly. Blocks are scans in our input files so we have to call the data in a loop. I need some clarification about the format of the data, I'll ask Kiyo specific questions once I get there.

            #PA = sp.zeros(len(Blocks))# Trying to get an array of zeros the size of the number of data block
            #theta = [] # Parallactic angle vector if entering manually
            
            for Data in Blocks :
                S_med = sp.zeros(6,range(0,Data.dims[3]),4)
                d_3 = sp.zeros(3,range(0,Data.dims[3]),3)
                freq_len = Data.dims[3]
                #Calculate Parallactic Angle
                #Data.calc_pointing()
                #RA = Data.ra
                #DEC = Data.dec
                #LST = #Still need to figure this one out.
                #LAT = 38
                #PA[Data] = ma.arctan((ma.sin(LST-RA))/(ma.cos(DEC)*ma.tan(LAT)-ma.sin(DEC)*ma.cos(LST-RA)))

                #Getting the Stokes medians for each scan. Need to figure out how to get on_scan - off_scan medians. 
                on_ind = 0
                off_ind = 1
                I_ind = 0
                Q_ind = 1
                U_ind = 2
                V_ind = 3
                
                Data.calc_freq()
                freq_val = Data.freq  
                for freq in range(0,Data.dims[3]):
                    I_med = ma.median(Data.data[:,I_ind,off_ind,freq])
                    Q_med = ma.median(Data.data[:,Q_ind,off_ind,freq])
                    U_med = ma.median(Data.data[:,U_ind,off_ind,freq]) 
                    V_med = ma.median(Data.data[:,V_ind,off_ind,freq])
                    #Need to pre build this S array. 
                    S_med[Data][freq] = [I_med,Q_med,U_med,V_med]
	            # Later will have to do on-off.

                    if Data == even :
#Note I'm making an assumption here that each of the onoff scans has as its first data block the on calibrator data, then as its second data block the off calibrator scan and only two data blocks total - Need to check with kevin, kiyo
                        I_onoff = S_med[Data][freq][0]-S_med[Data+1][freq][0]
                        Q_onoff_ratio = (S_med[Data][freq][1]-S_med[Data+1][freq][1])/I_onoff
                        U_onoff_ratio = (S_med[Data][freq][2]-S_med[Data+1][freq][2])/I_onoff
                        V_onoff_ratio = (S_med[Data][freq][3]-S_med[Data+1][freq][3])/I_onoff
                        d_3[Data/2][freq] = [Q_onoff_ratio, U_onoff_ratio, V_onoff_ratio]
                    else :


            d = sp.zeros(9,range(0,freq_len))
            p0 = sp.zeros(9,range(0,freq_len))           
            for freq in range(0,freq_len):
# Seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], Qsrc[5], Usrc[6] => the parameter vector is p
                p0[freq] = [] # preliminary values based on heiles generation.
                d[0][freq] = d_3[0][freq][0]
                d[1][freq] = d_3[0][freq][1]
                d[2][freq] = d_3[0][freq][2]
                d[3][freq] = d_3[1][freq][0]
                d[4][freq] = d_3[1][freq][1]
                d[5][freq] = d_3[1][freq][2]
                d[6][freq] = d_3[2][freq][0]
                d[7][freq] = d_3[2][freq][1]
                d[8][freq] = d_3[2][freq][2]
        
# t[:][freq] is the vector of equations that should equal d[:][freq] if the parameters are correct.
                value = [0,1,2,3,4,5,6,7,8]
    def t(p,value,freq):
        t[0][freq] = 0.5*p[0][freq] + p[5][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.cos(2*theta[0]*sp.pi/180)+p[6][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.sin(2*theta[0]*sp.pi/180)
        t[1][freq] = 2*p[4][freq]*ma.cos((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[0]*sp.pi/180) - ma.cos(p[2][freq]*sp.pi/180)*ma.sin(2*theta[0]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[0]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[0]*sp.pi/180))
        t[2][freq] = 2*p[4][freq]*ma.sin((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(-ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi/180)*ma.cos(2*theta[0]*sp.pi/180) - ma.sin(p[2][freq]*sp.pi/180)*ma.sin(2*theta[0]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[0]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi.180)*ma.cos(2*theta[0]*sp.pi/180))
        t[3][freq] = 0.5*p[0][freq] + p[5][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.cos(2*theta[1]*sp.pi/180)+p[6][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.sin(2*theta[1]*sp.pi/180)
        t[4][freq] = 2*p[4][freq]*ma.cos((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[1]*sp.pi/180) - ma.cos(p[2][freq]*sp.pi/180)*ma.sin(2*theta[1]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[1]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[1]*sp.pi/180)) 
        t[5][freq]= 2*p[4][freq]*ma.sin((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(-ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi/180)*ma.cos(2*theta[1]*sp.pi/180) - ma.sin(p[2][freq]*sp.pi/180)*ma.sin(2*theta[1]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[1]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi.180)*ma.cos(2*theta[1]*sp.pi/180))
        t[6][freq] = 0.5*p[0][freq] + p[5][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.cos(2*theta[2]*sp.pi/180)+p[6][freq]*ma.cos(2*p[1][freq]*sp.pi/180)*ma.sin(2*theta[2]*sp.pi/180)
        t[7][freq] = 2*p[4][freq]*ma.cos((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[2]*sp.pi/180) - ma.cos(p[2][freq]*sp.pi/180)*ma.sin(2*theta[2]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[2]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180)*ma.cos(2*theta[2]*sp.pi/180)) 
        t[8][freq] = 2*p[4][freq]*ma.sin((p[2][freq]+p[3][freq])*sp.pi/180)+p[5][freq]*(-ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi/180)*ma.cos(2*theta[2]*sp.pi/180) - ma.sin(p[2][freq]*sp.pi/180)*ma.sin(2*theta[2]*sp.pi/180))+p[6][freq]*(ma.sin(2*theta[2]*sp.pi/180)*ma.sin(p[2][freq]*sp.pi/180) - ma.sin(2*p[1][freq]*sp.pi/180)*ma.cos(p[2][freq]*sp.pi.180)*ma.cos(2*theta[2]*sp.pi/180))
        return t[:][freq]
    def residuals(p, d, value, errors, freq):
        err[:][freq] = (d[:][freq] - t(p,value,freq)/errors
        return err         
#Note that error can be used to weight the equations if not all set to one.
                plsq[freq] = leastsq(residuals,p0,args=(d,value,error,freq),full_output=1, maxfev=2000)
                p_val[freq] = plsq[0][freq]
                p_err[freq] = plsq[1][freq]

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()

