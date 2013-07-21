"""Procedure to calculate only the flux and differential gain between XX and YY for each frequency from on-off scans of a calibrator such as 3C286.

Currently setup so that it uses all the onoff scans from a session for a particular calibrator and gets the best fit. 
Run in analysis_IM: python cal/flux_diff_gain_cal_gen.py input/tcv/diff_gain_gen_guppi.ini
Note that the .ini file should indicate which session(s) and sourse you want to use. Script is run using data from a single source. The output is saved in my data directory under the folder diff_gain_params as a .txt file with three columns (freq, XXGain, YYGain).
Need to change script to use correct source data when changing sources.
Also, script currently also includes rotation measure correction to the rotation matrix (which only matters when generating using a polarized source such as 3C286).
 """
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np

from utils import misc

def calcGain(OnData,OffData,file_num,freq_len,src,beamwidth):
    def peval(data,p,f):
        d = data
        XG = p[0]
        YG = p[1]
        t = sp.zeros(len(data)*4)
        for i in range(0,len(t),4):
            t[i] = XG*d[i/4,0,f]
            t[i+1] = 0
            t[i+2] = 0
            t[i+3] =YG*d[(i+3)/4,3,f]
        return t 

    def residuals(p,errors, f,freq_val,src,theta,data,width):

        wavelength = 300.0/freq_val[f]
        BW = width[f]*180./sp.pi
        JtoK = (sp.pi*wavelength**2)/(8*1380.648*BW**2)
	Jsrc_name = ['3C286','3C48','3C67','3C147','3C295']
	Jsrc_val = [19.74748409*pow((750.0/freq_val[f]),0.49899785),
		    25.15445092*pow((750.0/freq_val[f]),0.75578842),
		    4.56303633*pow((750.0/freq_val[f]),0.59237327),
		    31.32846821*pow((750.0/freq_val[f]),0.52113534),
                    34.11187767*pow((750.0/freq_val[f]),0.62009421)]
	for i in range(0,len(Jsrc_name)):
	    if Jsrc_name[i]==src:
		src_ind = i
	PAsrc = [33.*sp.pi/180.,0.,0.,0.,0.,0.]
	Psrc = [0.07,0.,0.,0.,0.]
	Isrc = Jsrc_val[src_ind]*JtoK
	Qsrc = Isrc*Psrc[src_ind]*sp.cos(2*PAsrc[src_ind])
	Usrc = Isrc*Psrc[src_ind]*sp.sin(2*PAsrc[src_ind])
        Vsrc = 0
        XXsrc0 = Isrc-Qsrc
        YYsrc0 = Isrc+Qsrc
        source = sp.zeros(4*self.file_num)
        for i in range(0,len(source),4):
            source[i] = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc+0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
            source[i+1] = 0
            source[i+2] = 0
            source[i+3] = (0.5*(1-sp.cos(2*theta[i])*XXsrc0+sp.sin(2*theta[i])*Usrc+0.5*(1+sp.cos(2*theta[i])*YYsrc0)
        err = (source-self.peval(data,p,f))/errors
        return err
    
###################################################
# Setting labels for indices for later
    XX_ind = 0
    YY_ind = 3
    XY_ind = 1
    YX_ind = 2
 
    S_med_src = sp.zeros((file_num,4,freq_len))
    S_med = sp.zeros((file_num,4,freq_len))

    PA_on = []
    m=0
    for Data in OnBlocks:
        S_med_src[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med_src[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med_src[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med_src[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
	PA_on.append(ma.mean(Data.PA))
	Data.calc_freq()
	freq_val = Data.freq/1e6
        m+=1

    PA_off = []
    m=0
    for Data in OffBlocks:
        S_med[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
	PA_off.append(ma.mean(Data.PA))
        m+=1
 
    S_data = sp.zeros((file_num,4,freq_len))
    for i in range(0,len(S_med)):
        S_data[i,0,:] = S_med_src[i,0,:]-S_med[i,0,:]
        S_data[i,1,:] = S_med_src[i,1,:]-S_med[i,1,:]
	S_data[i,2,:] = S_med_src[i,2,:]-S_med[i,2,:]
	S_data[i,3,:] = S_med_src[i,3,:]-S_med[i,3,:]
#There are 2 parameters for this version p[0] is XX gain and p[1] is YY gain. 
    p0 = [1,1] # guessed preliminary values
    error = sp.ones(4*file_num)
    #Note that error can be used to weight the equations if not all set to one.

    p_val_out = sp.zeros((freq_len, 3))
    for f in range(0,freq_len):   
        plsq = leastsq(residuals,p0,args=(error,f,freq_val,src,PA_on,S_data,beamwidth),full_output=0, maxfev=5000)
        pval = plsq[0] # this is the 1-d array of results0

        p_val_out[f,0] = freq_val[f]
        p_val_out[f,1] = pval[0]
        p_val_out[f,2] = pval[1]

    out_path = output_root+sess+'_diff_gain_calc'+output_end
    np.savetxt(out_path,p_val_out,delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    DiffGainGen(str(sys.argv[1])).execute()

