"""Procedure to test the effectiveness of the flux and differential gain calibration.


Run in analysis_IM: python cal/flux_diff_gain_cal_test.py input/tcv/diff_gain_tes_guppi.ini
Note that the .ini file should indicate which session(s) and sourse you want to use. Script is run using data from a single source.
 """
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np
import pylab as pl

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
prefix = 'fgc_'


class MuellerGen(object) :
    """Calculates uncalibrated mueller matrix elements for each frequency using on-off scans at many separate parallactic angles. 
    """
   
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
# Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                          prefix=prefix)
        self.feedback = feedback

    def peval(self,p,f):
        d = self.d
        XG = p[0]
        YG = p[1]
        theta = self.theta
        t = self.function
        
        for i in range(0,len(t),4):
            t[i] = XG*d[i,f]
            t[i+1] = 0
            t[i+2] = 0
            t[i+3] =YG*d[i+3,f]
        return t 

    def residuals(self, p,errors, f,freq_val):
        theta = self.theta
#        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)*2 
#        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)*(2.28315426-0.000484307905*freq_val[f]) # Added linear fit for Jansky to Kelvin conversion.
#        Isrc = 19.74748409*pow((750.0/freq_val[f]),0.49899785)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for 3C286
        Isrc = 25.15445092*pow((750.0/freq_val[f]),0.75578842)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for  3C48
#        Isrc = 4.56303633*pow((750.0/freq_val[f]),0.59237327)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for 3C67
        PAsrc = 33.0*sp.pi/180.0 # for 3C286, doesn't matter for unpolarized. 
#        Psrc = 0.07 #for 3C286 
        Psrc = 0 #for #3C48,3C67 
        Qsrc = Isrc*Psrc*sp.cos(2*PAsrc) 
        Usrc = Isrc*Psrc*sp.sin(2*PAsrc) 
        Vsrc = 0
        XXsrc0 = Isrc-Qsrc
        YYsrc0 = Isrc+Qsrc
#        XXsrc = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc+0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
#        YYsrc = (0.5*(1-sp.cos(2*theta[i]))*XXsrc0+sp.sin(2*theta[i])*Usrc+0.5*(1+sp.cos(2*theta[i]))*YYsrc0)
        source = sp.zeros(4*self.file_num)
        for i in range(0,len(source),4):
            source[i] = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc+0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
            source[i+1] = 0
            source[i+2] = 0
            source[i+3] = (0.5*(1-sp.cos(2*theta[i]))*XXsrc0+sp.sin(2*theta[i])*Usrc+0.5*(1+sp.cos(2*theta[i]))*YYsrc0)
        err = (source-self.peval(p,f))/errors
        return err
    
    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        guppi_result = params['Guppi_test']
        output_root = params['output_root']
        output_end = params['output_end']
        file_name = params['file_middles'][0].split('/')[1]
#        print file_name
        sess = file_name.split('_')[0]
#        print sess

        self.file_num = len(params['file_middles']) # getting a variable for number of calibrator files being used

# Need to remove count for calibrator files that are not the right size.
        session_nums = sp.zeros(self.file_num)
        c = 0
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle + 
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname)
            n_scans = len(Reader.scan_set)
            Len_set = Reader.read(0,0,force_tuple=True)
#            session_nums[c] = file_middle.split('_')[0]
#            print session_nums[c]
            for Data in Len_set :
                freq_num = Data.dims[3] # Setting the frequency binning to match whatever it's been set to. 
            if guppi_result == True : 
                if n_scans != 2 :
                    self.file_num -=1
            elif guppi_result == False :
                if n_scans != 4 :
                    self.file_num -=1
            c+=1

# Need to know the general frequency binning (going to assume that it's 200 for guppi, 260 for spectrometer, aka 1 MHz binning)
#        if guppi_result == True :
#            freq_num = 200 
        if guppi_result == False :
#            freq_num = 260    
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
            XX_ind = 0
            YY_ind = 3
            XY_ind = 1
            YX_ind = 2
            
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
            
#Building the measured data into arrays (guppi version)         
            if guppi_result == True : 
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
                        S_med_caloff_src[:,0] = ma.median(Data.data[:,XX_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,1] = ma.median(Data.data[:,XY_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,2] = ma.median(Data.data[:,YX_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,3] = ma.median(Data.data[:,YY_ind,off_ind,:],axis=0)

                        S_med_calon_src[:,0] = ma.median(Data.data[:,XX_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,1] = ma.median(Data.data[:,XY_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,2] = ma.median(Data.data[:,YX_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,3] = ma.median(Data.data[:,YY_ind,on_ind,:],axis=0)
                    
                    for Data in OffBlocks:
                        S_med_caloff[:,0] = ma.median(Data.data[:,XX_ind,off_ind,:],axis=0)
                        S_med_caloff[:,1] = ma.median(Data.data[:,XY_ind,off_ind,:],axis=0)
                        S_med_caloff[:,2] = ma.median(Data.data[:,YX_ind,off_ind,:],axis=0)
                        S_med_caloff[:,3] = ma.median(Data.data[:,YY_ind,off_ind,:],axis=0)
 
                        S_med_calon[:,0] = ma.median(Data.data[:,XX_ind,on_ind,:],axis=0)
                        S_med_calon[:,1] = ma.median(Data.data[:,XY_ind,on_ind,:],axis=0)
                        S_med_calon[:,2] = ma.median(Data.data[:,YX_ind,on_ind,:],axis=0)
                        S_med_calon[:,3] = ma.median(Data.data[:,YY_ind,on_ind,:],axis=0)
 
                     
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

        print self.d[:,150]

# self is a directory of Tsrc data. I can use this as my corrected data for plotting. 
# It looks like at the end k should equal the number of on, off src sets (so if one set of onoff scans eg 6-9, k = 8)
       

# source data for use in ploting 
        XXsrc_3C286 = sp.zeros(freq_len)
        YYsrc_3C286 = sp.zeros(freq_len)
        XXsrc_3C48 = sp.zeros(freq_len)
        YYsrc_3C48 = sp.zeros(freq_len)
        Usrc_3C286 = sp.zeros(freq_len)
       
        for f in range(0,freq_num):
            Isrc_3C286 = 19.74748409*pow((750.0/freq_val[f]),0.49899785)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for 3C286
            Isrc_3C48 = 25.15445092*pow((750.0/freq_val[f]),0.75578842)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for  3C48 
#        Isrc_3C67 = 4.56303633*pow((750.0/freq_val[f]),0.59237327)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for 3C67
            PAsrc_3C286 = 33.0*sp.pi/180.0 # for 3C286, doesn't matter for unpolarized. 
            Psrc_3C286 = 0.07 #for 3C286 
            Psrc = 0 #for #3C48,3C67  
#        Qsrc = Isrc*Psrc*sp.cos(2*PAsrc) 
#        Usrc = Isrc*Psrc*sp.sin(2*PAsrc) 
#        Vsrc = 0
            XXsrc_3C286[f] = Isrc_3C286*(1-Psrc_3C286*sp.cos(2*PAsrc_3C286))
            YYsrc_3C286[f] = Isrc_3C286*(1+Psrc_3C286*sp.cos(2*PAsrc_3C286))
            Usrc_3C286[f] = Isrc_3C286*Psrc_3C286*sp.sin(2*PAsrc_3C286)
            XXsrc_3C48[f] = Isrc_3C48
            YYsrc_3C48[f] = Isrc_3C48
#            Usrc_3C48 = 0

        XX_compare = sp.zeros((freq_len,k/4))
        YY_compare = sp.zeros((freq_len,k/4))
        XX_PA_3C286 = sp.zeros((freq_len,k/4))
        YY_PA_3C286 = sp.zeros((freq_len,k/4))
        for c in range(0,k/4):
            XX_PA_3C286[:,c] = 0.5*(1+sp.cos(2*self.theta[4*c]))*XXsrc_3C286[:]-sp.sin(2*self.theta[4*c])*Usrc_3C286[:]+0.5*(1-sp.cos(2*self.theta[4*c]))*YYsrc_3C286[:]
            YY_PA_3C286[:,c] = 0.5*(1-sp.cos(2*self.theta[4*c]))*XXsrc_3C286[:]+sp.sin(2*self.theta[4*c])*Usrc_3C286[:]+0.5*(1+sp.cos(2*self.theta[4*c]))*YYsrc_3C286[:]
            XX_compare[:,c] = self.d[c*4,:]
            YY_compare[:,c] = self.d[c*4+3,:] 

#        pl.plot(freq_val,XXsrc_3C286,label='XX_3C286',color='b')
#        pl.plot(freq_val,YYsrc_3C286,label='YY_3C286')
        pl.plot(freq_val,XXsrc_3C48,label='XX_3C48',color='b')
      
        for d in range(0,k/4):
            if d == 0:
                col = 'g'
            elif d == 1:
                col = 'r'
            elif d == 2:
                col = 'c'
            elif d == 3:
                col = 'm'
            elif d == 4:
                col = 'y'
            else:
                col = 'k'
            pl.plot(freq_val,XX_compare[:,d],label='XX_'+str(d),color = col)
#            pl.plot(freq_val,YY_compare[:,d],label='YY_'+str(d),color = col)
#            pl.plot(freq_val,XX_PA_3C286[:,d],label='XXsrc_'+str(d), color = col)
#            pl.plot(freq_val,YY_PA-3C286[:,d],label='YYsrc_'+str(d), color = col)
        leg = pl.legend(fancybox='True')
        leg.get_frame().set_alpha(0.25)
        pl.ylim(20,55)
        pl.xlabel("Frequency (MHz)")
        pl.ylabel("Temperature (K)")
        title0 = sess+'_init_XX_3C48_comparison_Test.png'
        pl.savefig(title0)
        pl.clf()

        #There are 2 parameters for this version p[0] is XX gain and p[1] is YY gain. 
#        p0 = [1,1] # guessed preliminary values
#        error = sp.ones(4*self.file_num)
        #Note that error can be used to weight the equations if not all set to one.

#        p_val_out = sp.zeros((freq_len, 3))
 #       p_err_out = sp.zeros((freq_len, 17))
     
#        for f in range(0,freq_len):   
#            plsq = leastsq(self.residuals,p0,args=(error,f,freq_val),full_output=0, maxfev=5000)
#            pval = plsq[0] # this is the 1-d array of results0

#            p_val_out[f,0] = freq_val[f]
#            p_val_out[f,1] = pval[0]
#            p_val_out[f,2] = pval[1]

#        sess_num = int(session_nums[0])
#        print sess_num
#        np.savetxt(output_root+str(sess_num)+'_flux_mueller_matrix_calc'+output_end, p_val_out, delimiter = ' ')
#        out_path = output_root+sess+'_diff_gain_calc'+output_end
#        np.savetxt(out_path,p_val_out,delimiter = ' ')
#        out_path = sess+'_corrected_data_' + output_end
#        np.savetxt(out_path, self.d, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

