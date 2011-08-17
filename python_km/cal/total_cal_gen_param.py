"""Procedure to calculate the Mueller parameters for each frequency from on-off scans of a calibrator such as 3C286."""
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
        d = self.d
        dG = p[0]
        al = ((p[1]+180)%360-180)*sp.pi/180
        ps = ((p[2]+180)%360-180)*sp.pi/180
        ph = ((p[3]+180)%360-180)*sp.pi/180
        ep = p[4]
        ch = ((p[5]+180)%360-180)*sp.pi/180
        flux = p[6] #This is a variable for accounting for the difference in flux.
        be = ((p[7]+180)%360-180)*sp.pi/180
        theta = self.theta
        t = self.function
        m_tot = sp.zeros((4,4),float)
        m_tot[0,0] = flux
        m_IQ = flux*(0.5*dG*sp.cos(2.*al)-2.*ep*sp.cos(ph-ch)*sp.sin(2.*al))
        m_IU = flux*(0.5*dG*sp.sin(2.*al)*sp.cos(ch)+2*ep*(sp.cos(al)*sp.cos(al)*sp.cos(ph)-sp.sin(al)*sp.sin(al)*sp.cos(ph-2.*ch)))
        m_tot[0,1] = m_IQ
        m_tot[0,2] = m_IU
        m_tot[0,3] = flux*(0.5*dG*sp.sin(2.*al)*sp.sin(ch)+2*ep*(sp.cos(al)*sp.cos(al)*sp.sin(ph)+sp.sin(al)*sp.sin(al)*sp.sin(ph-2.*ch)))
        m_tot[1,0] = flux*0.5*dG
        m_QQ = flux*sp.cos(2.*al)
        m_QU = flux*sp.sin(2*al)*sp.cos(ch)
        m_tot[1,1] = m_QQ
        m_tot[1,2] = m_QU
        m_tot[1,3] = flux*sp.sin(2.*al)*sp.sin(ch)
        m_tot[2,0] = flux*2*ep*sp.cos(ph+ps)
        m_UQ = -flux*sp.sin(2*al)*sp.cos(ps+ch)
        m_UU = flux*(sp.cos(al)*sp.cos(al)*sp.cos(ps)-sp.sin(al)*sp.sin(al)*sp.cos(ps+2.*ch))
        m_tot[2,1] = m_UQ
        m_tot[2,2] = m_UU
        m_tot[2,3] = flux*(-sp.cos(al)*sp.cos(al)*sp.sin(ps)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2.*ch))
        m_tot[3,0] = flux*2*ep*sp.sin(ps+ph)
        m_VQ = -flux*sp.sin(2*al)*sp.sin(ps+ch)
        m_VU = flux*(sp.sin(ps)*sp.cos(al)*sp.cos(al)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2.*ch))
        m_tot[3,1] = m_VQ
        m_tot[3,2] = m_VU
        m_tot[3,3] = flux*(sp.cos(al)*sp.cos(al)*sp.cos(ps)+sp.sin(al)*sp.sin(al)*sp.cos(ps+2.*ch))
        M_total = sp.mat(m_tot)
        M_total = M_total.I

        for i in range(0,len(t),4):
            t[i] = M_total[0,0]*d[i,f]+M_total[0,1]*d[i+1,f]+M_total[0,2]*d[i+2,f]+M_total[0,3]*d[i+3,f]
            t[i+1] =(M_total[1,0]*sp.cos(2*theta[i])-M_total[2,0]*sp.sin(2*theta[i]))*d[i,f]+(M_total[1,1]*sp.cos(2*theta[i])-M_total[2,1]*sp.sin(2*theta[i]))*d[i+1,f]+(M_total[1,2]*sp.cos(2*theta[i])-M_total[2,2]*sp.sin(2*theta[i]))*d[i+2,f]+(M_total[1,3]*sp.cos(2*theta[i])-M_total[2,3]*sp.sin(2*theta[i]))*d[i+3,f]
            t[i+2] =(M_total[1,0]*sp.sin(2*theta[i])+M_total[2,0]*sp.cos(2*theta[i]))*d[i,f]+(M_total[1,1]*sp.sin(2*theta[i])+M_total[2,1]*sp.cos(2*theta[i]))*d[i+1,f]+(M_total[1,2]*sp.sin(2*theta[i])+M_total[2,2]*sp.cos(2*theta[i]))*d[i+2,f]+(M_total[1,3]*sp.sin(2*theta[i])+M_total[2,3]*sp.cos(2*theta[i]))*d[i+3,f]
            t[i+3] =M_total[3,0]*d[i,f]+M_total[3,1]*d[i+1,f]+M_total[3,2]*d[i+2,f]+M_total[3,3]*d[i+3,f]
        return t 

    def residuals(self, p,errors, f,freq_val):
#        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)*2 
#        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)*(2.28315426-0.000484307905*freq_val[f]) # Added linear fit for Jansky to Kelvin conversion.
        Isrc = 19.74748409*pow((750.0/freq_val[f]),0.49899785)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for 3C286
#        Isrc = 19.6*pow((750.0/freq_val[f]),0.495)*(2.28315426-0.000484307905*freq_val[f]) # My fit solution for  3C48
        PAsrc = 33.0*sp.pi/180.0 # for 3C286 
        Psrc = 0.07 #for 3C286 
#        Psrc = 0 #for #3C48 
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
        output_root = params['output_root']
        output_end = params['output_end']

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
            session_nums[c] = file_middle.split('_')[0]
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

        #The seven parameters are in order deltaG[0], alpha[1], psi[2], phi[3], epsilon[4], chi[5], flux[6]=> the parameter vector is p
        p0 = [0.3,90.0,170.0,10.0,0.016,0.00,2.0,0.0] # guessed preliminary values
        error = sp.ones(4*self.file_num)
        #Note that error can be used to weight the equations if not all set to one.

        p_val_out = sp.zeros((freq_len, 9))
 #       p_err_out = sp.zeros((freq_len, 17))
     
        for f in range(0,freq_len):   
            plsq = leastsq(self.residuals,p0,args=(error,f,freq_val),full_output=0, maxfev=5000)
            pval = plsq[0] # this is the 1-d array of results
#            perr = plsq[1] # this is a 2d array representing the estimated covariance of the results. - Not working properly. 

#            Mueller = sp.mat([[pval[0],pval[1],pval[2],pval[3]],[pval[4],pval[5],pval[6],pval[7]],[pval[8],pval[9],pval[10],pval[11]],[pval[12],pval[13],pval[14],pval[15]]])
#            Mueller = Mueller.I  
#            print Mueller
            pval[1]=(pval[1]+180)%360-180
            pval[2]=(pval[2]+180)%360-180
            pval[3]=(pval[3]+180)%360-180
            pval[5]=(pval[5]+180)%360-180
#            pval[7]=(pval[7]+180)%360-180

            p_val_out[f,0] = freq_val[f]
            p_val_out[f,1] = pval[0]
            p_val_out[f,2] = pval[1]
            p_val_out[f,3] = pval[2]
            p_val_out[f,4] = pval[3]
            p_val_out[f,5] = pval[4]
            p_val_out[f,6] = pval[5]
            p_val_out[f,7] = pval[6]
#            p_val_out[f,8] = pval[7]

#        sess_num = int(session_nums[0])
#        print sess_num
#        np.savetxt(output_root+str(sess_num)+'_flux_mueller_matrix_calc'+output_end, p_val_out, delimiter = ' ')
        np.savetxt('mueller_params_calc.txt',p_val_out,delimiter = ' ')
#        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

