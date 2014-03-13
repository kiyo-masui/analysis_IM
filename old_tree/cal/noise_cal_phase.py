"""Procedure to calculate the noise cal phase that causes some of the U signal to bleed into the V signal for the data. This expected to be at least session stable (if not more) as this phase is reset whenever guppi is started. Setting up the calcualtion of the phase using the point source data. 

Using the point source data, want to take the ratio of U/sqrt(U**2+V**2) as a function of frequency and cos(2*phi), where phi is the noise cal phase. 

Run in analysis_IM: python cal/noise_cal_phase.py input/tcv/noise_cal_gen.ini
Note that the .ini file should indicate which session(s) and sourse you want to use. Script is run using data from a single source. The output is saved in a file called point_source_differential_phase_compare.txt
 """
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np
from scipy import optimize
import pylab

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
    """Calculates a phase in frequency between the U and V data. 
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
#            session_nums[c] = file_middle.split('_')[0]
#            print session_nums[c]
            for Data in Len_set :
                freq_num = Data.dims[3] # Setting the frequency binning to match whatever it's been set to. 

            # Skipping any guppi data that doesn't have both an on and off scan
            if guppi_result == True : 
                if n_scans != 2 :
                    self.file_num -=1
            c+=1

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
            X_ind = 0
            Y_ind = 3
            U_ind = 1
            V_ind = 2
            
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
                        S_med_caloff_src[:,0] = ma.median(Data.data[:,X_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,1] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,2] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                        S_med_caloff_src[:,3] = ma.median(Data.data[:,Y_ind,off_ind,:],axis=0)

                        S_med_calon_src[:,0] = ma.median(Data.data[:,X_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,1] = ma.median(Data.data[:,U_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,2] = ma.median(Data.data[:,V_ind,on_ind,:],axis=0)
                        S_med_calon_src[:,3] = ma.median(Data.data[:,Y_ind,on_ind,:],axis=0)
                    
                    for Data in OffBlocks:
                        S_med_caloff[:,0] = ma.median(Data.data[:,X_ind,off_ind,:],axis=0)
                        S_med_caloff[:,1] = ma.median(Data.data[:,U_ind,off_ind,:],axis=0)
                        S_med_caloff[:,2] = ma.median(Data.data[:,V_ind,off_ind,:],axis=0)
                        S_med_caloff[:,3] = ma.median(Data.data[:,Y_ind,off_ind,:],axis=0)
 
                        S_med_calon[:,0] = ma.median(Data.data[:,X_ind,on_ind,:],axis=0)
                        S_med_calon[:,1] = ma.median(Data.data[:,U_ind,on_ind,:],axis=0)
                        S_med_calon[:,2] = ma.median(Data.data[:,V_ind,on_ind,:],axis=0)
                        S_med_calon[:,3] = ma.median(Data.data[:,Y_ind,on_ind,:],axis=0)
 
# We want Tcal as the source:
                    self.d[k,:] = S_med_calon_src[:,0]-S_med_caloff_src[:,0]
                    self.d[k+1,:] = S_med_calon_src[:,1]-S_med_caloff_src[:,1]
                    self.d[k+2,:] = S_med_calon_src[:,2]-S_med_caloff_src[:,2]
                    self.d[k+3,:] = S_med_calon_src[:,3]-S_med_caloff_src[:,3]
                    
#                    self.d[k,:] = 0.5*(S_med_calon_src[:,0]+S_med_caloff_src[:,0]-S_med_calon[:,0]-S_med_caloff[:,0])
#                    self.d[k+1,:] = 0.5*(S_med_calon_src[:,1]+S_med_caloff_src[:,1]-S_med_calon[:,1]-S_med_caloff[:,1])
#                    self.d[k+2,:] = 0.5*(S_med_calon_src[:,2]+S_med_caloff_src[:,2]-S_med_calon[:,2]-S_med_caloff[:,2])
#                    self.d[k+3,:] = 0.5*(S_med_calon_src[:,3]+S_med_caloff_src[:,3]-S_med_calon[:,3]-S_med_caloff[:,3])
                    k+=4

#        for a in range(0,4*self.file_num):
#            for b in range(0,freq_num):
#                print self.d[a,b]
#                if self.d[a,b] > 3 :
#                   self.d[a,b] = 3
#                elif self.d[a,b]<-3:
#                   self.d[a,b] = -3

#        fitfunc = lambda p,x: sp.arccos(sp.cos(p*x))
#        fitfunc = lambda p,x: sp.arccos(sp.cos(p[0]*x+p[1]))
        fitfunc = lambda p,x: sp.cos(p[0]*x +p[1])
#        fitfunc = lambda p,x:0.9*sp.cos(2*sp.pi/84*x+p[2]%(2*sp.pi))
        errfunc = lambda p, x, y: fitfunc(p,x)-y
        fitfunc2 = lambda q,x: q[0]*sp.sin(q[1]*x+q[2])
        errfunc2 = lambda q,x,y: fitfunc2(q,x)-y
        fitfunc3 = lambda s, x: sp.sin(2*(s[0]*x+s[1]%(2*sp.pi)))
        errfunc3 = lambda s,x,y: fitfunc3(s,x)-y
#        fitfunc3 = lambda s,x: sp.tan(2*sp.pi*s[0]+s[1])
#        errfunc3 = lambda s,x,y: sp.arctan(fitfunc3(s,x))-sp.arctan(y)

#        freqs = []
#        U_trunc = []
#        V_trunc = []
#        Mag_trunc = []
#        for i in range(245,150,-1):
#             freqs.append(freq_val[i])
#             U_trunc.append(self.d[2,i])
#             V_trunc.append(self.d[3,i])
#             Mag_trunc.append(sp.sqrt(self.d[2,i]**2+self.d[3,i]**2))
#        print freqs
#        print len(freqs)
#        print len(Mag_trunc)
#        print Mag_trunc
#        pylab.scatter(freqs,U_trunc,s=1,c='c',edgecolors='none')
#        pylab.scatter(freqs,V_trunc,s=1,c='b',edgecolors='none')
#        pylab.scatter(freqs,Mag_trunc,s=1,c='g',edgecolors='none')
#        m0 = [0.1,20.0,1.0]
#        M,success = optimize.leastsq(errfunc,m0[:],args=(freqs,Mag_trunc))
#       pylab.plot(freqs,fitfunc(M,freqs),c='r') 
#       print freq_val
        freqs = sp.zeros(len(freq_val))
        U_data = sp.zeros((self.file_num,len(freqs)))
        V_data = sp.zeros((self.file_num,len(freqs)))
        R_data = sp.zeros((self.file_num,len(freqs)))
        for i in range(0,len(freq_val)):
#            freqs[i] = freq_val[i]
            freqs[i] = freq_val[len(freq_val)-i-1]
            for j in range(0,self.file_num):
                U_data[j,i] = self.d[4*j+1,len(freq_val)-i-1]
                V_data[j,i] = self.d[4*j+2,len(freq_val)-i-1]
#               R_data[j,i] = sp.sin(2*sp.arctan(V_data[j,i]/U_data[j,i]))
#                R_data[j,i] = U_data[j,i]
#                R_data[j,i] = sp.arccos(U_data[j,i]/sp.sqrt(U_data[j,i]**2+ V_data[j,i]**2))
                R_data[j,i] = U_data[j,i]/sp.sqrt(U_data[j,i]**2+V_data[j,i]**2)

#        print R_data[:,0]                
#        print R_data[:,-1]
        print np.any(np.isnan(R_data))
        print np.any(np.isinf(R_data))
        R_fin = sp.zeros((2,self.file_num))
        print freqs
 
        for i in range(0,self.file_num):
            for j in range(0,len(freqs)):
                if int(freqs[j])==710:
                    mask_num = j
                if int(freqs[j])==740:
                    mask_num2 = j
#            mask_num = 7 #705.859 in frequency space
#            mask_num2 = 24 #719.140625 in frequency space
#            mask_num = 101 #775.390625 in frequency space
#            mask_num2 = 121 #790.234375 in frequency space
#            mask_num = 0
#            mask_num2 = -1

            Datain = R_data[i,mask_num:mask_num2]
            fin = freqs[mask_num:mask_num2]

            bad_pts = np.logical_or(np.isnan(Datain), np.isinf(Datain))
	    good_ind = np.where(np.logical_not(bad_pts))
            Datain = Datain[good_ind]
            fin = fin[good_ind]

# settings: 0.08, 0.10,0.12,0.16
            R0 = [0.18,1.0]
#            R0 = [0.08,]
            R,success = optimize.leastsq(errfunc,R0[:],args=(fin,Datain),maxfev=10000)
            R[1] = R[1]%(2*sp.pi)
            print success
#            print R
#            R = R%(2*sp.pi)
#            print R
#            title = 'Frequency is: %1.4f, Phase is: %1.4f' %(R[0],R[1])
#            pylab.suptitle(title)
            pylab.plot(freqs,fitfunc(R,freqs),label='Ratio_'+str(i),c='g')
#            R_avg = [0.1354,2.341]
#            pylab.plot(fin,fitfunc(R_avg,fin),label='test_'+str(i),c='b')
            pylab.scatter(freqs,R_data[i,:],label='Ratio_'+str(i),s=1,c='r',edgecolors='none')
            
            print 'Differential phase frequency for scan ',i,' is:',R[0]
            print 'Differential phase "phase" for scan ',i,' is:',R[1]
            R_fin[0,i] = R[0]
            R_fin[1,i] = R[1]
        
#        np.savetxt('3C286_diff_phase_slope.txt',R_fin[0],delimiter = ' ')
#        np.savetxt('3C286_diff_phase_intercept.txt',R_fin[1],delimiter = ' ')    
#        freq_avg = ma.mean(R_fin[0,:])
#        phase_avg = ma.mean(R_fin[1,:])
#        freq_avg = 0.1354
#        phase_avg = 2.341
#        print 'Average Differential phase frequency is:',freq_avg
#        print 'Average Differential phase "phase" is:',phase_avg
#        for i in range(0,self.file_num):        
#            pylab.scatter(freqs,U_data[i]*sp.cos(freq_avg*freqs+phase_avg)-V_data[i]*sp.sin(freq_avg*freqs+phase_avg),s=2,c='g',edgecolors='none',label='U_%d' %(i))
#            pylab.scatter(freqs,U_data[i]*sp.sin(freq_avg*freqs+phase_avg)+V_data[i]*sp.cos(freq_avg*freqs+phase_avg),s=2,c='b',edgecolor='none',label = 'V_%d' %(i))
#            pylab.scatter(freqs,U_data[i]*sp.cos(freq_avg*freqs+phase_avg)-V_data[i]*sp.sin(freq_avg*freqs+phase_avg),s=2,c='g',edgecolors='none',label='V_%1.2f_%1.2f_%d' %(R[0],R[1],i))
#            pylab.scatter(freqs,-U_data[i]*sp.sin(freq_avg*freqs+phase_avg)-V_data[i]*sp.cos(freq_avg*freqs+phase_avg),s=2,c='b',edgecolor='none',label = 'U_%1.2f_%1.2f_%d' %(R[0],R[1],i))

#        for i in range(0,4*self.file_num,4):
#            U0 = [1,0.05,0.1]
#            U,success = optimize.leastsq(errfunc,U0[:],args=(freqs,U_data[i,:]),maxfev=1000000)
#            U[2] = U[2]%(2*sp.pi)
#            V0 = [1,0.05,0.1]
#            V,success = optimize.leastsq(errfunc2,V0[:],args=(freqs,V_data[i,:]),maxfev=1000000)
#            V[2] = V[2]%(2*sp.pi)
#            pylab.plot(freqs,fitfunc(U,freqs),label='U_'+str(i),c='c')
#            pylab.plot(freqs,fitfunc2(V,freqs),label='V_'+str(i),c='b')
#            pylab.plot(freqs,fitfunc2(V,freqs)/fitfunc(U,freqs),label='Ratio_'+str(i),c='g')
#        i = 0
#            pylab.scatter(freqs,U_data[i,:],label='U_'+str(i),s=2,c='c',edgecolors='none')
#            pylab.scatter(freqs,V_data[i,:],label='V_'+str(i),s=2,c='y',edgecolors='none')
#            magnitude = sp.sqrt(self.d[i+3]**2+self.d[i+2]**2)
#            pylab.scatter(freqs,M_data[i,:], label = 'amplitude'+str(i),s=1,c='g', edgecolors='none')
#            m0 = [0.1,20.0,1.0]
#            M,success = optimize.leastsq(errfunc,m0[:],args=(freqs,M_data[i,:]))
#            pylab.plot(freqs,fitfunc(M,freqs),c='g')
#            pylab.scatter(freqs,U_data[i,:]/M_data[i,:],s=1,c='r',edgecolors='none')
#            test0 = [1.0,80.0,1.0]
#            UT,success = optimize.leastsq(errfunc,test0[:],args=(freqs,U_data[i,:]/M_data[i,:]))
#            UT[2] = UT[2]%(2*sp.pi)
#            pylab.plot(freqs,fitfunc(UT,freqs),c='r')
#            pylab.scatter(freqs,V_data[i,:]/M_data[i,:],s=1,c='m',edgecolors='none')
#            VT,success = optimize.leastsq(errfunc,test0[:],args=(freqs,V_data[i,:]/M_data[i,:]))
#            VT[2] = VT[2]%(2*sp.pi)
#            pylab.plot(freqs,fitfunc(VT,freqs),c='m')
#            pylab.scatter(freqs,(U_data[i,:]*V_data[i,:])/M_data[i,:]**2,s=1,c='y',edgecolors='none')
            
#            print 'U amplitude for scan ',i,' is:', U[0]
#            print 'U phase for scan ',i,' is:',U[2]
#            print 'U frequency for scan ',i,' is:',U[1]
#            print 'U offset for scan ',i,' is:',U[3]
#            print '-----------------------------------'
#            print 'V amplitude for scan ',i,' is:',V[0]
#            print 'V phase for scan ',i,' is:', V[2]
#            print 'V frequency for scan ',i,' is:',V[1]
#            print 'V offset for scan ',i,' is:',V[3]
#            print '-----------------------------------'
#            print 'M amplitude for scan ',i,' is:',M[0]
#            print 'M phase for scan ',i,' is:',M[2]
#            print 'M wavelength for scan ',i,' is:',M[1]
#            print 'M offset for scan ',i,' is:',M[3]
#            print '------------------------------------'
#            print 'U/M amplitude for scan ',i,' is:',UT[0]
#            print 'U/M phase for scan ',i,' is:',UT[2]
#            print 'U/M wavelength for scan ',i,' is:',UT[1]
#            print 'U/M offset for scan ',i,' is:',UT[3]
#            print '-----------------------------------'
#            print 'V/M amplitude for scan ',i,' is:',VT[0]
#            print 'V/M phase for scan ',i,' is:',VT[2]
#            print 'V/M wavelength for scan ',i,' is:',VT[1]
#            print 'V/M offset for scan ',i,' is:',VT[3]
#            print '==================================='
#            print 'The difference in phase for scan ',i,' between U and V is:',U[2]-V[2]
#            print 'The difference in phase for scan ',i,' between U/M and V/M is:',UT[2]-VT[2]
#            print '==================================='
        
#        leg = pylab.legend(fancybox='True')
#        leg.get_frame().set_alpha(0.25)    
        pylab.xlim(freqs[0],freqs[-1])
#        pylab.ylim(-1,1)
#        pylab.legend()
        pylab.savefig('point_source_differential_phase_compare.png')
        pylab.clf()
        
#        np.savetxt('mueller_params_calc.txt',p_val_out,delimiter = ' ')
#        np.savetxt('mueller_params_error.txt', p_err_out, delimiter = ' ')

#If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    MuellerGen(str(sys.argv[1])).execute()

