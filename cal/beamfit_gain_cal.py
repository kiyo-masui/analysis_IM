import os

import numpy as np
import numpy.ma as ma
import scipy as sp
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import *
import ephem
from numpy import ma
import matplotlib.animation as animation

from core import fitsGBT, dir_data
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal
from map import pol_beam
from utils import misc
import cal.source
from cal import beam_fit
#from cal import flux_diff_gain_gen_beamfit

def calcGain(OnData,OffData,file_num,freq_len,src,beamwidth):
    """ Perform gain calibration on dataset.

    """
    def peval(p,data):
        d = data
        XG = p[0]
        YG = p[1]
        act = sp.zeros(len(data)*4)
        for i in range(0,len(act),4):
            act[i] = XG*d[i/4,0]
            act[i+1] = 0
            act[i+2] = 0
            act[i+3] =YG*d[(i+3)/4,3]
        return act

    def residuals(p,errors,freq_val,src,theta,data,width,file_num):

        wavelength = 300.0/freq_val
        BW = width*sp.pi/180.
        JtoK = (sp.pi*wavelength**2)/(8*1380.648*BW**2)
        Jsrc_name = ['3C286','3C48','3C67','3C147','3C295']
        Jsrc_val = [19.74748409*pow((750.0/freq_val),0.49899785),
                    25.15445092*pow((750.0/freq_val),0.75578842),
                    4.56303633*pow((750.0/freq_val),0.59237327),
                    31.32846821*pow((750.0/freq_val),0.52113534),
                    34.11187767*pow((750.0/freq_val),0.62009421)]
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
        expec =sp.zeros(4*file_num)
        for i in range(0,len(source),4):
            expec[i] = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc+0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
            expec[i+1] = 0
            expec[i+2] = 0
            expec[i+3] = (0.5*(1-sp.cos(2*theta[i]))*XXsrc0+sp.sin(2*theta[i])*Usrc+0.5*(1+sp.cos(2*theta[i]))*YYsrc0)
        err = (expec-peval(p,data))/errors
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
    for Data in OnData:
        S_med_src[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med_src[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med_src[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med_src[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
        for i in range(0,4):
            PA_on.append(ma.mean(Data.PA))
        Data.calc_freq()
        freq_val = Data.freq/1e6
        m+=1

    PA_off = []
    m=0
    for Data in OffData:
        S_med[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
        for i in range(0,4):
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
        plsq = leastsq(residuals,p0,args=(error,freq_val[f],src,PA_on,S_data[:,:,f],beamwidth[f],file_num),full_output=0, maxfev=5000)
        pval = plsq[0] # this is the 1-d array of results0

        p_val_out[f,0] = freq_val[f]
        p_val_out[f,1] = pval[0]
        p_val_out[f,2] = pval[1]

    JtoK = sp.pi*(300./freq_val)**2/(8*1380.648*(beamwidth*sp.pi/180.)**2)

#    out_path = output_root+sess+'_diff_gain_calc'+output_end
#    np.savetxt(out_path,p_val_out,delimiter = ' ')
    return p_val_out,JtoK

data_root = '/home/scratch/kmasui/converted_fits/GBT12A_418/'
end = '.fits'

source = '3C295'
gain_source = '3C295'

#beam_cal_files = ['21_3C286_track_'+str(ii) for ii in range(18,26)]
beam_cal_files = ['22_3C295_track_'+str(ii) for ii in range(59,67)]

gain_cal_files = ['22_3C295_onoff_76-77','22_3C295_onoff_78-79']
#gain_cal_files = ['22_3C147_onoff_50-51','22_3C147_onoff_52-53','22_3C147_onoff_6-7','22_3C147_onoff_8-9']

beam_cal_Blocks = []
for fname in beam_cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    Data = Reader.read(0,0)
    beam_cal_Blocks.append(Data)

for Data in beam_cal_Blocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)
    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

gain_cal_OnBlocks = []
gain_cal_OffBlocks = []
for fname in gain_cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    OnData = Reader.read(0,0)
    OffData = Reader.read(1,0)
    gain_cal_OnBlocks.append(OnData)
    gain_cal_OffBlocks.append(OffData)

for Data in gain_cal_OnBlocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)
    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

for Data in gain_cal_OffBlocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)

    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

Data.calc_freq()

BeamData = beam_fit.FormattedData(beam_cal_Blocks)
#GainCalOnData = beam_fit.FormattedData(gain_cal_OnBlocks)
#GainCalOffData = beam_fit.FormattedData(gain_cal_OffBlocks)

# Source object.  This just calculates the ephemeris of the source compared to
# where the telescope is pointing.
S = cal.source.Source(source)

# Do a preliminary fit to just the XX and YY polarizations.  This is a
# non-linear fit to the Gaussian and gets things like the centriod and the
# Gaussian width.  All fits are channel-by-channel (independantly).
center_offset, width, amps, Tsys = beam_fit.fit_simple_gaussian(BeamData, S)

# Basis basis functions to be used in the fit.
HermiteBasis = pol_beam.HermiteBasis(Data.freq, center_offset, width)
# Perform the fit.
beam_params, scan_params, model_data = beam_fit.linear_fit(BeamData, HermiteBasis,
                                                           S, 3, 2)
p_val_out,JtoK = calcGain(gain_cal_OnBlocks,gain_cal_OffBlocks,len(gain_cal_files),len(Data.freq),gain_source,width)

#out_path = output_root+sess+'_diff_gain_calc'+output_end
np.savetxt(source+'_test_cal',p_val_out,delimiter = ' ')

# Make a beam object from the basis funtions and the fit parameters (basis
# coefficients).
Beam = pol_beam.LinearBeam(HermiteBasis, beam_params)

# Some plots.
beam_map = Beam.get_full_beam(100, 1.)

freq = Data.freq/1e6
plt.figure()
plt.plot(freq,center_offset[:,0])
plt.plot(freq,center_offset[:,1])
plt.xlabel('Frequency (MHz)')
plt.ylabel('Pointing Offset (Degrees)')
plt.legend(('Az','El'))
plt.xlim(700,900)
plt.ylim(-0.05,0.05)
plt.grid()
plt.savefig(source+'_pointing_offsets',dpi=300)
plt.clf()

plt.figure()
plt.plot(freq,width)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Beam Width (Degrees)')
plt.xlim(700,900)
plt.grid()
plt.savefig(source+'_beamwidth',dpi=300)
plt.clf()

plt.figure()
plt.plot(freq,JtoK)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Jansky to Kelvin Conversion')
plt.xlim(700,900)
plt.grid()
plt.savefig(source+'_JytoK_convert',dpi=300)
plt.clf()

n_chan = len(freq)
for i in range(0,20):
    f_ind = i*n_chan/20.
    plt.figure()
#    color_map=1
#    spec_beam_map = beam_map[f_ind,...]
#    n_pol = spec_beam_map.shape[0]
#    n_side = spec_beam_map.shape[1]
#    norm_x = spec_beam_map[0,n_side//2,n_side//2]
#    norm_y = spec_beam_map[3,n_side//2,n_side//2]
#    print np.where(np.isnan(norm_x))[0],np.where(np.isinf(norm_x))[0],np.where(norm_x==0)[0]
#    print np.where(np.isnan(norm_y))[0],np.where(np.isinf(norm_y))[0],np.where(norm_y==0)[0]
#    badx = 0
#    bady = 0
#    if norm_x==0:
#        if norm_y==0:
#            norm_x =1.
#        else:
#            norm_x = norm_y 
#        print 'XX Normalization failed'
#        badx = 1
#    if norm_y==0:
#        if norm_x==0:
#            norm_y = 1.
#        else:
#            norm_y = norm_x
#        bady = 1
#        print 'YY Normalization failed'
#    norm_cross = np.sqrt(norm_x*norm_y)
#    spec_beam_map[0] /=norm_x
#    spec_beam_map[3] /=norm_y
#    spec_beam_map[[1,2]]/=norm_cross
#    neg = spec_beam_map<0
#    spec_beam_map = abs(spec_beam_map)**color_map
#    spec_beam_map.shape = (n_pol*n_side,n_side)
#    spec_beam_map.shape = (n_pol*1.,1.)
#    plt.imshow(spec_beam_map.T)
#    pol_beam.plot_beam_map(beam_map[f_ind,...],color_map=0.5)
#    cbar = plt.colorbar()
#    if badx==0:
#        if bady==0:           
#           cbar.set_label(r"Beam Center Normalized Intensity")
#        elif bady==0:
#           cbar.set_label(r"Beam Center Normalized Intensity (Y-norm failed)")
#    elif badx==1:
#        if bady==0:
#           cbar.set_label(r"Beam Center Normalized Intensity (X-norm failed)")
#        elif bady==1:
#           cbar.set_label(r"Beam Center Un-Normalized Intensity")
#    plt.xlabel(r"XX'XY'YX'YY, Azimuth")
#    plt.ylabel(r"Elevation")
#    plt.xlabel(r"XX'XY'YX'YY, Azimuth (degrees)")
#    plt.ylabel(r"Elevation (degrees)")
    pol_beam.plot_beam_map(beam_map[f_ind,...],color_map=0.5,side=1.,normalize='max03',rotate='XXYYtoIQ')
    cbar = plt.colorbar()
    cbar.set_label(r"Root Intensity (Normalized to I beam center with Sign)")
    plt.xlabel(r"IQUV, Azimuth (degrees)")
    plt.ylabel(r"Elevation (degrees)")
    plt.title("Beam Patterns for %0.1f MHz" %freq[f_ind])
    raw_title = source+'_'+str(int(freq[f_ind]))+'_MHz_beam_pattern'
    plt.savefig(raw_title,dpi=300)
    plt.clf()

plt.figure()
Tsys_Kel = Tsys
Tsys_Kel[:,0] = Tsys[:,0]*p_val_out[:,1]
Tsys_Kel[:,1] = Tsys[:,1]*p_val_out[:,2]
plt.plot(freq,Tsys_Kel[:,0])
plt.plot(freq,Tsys_Kel[:,1])
plt.xlabel('Frequency (MHz)')
plt.ylabel('Tsys (Kelvin)')
plt.grid()
plt.legend(('XX','YY'))
plt.savefig(gain_source+'_Tsys_Kelvin',dpi=300)
plt.clf()

plt.figure()
plt.plot(freq,p_val_out[:,1])
plt.plot(freq,p_val_out[:,2])
plt.xlabel('Frequency (MHz)')
plt.ylabel('Tcal (Kelvin)')
plt.grid()
plt.legend(('XX','YY'))
plt.savefig(gain_source+'_Tcal_Kelvin',dpi=300)
plt.clf()


#plt.figure()
#this_data, this_weight = BeamData.get_data_weight_chan(35)
#plt.plot(this_data[0,:])
#plt.plot(model_data[35,0,:])

#pol_beam.plot_beam_map(beam_map[35,...])






