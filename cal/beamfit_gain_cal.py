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
from time_stream import rebin_time, combine_cal, stitch_windows_crude
from time_stream import split_bands
from map import pol_beam
from utils import misc
import cal.source
from cal import beam_fit
#from cal import flux_diff_gain_gen_beamfit

def calcGain(OnData,OffData,file_num,freq_val,src,beamwidth):
    """ Perform gain calibration on dataset. 
    Should be in a separate file, but I'm still figuring out how to do this. 
    """
    
    # Model for XX and YY Gain 
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

    #Setting up Residuals for Least Squares Fitting
    def residuals(p,errors,freq_val,src,theta,data,width,file_num):
        wavelength = 300./freq_val
        BW = width*sp.pi/180.
        JtoK = sp.pi*wavelength**2/(8*1380.648*BW**2)
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
        Isrc = Jsrc_val[src_ind]
        Qsrc = Isrc*Psrc[src_ind]*sp.cos(2*PAsrc[src_ind])
        Usrc = Isrc*Psrc[src_ind]*sp.sin(2*PAsrc[src_ind])
        Vsrc = 0
        XXsrc0 = (Isrc-Qsrc)*JtoK[0]
        YYsrc0 = (Isrc+Qsrc)*JtoK[1]
        expec =sp.zeros(4*file_num)
        for i in range(0,len(source),4):
            expec[i] = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc
                        +0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
            expec[i+1] = 0
            expec[i+2] = 0
            expec[i+3] = (0.5*(1-sp.cos(2*theta[i]))*XXsrc0+sp.sin(2*theta[i])*Usrc
                          +0.5*(1+sp.cos(2*theta[i]))*YYsrc0)
        err = (expec-peval(p,data))/errors
        return err

# Main Code
# Setting labels for indices for later
    XX_ind = 0
    YY_ind = 3
    XY_ind = 1
    YX_ind = 2
    freq_len = len(freq_val)
    S_med_src = sp.zeros((file_num,4,freq_len))
    S_med = sp.zeros((file_num,4,freq_len))

#Loading On/Off Data
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

#Setting actual data arrays.
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
        plsq = leastsq(residuals,p0,args=(error,freq_val[f],src,PA_on,S_data[:,:,f],
                       beamwidth[f],file_num),full_output=0, maxfev=5000)
        pval = plsq[0] # this is the 1-d array of results0
        p_val_out[f,0] = freq_val[f]
        p_val_out[f,1] = pval[0]
        p_val_out[f,2] = pval[1]


    wavelength = sp.ones((len(freq_val),2))
    for i in range(0,len(freq_val)):
        wavelength[i] = [300./freq_val[i],300./freq_val[i]]
    JtoK = sp.pi*(wavelength)**2/(8*1380.648*(beamwidth*sp.pi/180.)**2)
    return p_val_out,JtoK

##################################################################################
#Main Code:

#Data Info for Old Guppi Data:
#data_root = '/home/scratch/kmasui/converted_fits/GBT12A_418/'
#end = '.fits'
#source = '3C295'
#beam_cal_files = ['21_3C286_track_'+str(ii) for ii in range(18,26)]
#beam_cal_files = ['22_3C295_track_'+str(ii) for ii in range(59,67)]
#gain_cal_files = ['22_3C295_onoff_76-77','22_3C295_onoff_78-79']
#gain_cal_files = ['22_3C147_onoff_50-51','22_3C147_onoff_52-53','22_3C147_onoff_6-7','22_3C147_onoff_8-9']

#Data Info for Guppi Data:
data_root = '/home/scratch/tvoytek/converted_fits/GBT13A_510/'
end ='.fits'
#beam_cal_files = ['06_3C48_track_'+str(ii) for ii in range(26,34)]
#beam_cal_files = ['06_3C48_track_'+str(ii) for ii in range(42,50)]
beam_cal_files = ['06_3C48_track_'+str(ii) for ii in range(58,66)]
gain_cal_files = ['06_3C48_onoff_54-55','06_3C48_onoff_56-57']
source = '3C48'
out_dir = '/users/tvoytek/beamcal_results_guppi/'

#Data Info for Spectrometer Data:
#data_root = '/users/chanders/sdfits_files/July25/'
#end = '.raw.acs.fits'
#beam_cal_files = ['70:77/TGBT13A_510_06',]#Error because one time step is exactly at second 60.00
#beam_cal_files = ['82:89/TGBT13A_510_06',]#Noise Cal "off" always
#beam_cal_files = ['98:105/TGBT13A_510_06',]
#gain_cal_files = ['66:67/TGBT13A_510_06','68:69/TGBT13A_510_06']
#gain_cal_files = ['78:79/TGBT13A_510_06','80:81/TGBT13A_510_06']#Noise Cal "off" always
#gain_cal_files = ['94:95/TGBT13A_510_06','96:97/TGBT13A_510_06']
#source = '3C48'
#out_dir = '/users/tvoytek/beamcal_results/'

#Program Settings 
# Number of IFs corresponds to number of freq windows
# GUPPI=True means working with guppi data
# STITCH=True means multiple freq windows
# SPIDER=True only for spectrometer files (Guppi files stored separately).
# Beam_cal=True means running basic beam fit using spider scans
# Hermite=True means running polarized beam fit using hermite basis functions
# Gain_cal=True means running forward gain calibration using on-off scans
# Plotting means output plots produced for the data
# Scans is for the spectrometer spider scans (which are stored in one file).
IFs = tuple(np.arange(0,8))
STITCH=False
GUPPI=True
SPIDER=False
Scans=np.arange(0,8)
Beam_cal = True
Hermite = True
Gain_cal = True
Plotting = True

#File Prep/Preproccessing for Spider Scans
beam_cal_Blocks = []
for fname in beam_cal_files:
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    if SPIDER:
        if STITCH:
            for i in range(0,len(Scans)):
                Data = Reader.read(i,IFs)
                Stitch = stitch_windows_crude.stitch(Data)
                beam_cal_Blocks.append(Stitch)
        else:
            for i in range(0,len(Scans)):
                Data = Reader.read(i,0)
                beam_cal_Blocks.append(Data)
    elif STITCH:
        Data = Reader.read(0,IFs)
        Stitch = stitch_windows_crude.stitch(Data) 
        beam_cal_Blocks.append(Stitch)
    else:
        Data = Reader.read(0,0)
        trunc_Data=split_bands.split(Data,n_bands=1,bins_band=1026,offset=2559)
        beam_cal_Blocks.append(trunc_Data[0])

for Data in beam_cal_Blocks:
    if GUPPI:
        rotate_pol.rotate(Data, (-5, -7, -8, -6))
        cal_scale.scale_by_cal(Data, True, False, False, False, True,False)
        flag_data.flag_data(Data, 5, 0.1,40)
        rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (0., 1.), False, True) 
    else: 
        cal_scale.scale_by_cal(Data,True,False,False,False,True,True)
#        flag_data.flag_data(Data, 5, 0.1,1)
        rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (1., 0.), False, True) 
        #weights inverted due to cal inverted

beam_cal_Blocks[0].calc_freq()
beam_cal_freq = beam_cal_Blocks[0].freq

#File Prep/Processing for OnOff Scans
gain_cal_OnBlocks = []
gain_cal_OffBlocks = []
for fname in gain_cal_files:
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    if STITCH:
        OnData = Reader.read(0,IFs)
        OffData = Reader.read(1,IFs)
        OnStitch = stitch_windows_crude.stitch(OnData)
        OffStitch = stitch_windows_crude.stitch(OffData)
        gain_cal_OnBlocks.append(OnStitch)
        gain_cal_OffBlocks.append(OffStitch)
    else:
        OnData = Reader.read(0,0)
        OffData = Reader.read(1,0)
        trunc_OnData=split_bands.split(OnData,n_bands=1,bins_band=1026,offset=2559)
        trunc_OffData=split_bands.split(OffData,n_bands=1,bins_band=1026,offset=2559)
        gain_cal_OnBlocks.append(trunc_OnData[0])
        gain_cal_OffBlocks.append(trunc_OffData[0])

for Data in gain_cal_OnBlocks:
    if GUPPI:
        rotate_pol.rotate(Data, (-5, -7, -8, -6))
        cal_scale.scale_by_cal(Data, True, False, False, False, True,False)
        flag_data.flag_data(Data,5,0.1,40)
        rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (0., 1.), False, True)
    else:
        cal_scale.scale_by_cal(Data,True,False,False,False,True,True)
#        flag_data.flag_data(Data, 5, 0.1, 1)
        rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (1., 0.), False, True)

for Data in gain_cal_OffBlocks:
    if GUPPI:
        rotate_pol.rotate(Data, (-5, -7, -8, -6))
    	cal_scale.scale_by_cal(Data, True, False, False, False, True,False)
	flag_data.flag_data(Data, 5, 0.1,40)
	rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (0., 1.), False, True)
    else:
        cal_scale.scale_by_cal(Data,True,False,False,False,True,True)
#        flag_data.flag_data(Data, 5, 0.1, 1)
        rebin_freq.rebin(Data, 16, True, True)
        combine_cal.combine(Data, (1., 0.), False, True)

Data.calc_freq()
if GUPPI:
    gain_cal_freq = Data.freq/1e6
else: 
    gain_cal_freq = Data.freq/1e6-280

#Beam Cal Basic Fitting
if Beam_cal:
    BeamData = beam_fit.FormattedData(beam_cal_Blocks)

# Source object.  This just calculates the ephemeris of the source compared to
# where the telescope is pointing.
    S = cal.source.Source(source)

# Do a preliminary fit to just the XX and YY polarizations.  This is a
# non-linear fit to the Gaussian and gets things like the centriod and the
# Gaussian width.  All fits are channel-by-channel (independantly).
    center_offset, width, amps, Tsys = beam_fit.fit_simple_gaussian(BeamData, S)
    print width[0]
    elevation = ma.mean(BeamData.get_all_field('el'))

#Hermite Fitting
if Hermite:
# Basic basis functions to be used in the fit.
    HermiteBasis = pol_beam.HermiteBasis(beam_cal_freq, center_offset, width[:,0])
# Perform the fit.
    beam_params, scan_params, model_data = beam_fit.linear_fit(BeamData, HermiteBasis,S, 3, 2)

# Make a beam object from the basis funtions and the fit parameters (basis coefficients).
    Beam = pol_beam.LinearBeam(HermiteBasis, beam_params)

#Forward Gain Fitting
if Gain_cal:
    p_val_out,JtoK = calcGain(gain_cal_OnBlocks,gain_cal_OffBlocks,len(gain_cal_files),
                              gain_cal_freq,source,width)
    np.savetxt(out_dir+source+'_'+str(int(elevation))+'_test_cal.txt',p_val_out,delimiter = ' ')

#Plot Intermediate Results
if Plotting:
    if GUPPI:
        freq = beam_cal_freq/1e6
    else:
        freq = beam_cal_freq/1e6-280

#Pointing Offset Plot
    plt.figure()
    plt.plot(freq,center_offset[:,0])
    plt.plot(freq,center_offset[:,1])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Pointing Offset (Degrees)')
    plt.legend(('Az','El'))
    plt.title('Pointing Offset from '+source+' at %0.1f degrees El' %elevation)
    plt.xlim(700,900)
#    plt.ylim(-0.05,0.05)
    plt.grid()
    plt.savefig(out_dir+source+'_'+str(int(elevation))+'_pointing_offsets',dpi=300)
    plt.clf()

#Beam Width Plot
    plt.figure()
    plt.plot(freq,width[:,0])
    plt.plot(freq,width[:,1])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Beam Width (Degrees)')
    plt.legend(('XX','YY'))
    plt.xlim(700,900)
    plt.title('Symmetric Gaussian FWHM from '+source+' at %0.1f degrees El' %elevation)
    plt.grid()
    plt.savefig(out_dir+source+'_'+str(int(elevation))+'_beamwidth',dpi=300)
    plt.clf()

#Jansky to Kelvin Conversion Factor Plot
    if Gain_cal:
        plt.figure()
        plt.plot(freq,JtoK[:,0])
        plt.plot(freq,JtoK[:,1])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Jansky to Kelvin Conversion Factor')
        plt.xlim(700,900)
        plt.legend(('XX','YY'))
        plt.title('Jy to K Conversion Factor from '+source+' at %0.1f degrees El' %elevation)
        plt.grid()
        plt.savefig(out_dir+source+'_'+str(int(elevation))+'_JytoK_convert',dpi=300)
        plt.clf()

#Beam Map Plots
    if Hermite:
        n_chan = len(freq)
        beam_map = Beam.get_full_beam(100,1.)
        for i in range(0,20):
            f_ind = i*n_chan/20.
            plt.figure()
            pol_beam.plot_beam_map(beam_map[f_ind,...],color_map=1,side=1.,
                                   normalize='max03',rotate='XXYYtoIQ')
            cbar = plt.colorbar()
            cbar.set_label(r"Intensity (Normalized to I beam center with Sign)")
            plt.xlabel(r"IQUV, Azimuth (degrees)")
            plt.ylabel(r"Elevation (degrees)")
            plt.title("Beam Patterns for %0.1f MHz at %0.1f degrees El" %(freq[f_ind],elevation))
            raw_title = out_dir+source+'_'+str(int(elevation))+'_'+str(int(freq[f_ind]))+'_MHz_beam_pattern'
            plt.savefig(raw_title,dpi=300)
            plt.clf()

#Tsys Plot
    if Beam_cal:
        plt.figure()
        Tsys_Kel = Tsys
        if Gain_cal:
            Tsys_Kel[:,0] = Tsys[:,0]*p_val_out[:,1]
            Tsys_Kel[:,1] = Tsys[:,1]*p_val_out[:,2]
        plt.plot(freq,Tsys_Kel[:,0])
        plt.plot(freq,Tsys_Kel[:,1])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Tsys')
        if Gain_cal:
            plt.ylabel('Tsys (Kelvin)')
        plt.grid()
        plt.xlim(700,900)
        plt.title('Tsys from '+source+' at %0.1f degrees El' %elevation)
        plt.legend(('XX','YY'))
        if Gain_cal:
            plt.savefig(out_dir+source+'_'+str(int(elevation))+'_Tsys_Kelvin',dpi=300)
        else:
            plt.savefig(out_dir+source+'_'+str(int(elevation))+'_Tsys',dpi=300)
        plt.clf()

#Tcal Plot
    if Gain_cal:
        plt.figure()
        plt.plot(freq,p_val_out[:,1])
        plt.plot(freq,p_val_out[:,2])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Tcal (Kelvin)')
        plt.grid()
        plt.legend(('XX','YY'))
        plt.xlim(700,900)
        plt.title('Tcal from '+source+' at %0.1f degrees El' %elevation)
        plt.savefig(out_dir+source+'_'+str(int(elevation))+'_Tcal_Kelvin',dpi=300)
        plt.clf()

