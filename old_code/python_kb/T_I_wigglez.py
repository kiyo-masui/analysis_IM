import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *
import gbt_cal

cal_filename = '/home/scratch/kbandura/GBT10B_036/02_3c48_onoff_161-164.raw.acs.fits'

T_noise,T_sys =  gbt_cal.getcal_3c48(cal_filename)

print T_noise.shape
print T_sys.shape

print "finished loading cal"

day = "d2"

filenames = ['/home/scratch/kbandura/GBT08B_037/aug30/d2_field4_drift_67-72.raw.acs.fits','/home/scratch/kbandura/GBT08B_037/aug30/d2_field4_drift_61-66.raw.acs.fits']
for filename in filenames:
 spectra_clean,ra_gbt,dec_gbt,scan_n,crpix1,cdelt1,crval1,index,cal,crval4 = gbt_cal.getdata(filename)
 print "data loaded"
 unique_scan_n = unique(scan_n)
 # define windows 0 through 7
 windows = [((crval1 < 6.95e8 +2e6) & (crval1 > 6.95e8 -2e6)), ((crval1 < 7.25e8 +2e6) & (crval1 > 7.25e8 -2e6)), ((crval1 < 7.55e8 +2e6) & (crval1 > 7.55e8 -2e6)), ((crval1 < 7.85e8 +2e6) & (crval1 > 7.85e8 -2e6)), ((crval1 < 8.15e8 +2e6) & (crval1 > 8.15e8 -2e6)), ((crval1 < 8.45e8 +2e6) & (crval1 > 8.45e8 -2e6)), ((crval1 < 8.75e8 +2e6) & (crval1 > 8.75e8 -2e6)), ((crval1 < 9.05e8 +2e6) & (crval1 > 9.05e8 -2e6))]
 window_names = ['695','725','755','785','815','845','875','905']
 #polarization is crval4
 #xx=-5, yy=-6, xy=-7,yx=-8 
 mask_polx = crval4 == -5
 mask_poly = crval4 == -6
 mask_polxy = crval4 == -7
 mask_polyx = crval4 == -8
 mask_cal_off = cal == 'F'
 mask_cal_on = cal == 'T'
 #masks = [(mask_polx),(mask_poly),(mask_polxy),(mask_polyx)]
 masks = [(mask_polx),(mask_poly)]
 mask_names = ['x','y','xy', 'yx']
 #mask = mask_poly
 cut = 3
 for a,scan in enumerate(unique_scan_n):
  scan_mask = ( scan_n == scan )
  T_x_y = ones((2,spectra_clean[windows[0]*masks[0]*mask_cal_on*scan_mask].shape[0],spectra_clean[windows[0]*masks[0]*mask_cal_on*scan_mask].shape[1]))
  mask_x_y = T_x_y > 2
  for i,window in enumerate(windows):
   for j,mask in enumerate(masks):
    cal_on_temp = spectra_clean[window*mask*mask_cal_on*scan_mask]
    if ( isnan(cal_on_temp).sum() > 0 ):
       print "problem with cal on spectrum at:"
       nanmask = isnan(cal_on_temp)
       print where(nanmask)
       print "masking out problem for calibration"
       cal_on_temp = cal_on_temp[~nanmask[:,2],:]

    cal_on = cal_on_temp.mean(axis=0)
    cal_off_temp = spectra_clean[window*mask*mask_cal_off*scan_mask]
    if ( isnan(cal_off_temp).sum() > 0 ):
       print "problem with cal off spectrum at:"
       nanmask = isnan(cal_off_temp)
       print where(nanmask)
       print "masking out problem for calibration"
       cal_off_temp = cal_off_temp[~nanmask[:,2],:]

    cal_off = cal_off_temp.mean(axis=0)
    off = (cal_on + cal_off)/2.0
    on = (spectra_clean[window*mask*mask_cal_on*scan_mask] + spectra_clean[window*mask*mask_cal_off*scan_mask])/2.0
    #print cal_on.shape
    cal_ratio = cal_off/(cal_on-cal_off)
    #T_sys1 = ones(cal_ratio.shape)
    #for k in arange(ratio.shape[0]):
    T_sys1 =  T_noise[j*len(windows)+i] * cal_ratio  + T_noise[j*len(windows)+i]/2
    #print T_sys1.shape
    T_x_y[j] = (on - off)/off * T_sys1
    #print T_x_y.shape
    mask_x_y[j] = ((T_x_y[j] > cut) | (T_x_y[j] < -cut))
    mask_x_y[j][:,:100] = True
    mask_x_y[j][:,-100:] = True
    mask_x_y[j][isnan(T_x_y[j])] = True
    mask_rfi1d = ( mask_x_y[j].sum(axis=0) > mask_x_y[j].shape[0]/2 )
    mask_x_y[j][:,mask_rfi1d] = True
    print window_names[i],  mask_names[j], scan
    if (mask_x_y[j].sum()/(mask_x_y[j].shape[0]*mask_x_y[j].shape[1]) == 1):
       print " masked everything!  getting rid of mask to look at"
       mask_x_y[j] = T_x_y[j] > 1000
       print "T",T_x_y[j]
       print "cal on",cal_on
       print "cal off", cal_off
       print "on",on
       print "off",off
       print "T_sys", T_sys1
       mask_x_y[j][isnan(T_x_y[j])] = True 
       if (mask_x_y[j].sum()/(mask_x_y[j].shape[0]*mask_x_y[j].shape[1]) == 1):
         print "T is all nan, something really wrong!"
         continue

    print T_x_y[j][~mask_x_y[j]].max()
    print T_x_y[j][~mask_x_y[j]].min()
    print T_x_y[j][~mask_x_y[j]].std()
    #T_plot = T_x_y[j].copy()
    #T_plot[T>cut] = -cut
    #T_plot[T<-cut] = -cut
    #T_plot[mask_x_y[j]] = -cut
    #pylab.figure(1)
    #pylab.imshow(T_plot,cmap=pylab.cm.hot)
    #pylab.xlabel('frequency')
    #pylab.ylabel('time')
    #pylab.savefig('field4_'+str(cut)+'K_cut_tnoisecal'+window_names[i]+'_'+mask_names[j]+'_'+str(scan)+'.png', dpi=300)
    #pylab.clf()
    #ra = ra_gbt[window*mask*mask_cal_on*scan_mask]
    #print ra.shape
    #pylab.figure(4)
    #pylab.plot(ra,T_plot[:,900:1100].mean(axis=1))
    #pylab.xlabel('ra')
    #pylab.ylabel('temperature')
    #pylab.savefig('field4_'+window_names[i]+"_"+mask_names[j]+"_"+str(scan)+'.png', dpi=300)
    #pylab.clf()
    #freq = (index - crpix1[window*mask*scan_mask].mean())*cdelt1[window*mask*scan_mask].mean() + crval1[window*mask*scan_mask].mean()

 #pylab.ylim(-2,2)
   T_I = T_x_y[0] + T_x_y[1]
   mask_T_I = (mask_x_y[0] | mask_x_y[1])
   #T_I[mask_T_I] = -cut
   #pylab.figure(2)
   #pylab.imshow(T_I, cmap=pylab.cm.hot)
   #pylab.xlabel("freq channel")
   #pylab.ylabel("RA channel")
   #pylab.savefig("field4_T_intensity_"+window_names[i]+"_"+str(cut)+"K_cut_tnoisecal2_"+str(scan)+".png", dpi=300)
   #pylab.clf()
   #pylab.close()
   ra = ra_gbt[window*masks[0]*mask_cal_off*scan_mask]
   dec = dec_gbt[window*masks[0]*mask_cal_off*scan_mask]
   freq = (index - crpix1[window*masks[0]*scan_mask*mask_cal_off*scan_mask].mean())*cdelt1[window*masks[0]*scan_mask*mask_cal_off*scan_mask].mean() + crval1[window*masks[0]*scan_mask*mask_cal_off*scan_mask].mean()
   outfilename = "/home/scratch/kbandura/GBT08B_037/np_data/field4_"+day+"_" + window_names[i] +'_' + str(scan) +'_' +str(a) 
   save((outfilename+"_ra.npy"), ra)
   save((outfilename+"_dec.npy"), dec)
   save((outfilename+"_mask_T_I.npy"), mask_T_I)
   save((outfilename+"_T_I.npy"), T_I)
   save((outfilename+"_freq.npy"), freq)
 
 del spectra_clean   
   #print ra.shape
   #T_plotting = T_I[:,900:1100].mean(axis=1)
   #print T_plotting.shape
   #pylab.figure(4)
   #pylab.plot(ra,T_I[:,900:1100].mean(axis=1))
   #pylab.xlabel('ra')
   #pylab.ylabel('temperature')
   #pylab.savefig('field4_T_I_total_power2_'+window_names[i]+"_"+mask_names[j]+"_"+str(scan)+'.png', dpi=300)
   #pylab.clf()
 

