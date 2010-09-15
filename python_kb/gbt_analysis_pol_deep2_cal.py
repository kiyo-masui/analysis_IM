import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *
import gbt_cal_3c48

cal_filename = '/home/scratch/kbandura/GBT08B_037/aug30/d2_3c48_onoff_51-54.raw.acs.fits'

T_noise,T_sys =  gbt_cal_3c48.getcal_3c48(cal_filename)

print T_noise.shape
print T_sys.shape

print "finished loading cal"

filename = '/home/scratch/kbandura/GBT08B_037/aug30/d2_field4_drift_55-60.raw.acs.fits'

spectra_clean,ra_gbt,dec_gbt,scan_n,crpix1,cdelt1,crval1,index,cal,crval4 = gbt_cal_3c48.getdata(filename)

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
masks = [(mask_polx),(mask_poly),mask_polxy,mask_polyx]
mask_names = ['x','y','xy', 'yx']
#mask = mask_poly

cut = 20



for scan in unique_scan_n:
 scan_mask = ( scan_n == scan )
 T_x_y = ones((len(masks),spectra_clean[windows[0]*masks[0]*mask_cal_on*scan_mask].shape[0],spectra_clean[windows[0]*masks[0]*mask_cal_on*scan_mask].shape[1]))
 mask_x_y = T_x_y > 2
 G2 = ones((2,spectra_clean[windows[0]*masks[0]*mask_cal_on*scan_mask].shape[1]))
 for i,window in enumerate(windows):
   for j,mask in enumerate(masks):
    cal_on = spectra_clean[window*mask*mask_cal_on*scan_mask].mean(axis=0)
    cal_off = spectra_clean[window*mask*mask_cal_off*scan_mask].mean(axis=0)
    off = (cal_on + cal_off)/2.0
    on = (spectra_clean[window*mask*mask_cal_on*scan_mask] + spectra_clean[window*mask*mask_cal_off*scan_mask])/2.0
    #print cal_on.shape
    cal_ratio = cal_off/(cal_on-cal_off)
    #T_sys1 = ones(cal_ratio.shape)
    #for k in arange(ratio.shape[0]):
    T_sys1 =  T_noise[j*len(windows)+i] * cal_ratio  + T_noise[j*len(windows)+i]/2
    if ( j < 2 ):
      G2[j] = off/(T_sys1 + T_noise[j*len(windows)+i]/2)
      T_x_y[j] = on/G2[j] - T_sys1 - T_noise[j*len(windows)+i]/2
    else:
      T_x_y[j] = on/(sqrt(G2[0])*sqrt(G2[1])) - T_noise[j*len(windows)+i]/2
    #print T_x_y.shape
    mask_x_y[j] = ((T_x_y[j] > cut) | (T_x_y[j] < -cut))
    mask_x_y[j][:,:100] = True
    mask_x_y[j][:,-100:] = True
    mask_x_y[j][isnan(T_x_y[j])] = True
    mask_rfi1d = ( mask_x_y[j].sum(axis=0) > mask_x_y[j].shape[0]/2 )
    mask_x_y[j][:,mask_rfi1d] = True
    if (mask_x_y[j].sum()/(mask_x_y[j].shape[0]*mask_x_y[j].shape[1]) == 1):
       print " masked everything!  getting rid of mask to look at"
       mask_x_y[j] = T_x_y[j] > 1000
       print T_x_y[j]
       print cal_on
       print cal_off
       print on
       print off
       print T_sys1
       mask_x_y[j][isnan(T_x_y[j])] = True 
       if (mask_x_y[j].sum()/(mask_x_y[j].shape[0]*mask_x_y[j].shape[1]) == 1):
         print " all nan, something really wrong!"
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
    

   freq = (index - crpix1[window*mask*scan_mask].mean())*cdelt1[window*mask*scan_mask].mean() + crval1[window*mask*scan_mask].mean()
   #pylab.ylim(-2,2)
   T_I = T_x_y[0] + T_x_y[1]
   T_Q = T_x_y[0] - T_x_y[1]
   T_U = T_x_y[2] + T_x_y[3]
   ml = sqrt( T_Q**2 + T_U**2 )/T_I
   #mask_T_I = (mask_x_y[0] | mask_x_y[1])
   mask_pol = ml > 0.05
   T_I[mask_pol] = -cut
   T_I[isnan(T_I)] = -cut
   T_I[:,:100] = -cut
   T_I[:,-100:] = -cut
   T_I[((T_I > cut ) | (T_I < - cut ))] = -cut
   pylab.figure(2)
   pylab.imshow(T_I, cmap=pylab.cm.hot)
   pylab.xlabel("freq channel")
   pylab.ylabel("RA channel")
   pylab.savefig("field4_T_intensity_"+window_names[i]+"_"+str(cut)+"K_cut_polcut_tnoisecal_"+str(scan)+".png", dpi=300)
   pylab.clf()
   pylab.close()
   ra = ra_gbt[window*masks[0]*mask_cal_on*scan_mask]
   print ra.shape
   T_plotting = T_I[:,900:1100].mean(axis=1)
   print T_plotting.shape
   pylab.figure(4)
   pylab.plot(ra,T_I[:,900:1100].mean(axis=1))
   pylab.xlabel('ra')
   pylab.ylabel('temperature')
   pylab.savefig('field4_T_I_total_power'+window_names[i]+"_"+mask_names[j]+"_"+str(scan)+'.png', dpi=300)
   pylab.clf()
    

