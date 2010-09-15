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
#mask_scan_on = (scan_n == 11) | (scan_n==13)
#mask_scan_off = (scan_n == 12)  | (scan_n==14)
#mask_scan_on = (scan_n == 5) | (scan_n==7)
#mask_scan_off = (scan_n == 6)  | (scan_n==8)


#masks = [(mask_polx*mask_cal_off*mask_scan_on),(mask_polx*mask_cal_on*mask_scan_on),(mask_polx*mask_cal_off*mask_scan_off),(mask_polx*mask_cal_on*mask_scan_off)]
#masks = [(mask_polx),(mask_poly),(mask_polxy),(mask_polyx)]
masks = [(mask_polx),(mask_poly)]
mask_names = ['x','y','xy', 'yx']
#mask = mask_poly

#T_sky = zeros((len(masks),len(windows),spectra_clean[windows[0]*masks[0]*mask_cal_on].shape[0],spectra_clean[windows[0]*masks[0]*mask_cal_on].shape[1]))
#print T_sky.shape
T_sky=[]
T_sky_mask_rfi = []

cut = 5


for scan in unique_scan_n:
 scan_mask = ( scan_n == scan )
 for j,mask in enumerate(masks):
  for i,window in enumerate(windows):
    on = spectra_clean[window*mask*mask_cal_on*scan_mask]
    off = spectra_clean[window*mask*mask_cal_off*scan_mask]
    ratio = off/(on-off)
    T = ones(ratio.shape)
    for k in arange(ratio.shape[0]):
      T[k] =  T_noise[j*len(windows)+i] * ratio[k] - T_sys[j*len(windows)+i]
    
    print T.shape
    mask_rfi = ((T > cut) | (T < -cut))
    mask_rfi[:,:100] = True
    mask_rfi[:,-100:] = True
    mask_rfi[isnan(T)] = True
    print T[~mask_rfi].max()
    print T[~mask_rfi].min()
    T_plot = T.copy()
    #T_plot[T>cut] = -cut
    #T_plot[T<-cut] = -cut
    T_plot[mask_rfi] = -cut
    pylab.figure(1)
    pylab.imshow(T_plot,cmap=pylab.cm.hot)
    pylab.xlabel('frequency')
    pylab.ylabel('time')
    pylab.savefig('field4_'+str(cut)+'K_cut'+window_names[i]+'_'+mask_names[j]+'_'+str(scan)+'.png', dpi=300)
    pylab.clf()
    #pylab.close()
    #T_sky.append(T)
    #T_sky_mask_rfi.append(mask_rfi)
    freq = (index - crpix1[window*mask*scan_mask].mean())*cdelt1[window*mask*scan_mask].mean() + crval1[window*mask*scan_mask].mean()
    #T = T_sky[i+j*len(windows)]
    #mask_rfi = T_sky_mask_rfi[i+j*len(windows)]
    print mask_rfi.shape
    # is rfi if half the time is greater than cut temperature
    mask_rfi1d = ( mask_rfi.sum(axis=0) > mask_rfi.shape[0]/2 )
    print mask_rfi1d.shape
    mask_rfi[:,mask_rfi1d] = True
    T_avg = ones(T.shape[1])
    #T_temp = ones(T.shape[0])
    for k in arange(0, T.shape[1]):
         T_temp2 = T[:,k]
         mask1chan = mask_rfi[:,k]
         T_avg[k] = T_temp2[~mask1chan].mean()
    pylab.figure(2)
    pylab.plot(freq[~mask_rfi1d],T_avg[~mask_rfi1d])
    #deltaT = T/T_avg - 1.0
    deltaT = ones(T.shape)
    for k in arange(T.shape[0]):
       deltaT[k] = T[k] - T_avg

    #mask_rfi2 = ( deltaT > 1.0) | (deltaT < -1.0 )
    #mask_rfi = mask_rfi2 + mask_rfi
    print deltaT[~mask_rfi].min()
    print deltaT[~mask_rfi].max()
    print deltaT[~mask_rfi].std()
    deltaT[mask_rfi] = deltaT[~mask_rfi].min()
    #deltaT[:,mask_rfi1d] = 0
    pylab.figure(3)
    pylab.imshow(deltaT, cmap=pylab.cm.hot)
    pylab.savefig('field4_deltaT_'+str(cut)+'K_cut'+window_names[i]+'_'+mask_names[j]+'_'+str(scan)+'.png', dpi=300)
    pylab.clf()

  #pylab.ylim(-2,2)
  pylab.figure(2)
  pylab.savefig("field4_t_sky_" + mask_names[j] +'_'+str(cut)+"K_cut_"+str(scan)+".png", dpi=300)
  pylab.clf()
  pylab.close()
    

#pylab.imshow(T_sky[3], extent=(min(x_plot), max(x_plot), min(y_plot), max(y_plot)),aspect=0.02)
#pylab.imshow(T_sky[3],aspect=0.02)
#pylab.xlabel('frequency')
#pylab.ylabel('time')
#pylab.savefig('field4_first_try.png', dpi=300)
#pylab.clf()

'''for scan in unique_scan_n:
 scan_mask = ( scan_n == scan )
 for j,mask in enumerate(masks):
  for i,window in enumerate(windows):
      freq = (index - crpix1[window*mask*scan_mask].mean())*cdelt1[window*mask*scan_mask].mean() + crval1[window*mask*scan_mask].mean()
      T = T_sky[i+j*len(windows)]
      mask_rfi = T_sky_mask_rfi[i+j*len(windows)]
      print mask_rfi.shape
      # is rfi if half the time is greater than cut temperature
      mask_rfi1d = ( mask_rfi.sum(axis=0) > mask_rfi.shape[0]/2 )
      T_plot = ones(T.shape[1])
      T_temp = ones(T.shape[0])
      for k in arange(0, T.shape[1]):
         T_temp = T[:,k]
         mask1chan = mask_rfi[:,k]
         T_plot[k] = T_temp[~mask1chan].mean() 
      pylab.plot(freq[~mask_rfi1d],T_plot[~mask_rfi1d])
  
  #pylab.ylim(-2,2)
  pylab.savefig("field4_t_sky_" + mask_names[j] +'_'+str(cut)+"K_cut.png", dpi=300)
  pylab.clf()
  pylab.close()



#pylab.ylim(-2,2)
#pylab.show()'''
