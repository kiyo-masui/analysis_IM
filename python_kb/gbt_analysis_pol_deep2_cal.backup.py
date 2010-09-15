import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *
import gbt_cal_3c48

T_noise,T_sys =  gbt_cal_3c48.getcal()

#file = pyfits.open('/home/scratch/kbandura/GBT09C_075/00_3c218_onoff_scan5-8.raw.acs.fits')
file = pyfits.open('/home/scratch/kbandura/GBT08B_037/aug29/d1_field3_drift_21-26.raw.acs.fits')

gbt_data = file[1].data

ctype2 = gbt_data.field('CTYPE2').copy()

if ( ctype2[0] == "AZ"):
  #for drift scan crval2 is az, crval3 is el
  #for ra-long maps, crval2 is ra, crval3 is dec
  az_gbt = gbt_data.field('crval2').copy()
  el_gbt = gbt_data.field('crval3').copy()
  times = gbt_data.field('DATE-OBS').copy()
  dur = gbt_data.field('DURATION').copy()
  GBT = ephem.Observer()
  GBT.long = '-79:50:23.4'
  GBT.lat = '38:25:59.23'
  GBT.pressure = 0
  GBT.temp = 0
  az_r = az_gbt*pi/180.0
  el_r = el_gbt*pi/180.0
  max_times = times.shape[0]
  rag = zeros(max_times)
  decg = zeros(max_times)
  for i in range(0,max_times):
     t1, t2 = times[i].split(".",1)
     t21 = str(float("0."+t2)+dur[i]/2)
     t23, t22 = t21.split(".",1)
     t3 = time.strptime(t1, "%Y-%m-%dT%H:%M:%S")
     t4 = time.strftime("%Y/%m/%d %H:%M:%S", t3)
     GBT.date = t4 + "." + t22
     rag[i], decg[i] = GBT.radec_of(az_r[i],el_r[i])

  ra_gbt = rag*180.0/pi
  dec_gbt = decg*180.0/pi

else:
  ra_gbt = gbt_data.field('CRVAL2').copy()
  dec_gbt = gbt_data.field('CRVAL3').copy()


scan_n = gbt_data.field('SCAN').copy()

spectra = gbt_data.field('DATA')

spectra_clean = spectra.copy()

hanning = array( [0.25, 0.5, 0.25] )

#hanning filter for all spectra
for i in arange(spectra.shape[0]):
   spectra_clean[i] = convolve( spectra[i], hanning)[1:-1]
   
spectra = 1

# crpix1 is the reference pixel for the frequency given
#crval1 is the frequency of the reference pixel
#cdelt1 is the difference in frequency between pixels
crpix1 = gbt_data.field('CRPIX1').copy()  
cdelt1 = gbt_data.field('CDELT1').copy()
crval1 = gbt_data.field('CRVAL1').copy()
index = arange(0, len(spectra_clean[1]))

#cheating here since the frequencies change a bit for each
#freq = (index - crpix1[1])*cdelt1[1] + crval1[1]

# T F if cal is on or off
cal = gbt_data.field('cal').copy()
crval4 = gbt_data.field('crval4').copy() #polarization

file.close()
#gbt_data.close()

gbt_data = 1 

# define windows 0 through 7
windows = [((crval1 < 6.95e8 +1e6) & (crval1 > 6.95e8 -1e6)), ((crval1 < 7.25e8 +1e6) & (crval1 > 7.25e8 -1e6)), ((crval1 < 7.55e8 +1e6) & (crval1 > 7.55e8 -1e6)), ((crval1 < 7.85e8 +1e6) & (crval1 > 7.85e8 -1e6)), ((crval1 < 8.15e8 +1e6) & (crval1 > 8.15e8 -1e6)), ((crval1 < 8.45e8 +1e6) & (crval1 > 8.45e8 -1e6)), ((crval1 < 8.75e8 +1e6) & (crval1 > 8.75e8 -1e6)), ((crval1 < 9.05e8 +1e6) & (crval1 > 9.05e8 -1e6))]

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
mask = mask_polx

T_sky = []
#for mask in masks:
for i,window in enumerate(windows):
    on = spectra_clean[window*mask*mask_cal_on]
    off = spectra_clean[window*mask*mask_cal_off]
    ratio = off/(on-off)
    Ta =  T_noise[i] * ratio
    T =  Ta - T_sys[i]
    T_sky.append(T)
    
j=0
#for j,mask in enumerate(masks):
for i,window in enumerate(windows):
      freq = (index - crpix1[window*mask].mean())*cdelt1[window*mask].mean() + crval1[window*mask].mean()
      T = T_sky[i+j*len(windows)]
      pylab.plot(freq,T.mean(axis=0))

#pylab.ylim(0,30)
pylab.show()
