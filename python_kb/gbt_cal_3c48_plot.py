import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *

file = pyfits.open('/home/scratch/kbandura/GBT08B_037/aug30/d2_3c48_onoff_51-54.raw.acs.fits')

gbt_data = file[1].data

#for drift scan crval2 is az, crval3 is el
#for ra-long maps, crval2 is ra, crval3 is dec
#az_gbt = gbt_data.field('crval2')
#el_gbt = gbt_data.field('crval3')
#times = gbt_data.field('DATE-OBS')
#dur = gbt_data.field('DURATION')

#GBT = ephem.Observer()
#GBT.long = '-79:50:23.4'
#GBT.lat = '38:25:59.23'
#GBT.pressure = 0
#GBT.temp = 0
#az_r = az_gbt*pi/180.0
#el_r = el_gbt*pi/180.0

#max_times = times.shape[0]
#rag = zeros(max_times)
#decg = zeros(max_times)

#for i in range(0,max_times):
#     t1, t2 = times[i].split(".",1)
#     t21 = str(float("0."+t2)+dur[i]/2)
#     t23, t22 = t21.split(".",1)
#     t3 = time.strptime(t1, "%Y-%m-%dT%H:%M:%S")
#     t4 = time.strftime("%Y/%m/%d %H:%M:%S", t3)
#     GBT.date = t4 + "." + t22
#     rag[i], decg[i] = GBT.radec_of(az_r[i],el_r[i])

#ra_gbt = rag*180.0/pi
#dec_gbt = decg*180.0/pi

ra_gbt = gbt_data.field('CRVAL2')
dec_gbt = gbt_data.field('CRVAL3')


scan_n = gbt_data.field('SCAN')

unique_scan_n = unique(scan_n)

spectra = gbt_data.field('DATA')

spectra_clean = spectra.copy()

hanning = array( [0.25, 0.5, 0.25] )

#hanning filter for all spectra
for i in arange(spectra.shape[0]):
   spectra_clean[i] = convolve( spectra[i], hanning)[1:-1]
   

# crpix1 is the reference pixel for the frequency given
#crval1 is the frequency of the reference pixel
#cdelt1 is the difference in frequency between pixels
crpix1 = gbt_data.field('CRPIX1')  
cdelt1 = gbt_data.field('CDELT1')
crval1 = gbt_data.field('CRVAL1')
index = arange(0, len(spectra[1]))
#cheating here since the frequencies change a bit for each
#freq = (index - crpix1[1])*cdelt1[1] + crval1[1]

# T F if cal is on or off
cal = gbt_data.field('cal')
crval4 = gbt_data.field('crval4') #polarization

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
mask_scan_on = (scan_n == unique_scan_n[0]) | (scan_n==unique_scan_n[2])
mask_scan_off = (scan_n == unique_scan_n[1])  | (scan_n==unique_scan_n[3])

#masks = [(mask_polx*mask_cal_off*mask_scan_on),(mask_polx*mask_cal_on*mask_scan_on),(mask_polx*mask_cal_off*mask_scan_off),(mask_polx*mask_cal_on*mask_scan_off)]

#for mask in masks:
#  for window in windows:
#    freq = (index - crpix1[window*mask].mean())*cdelt1[window*mask].mean() + crval1[window*mask].mean()
#    pylab.plot(freq, spectra_clean[window*mask].mean(axis=0))

#pylab.ylim(0,30)
#pylab.show()

def F_3c48_model(freq):
   return( (freq/8.84207705e+04)**(-6.70431772e-01) )

ep = 390e-6
c = 299792458.0
nl = 0.99

## get T_sys use T_3c48 * (off/(on-off)) all with cal off
mask_on = (mask_poly*mask_cal_off*mask_scan_on)
mask_off = (mask_poly*mask_cal_off*mask_scan_off)
mask_noise_on = (mask_poly*mask_cal_on*mask_scan_off)
mask_noise_off = (mask_poly*mask_cal_off*mask_scan_off)

T_sys = []
T_noise = []

for window in windows:
   freq = (index - crpix1[window*mask_on].mean())*cdelt1[window*mask_on].mean() + crval1[window*mask_on].mean()
   F_3c48 = F_3c48_model(freq/1000000)
   na = 0.71*exp(-(4*pi*ep*freq/c)**2)
   to = 0.008 + exp(sqrt(freq/1e9))/8000.0
   KpJy = 2.85 * na * nl *exp(-to)
   #pylab.plot(freq, na)
   T_3c48 = KpJy*F_3c48
   on = spectra_clean[window*mask_on].mean(axis=0)
   off = spectra_clean[window*mask_off].mean(axis=0)
   T = T_3c48 * off/(on-off)
   T_sys.append(T)
   on = spectra_clean[window*mask_noise_on].mean(axis=0)
   off = spectra_clean[window*mask_noise_off].mean(axis=0)
   T_n = (on-off)/off * T
   T_noise.append(T_n)
   pylab.plot(freq, T_n)

pylab.title("y Noise Diode Temperature")
pylab.xlabel("freq(Hz)")
pylab.ylabel("Temperature (K)")
pylab.ylim(0,5)
pylab.savefig("y_noise_diode_temp.png", dpi=300)
pylab.show() 

pylab.clf()
for i,window in enumerate(windows):
   freq = (index - crpix1[window*mask_on].mean())*cdelt1[window*mask_on].mean() + crval1[window*mask_on].mean()
   pylab.plot(freq,T_sys[i])

pylab.title("y System Temperature")
pylab.xlabel("freq(Hz)")
pylab.ylabel("Temperature (K)")
pylab.ylim(0,40)
pylab.savefig("y_system_temp.png", dpi=300)

