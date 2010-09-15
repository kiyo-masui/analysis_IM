import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *

file = pyfits.open('/home/scratch/kbandura/GBT09C_075/00_zCOSMOS_drift_9-17.raw.acs.fits')

gbt_data = file[1].data

#for drift scan crval2 is az, crval3 is el
#for ra-long maps, crval2 is ra, crval3 is dec
az_gbt = gbt_data.field('crval2')
el_gbt = gbt_data.field('crval3')
times = gbt_data.field('DATE-OBS')
dur = gbt_data.field('DURATION')

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

#az_gbt = gbt_data.field('AZIMUTH')
#el_gbt = gbt_data.field('ELEVATIO')

scan_n = gbt_data.field('SCAN')

spectra = gbt_data.field('DATA')

hanning = array( [0.25, 0.5, 0.25] )

#az_f =[]
#el_f = []
#az_b = []
#el_b = []

####  when moving foward and backward (not drift)
#for i in range(0, (len(az_gbt)-1)):
#    if( ra_gbt[i+1] == ra_gbt[i]):
#       nothing = 1
#    elif( ra_gbt[i+1] > ra_gbt[i]):
#       az_f.append(az_gbt[i])
#       el_f.append(el_gbt[i])
#    else:
#       az_b.append(az_gbt[i])
#       el_b.append(el_gbt[i])


spectrum1 = convolve( spectra[1], hanning)[1:-1]

# crpix1 is the reference pixel for the frequency given
#crval1 is the frequency of the reference pixel
#cdelt1 is the difference in frequency between pixels
crpix1 = gbt_data.field('crpix1')  
cdelt1 = gbt_data.field('cdelt1')
crval1 = gbt_data.field('crval1')
index = arange(0, len(spectra[1]))
#cheating here since the frequencies change a bit for each
#freq = (index - crpix1[1])*cdelt1[1] + crval1[1]

cal = gbt_data.field('cal')
crval4 = gbt_data.field('crval4') #polarization

mask_freq = ((crval1 < 785000000 +1000000) & (crval1 > 785000000 -1000000))
#polarization is crval4
#xx=-5, yy=-6, xy=-7,yx=-8 
mask_polx = crval4 == -5
mask_poly = crval4 == -6
mask_polxy = crval4 == -7
mask_polyx = crval4 == -8
mask_cal = cal == 'F'
#mask_cal = cal == 'T'

#add in
subpix = crpix1[mask_cal*mask_polx*mask_freq]
subdelt = cdelt1[mask_cal*mask_polx*mask_freq]
subval = crval1[mask_cal*mask_polx*mask_freq]
index = arange(0, len(spectra[1]))
#each frequency would be this for each spectrum
freq = (index - subpix[0])*subdelt[0] + subval[0]

pol_set = array([spectra[mask_cal*mask_polx*mask_freq], spectra[mask_cal*mask_poly*mask_freq], spectra[mask_cal*mask_polxy*mask_freq], spectra[mask_cal*mask_polyx*mask_freq]])
el_set = array([el_gbt[mask_cal*mask_polx*mask_freq], el_gbt[mask_cal*mask_poly*mask_freq],el_gbt[mask_cal*mask_polxy*mask_freq],el_gbt[mask_cal*mask_polyx*mask_freq]])
az_set = array([az_gbt[mask_cal*mask_polx*mask_freq],az_gbt[mask_cal*mask_poly*mask_freq],ra_gbt[mask_cal*mask_polxy*mask_freq],az_gbt[mask_cal*mask_polyx*mask_freq]])
ra_set = array([ra_gbt[mask_cal*mask_polx*mask_freq],ra_gbt[mask_cal*mask_poly*mask_freq],ra_gbt[mask_cal*mask_polxy*mask_freq],ra_gbt[mask_cal*mask_polyx*mask_freq]])
dec_set = array([dec_gbt[mask_cal*mask_polx*mask_freq],dec_gbt[mask_cal*mask_poly*mask_freq],dec_gbt[mask_cal*mask_polxy*mask_freq],dec_gbt[mask_cal*mask_polyx*mask_freq]])

Tout = zeros((4,len(pol_set[1]), len(pol_set[1][0])))

Tsys = 29  #from GBT, need to get a better number

for j in arange(0,4):
  for i in arange(0, len(pol_set[j])):
    pol_set[j][i] = convolve(pol_set[j][i], hanning)[1:-1]

for j in arange(0,4):
  for i in arange(0, len(pol_set[j])):
    mset1 = sqrt((el_set[j][i]-el_set[j])**2+(az_set[j][i]-az_set[j])**2) < 0.0325
    mset2 = sqrt((ra_set[j][i]-ra_set[j])**2+(dec_set[j][i]-dec_set[j])**2) > 0.0325
    mset1[i] = False
    #print (mset1*mset2).sum()
    if ((mset1*mset2).sum() == 0):
        loc_avg = pol_set[j][mset1].mean(axis=0)
    else:
        loc_avg = pol_set[j][mset1*mset2].mean(axis=0)
    loc_mask = loc_avg == 0
    loc_avg[loc_mask] = 0.0000001
    Tout[j][i] = (pol_set[j][i] - loc_avg)/loc_avg * Tsys
    #print ra_set1[i], dec_set1[i]

pol_removed = zeros((len(pol_set[3]), len(pol_set[3][0])))

#ml is polarization fraction [0] is xx, [1] is yy, [2] is xy, [3] is yx

#ml = sqrt((out1 - out2)**2 + (out3 + out4)**2)/(out1+out2)
ml = sqrt((Tout[0] - Tout[1])**2 + (Tout[2] + Tout[3])**2 - (Tout[3] - Tout[2])**2)/(Tout[0]+Tout[1])

pol_mask = (ml > 0.05) #| (ml < -0.5)

#cheat and make everything a positive temperature
T_final = sqrt((Tout[0] + Tout[1])**2)

#pol_removed[pol_mask] = 0

#only use frequencies with smaller temperatures
mask_rfi = (T_final < 5) & (T_final > -5)

#bin up location based and frequency based
fe = 1420405750.0
c = 299792.458
H = 70.0

zed = fe/freq - 1.0

D = c/H *((zed+1.0)**2-1.0)/((zed+1.0)**2+1.0)

#make histogram, iterate this somehow
#pix1 = (subra < 36.8185318) & (subdec < 0.61790279 ) & redshiftbin
#gbt_hist = ones((2,5,10))
#gbt_hist[0,0] = T_final[pix1].mean()



#To get the correlation function need to look up a better way
#need to have imported deep2 histogram already, and have made one for gbt
#check if 3d works

def distance(x,y,z):
     return sqrt((x-centerx)**2  + (y-centery)**2 + (z-centerz)**2)

radii = array([1,2,3,4,5])
total = zeros(len(radii))
total_count = zeros(len(radii))

#for k in arange(len(radii)):
#    for i in arange(gbt_hist.shape[0]):
#       for j in arange(gbt_hist.shape[1]):
#         for l in arange(gbt_hist.shape[2]):
#          centerx = i
#          centery = j
#          centerz = l
#          dist = fromfunction(distance, gbt_hist.shape)
#          mask = (dist < radii[k]) & (dist >= radii[k]-1)
#          out = hist[i,j] * gbt_hist[mask]
#          total[k] += out.sum()
#          total_count[k] += len(out)




## for plotting set large lines to zero
mask_rfi_all = (T_final > 2) | (T_final < -2)
T_final[mask_rfi_all] = 0
#set1[mask_rfi_all] = 0

img_buf = ((T_final - T_final.min()) * 255.0/T_final.ptp()).astype(uint8)
#img_buf = ((set1 - set1.min()) * 255.0/set1.ptp()).astype(uint8) 
img_file = Image.fromarray(img_buf, 'L')
img_file.save('spectra.png')

x_plot = freq/1000000.0
y_plot = arange(0,T_final.shape[0])

pylab.imshow(T_final, extent=(min(x_plot), max(x_plot), min(y_plot), max(y_plot)),aspect=0.02)
pylab.xlabel('frequency')
pylab.ylabel('time') 
pylab.savefig('zcosmos_first_try.pdf')
pylab.clf()

U, s, Vh = linalg.svd(T_final, full_matrices=False)

s[:10] = 0

S = diag(s)

T_final_svd = dot(U,dot(S,Vh))
pylab.imshow(T_final_svd, extent=(min(x_plot), max(x_plot), min(y_plot), max(y_plot)),aspect=0.02)
pylab.xlabel('frequency')
pylab.ylabel('time') 
pylab.savefig('zcosmos_first_try_svd.pdf')
pylab.clf()

#pylab.plot(freq,ml[3])
pylab.plot(freq, T_final[3])
pylab.savefig('zcosmos_single_time_first.pdf')
pylab.clf()
