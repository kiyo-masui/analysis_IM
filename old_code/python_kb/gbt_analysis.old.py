import pyfits
import pylab
from PIL import Image
from numpy import *

file = pyfits.open('/home/scratch/kbandura/nov06_08/d3_field3cdec2_fs60sec_12.5MHzBW_32-52.raw.acs.fits')

gbt_data = file[1].data

ra_gbt = gbt_data.field('crval2')
dec_gbt = gbt_data.field('crval3')

az_gbt = gbt_data.field('AZIMUTH')
el_gbt = gbt_data.field('ELEVATIO')

scan_n = gbt_data.field('SCAN')

spectra = gbt_data.field('DATA')

hanning = array( [0.25, 0.5, 0.25] )

az_f =[]
el_f = []
az_b = []
el_b = []


for i in range(0, (len(az_gbt)-1)):
    if( ra_gbt[i+1] == ra_gbt[i]):
       nothing = 1
    elif( ra_gbt[i+1] > ra_gbt[i]):
       az_f.append(az_gbt[i])
       el_f.append(el_gbt[i])
    else:
       az_b.append(az_gbt[i])
       el_b.append(el_gbt[i])


spectrum1 = convolve( spectra[1], hanning)[1:-1]

crpix1 = gbt_data.field('crpix1')  
cdelt1 = gbt_data.field('cdelt1')
crval1 = gbt_data.field('crval1')
index = arange(0, len(spectra[1]))
freq = (index - crpix1[1])*cdelt1[1] + crval1[1]

cal = gbt_data.field('cal')
crval4 = gbt_data.field('crval4') #polarization

mask_freq = ((crval1 < crval1[2] +1000000) & (crval1 > crval1[2] -1000000))
mask_pol = crval4 == -5
mask_cal = cal == 'T'

set1 = spectra[mask_cal*mask_pol*mask_freq]

el_set1 = el_gbt[mask_cal*mask_pol*mask_freq]
az_set1 = az_gbt[mask_cal*mask_pol*mask_freq]
ra_set1 = ra_gbt[mask_cal*mask_pol*mask_freq]
dec_set1 = dec_gbt[mask_cal*mask_pol*mask_freq]


set2 = zeros((len(set1), len(set1[0])))

for i in arange(0, len(set1)):
    set1[i] = convolve(set1[i], hanning)[1:-1]

for i in arange(0, len(set1)):
    mset1 = sqrt((el_set1[i]-el_set1)**2+(az_set1[i]-az_set1)**2) < 0.0325
    mset2 = sqrt((ra_set1[i]-ra_set1)**2+(dec_set1[i]-dec_set1)**2) > 0.0325
    mset1[i] = False
    #print (mset1*mset2).sum()
    if ((mset1*mset2).sum() == 0):
        loc_avg = set1[i]
    else:
        loc_avg = set1[mset1*mset2].mean(axis=0)
    loc_mask = loc_avg == 0
    loc_avg[loc_mask] = 0.0000001
    set2[i] = (set1[i] - loc_avg)/loc_avg
    #print ra_set1[i], dec_set1[i]

mask_rfi = (set2[3] < 0.6) & (set2[3] > -0.6)
mask_rfi_all = (set2 > 0.6) | (set2 < -0.6)

set2[mask_rfi_all] = 0

img_buf = ((set2 - set2.min()) * 255.0/set2.ptp()).astype(uint8)
img_file = Image.fromarray(img_buf, 'L')
img_file.save('spectra.png')

pylab.plot(freq[mask_rfi],set2[3][mask_rfi])
pylab.show()
