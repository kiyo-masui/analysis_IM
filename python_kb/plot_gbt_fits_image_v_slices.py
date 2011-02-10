import pylab
import pyfits
import sys
from numpy import *

filename1 = sys.argv[1]

filename1 = filename1.split('.')[0]

#filename1 = '12_wigglez1hr_azel_map'

#try:
#   slice = int(sys.argv[2])
#except:
#   slice = 50

hdudata = pyfits.open(filename1+'.fits')
# 0 is nothing, 1 is I, 2 is Q, 3 is U, 4 is V
hdr = hdudata[4].header
#RA I think, 3rd index in data
crval1 = hdr['CRVAL1']
crpix1 = hdr['CRPIX1']
cdelt1 = hdr['CDELT1']

#dec 2nd
crval2 = hdr['CRVAL2']
crpix2 = hdr['CRPIX2']
cdelt2 = hdr['CDELT2']
#freq 1st index in file
crval3 = hdr['CRVAL3']
crpix3 = hdr['CRPIX3']
cdelt3 = hdr['CDELT3']

image_cube = hdudata[4].data

slices = image_cube.shape[0]
print slices

ra = ( arange(image_cube.shape[2]) + 1 - crpix1)*cdelt1 + crval1
dec = ( arange(image_cube.shape[1]) +1 - crpix2)*cdelt2 + crval2
freqs = ( arange(image_cube.shape[0]) +1 - crpix3)*cdelt3 + crval3
freqs = freqs/1e6

for slice, freq in enumerate(freqs):
   nancut = (image_cube[slice] < 10e10) & ( image_cube[slice] != NaN )
   #print nancut
   print image_cube[slice][nancut].std()
   print median(image_cube[slice][nancut])
   cut = ( image_cube[slice] > 3.0*image_cube[slice][nancut].std() ) 
   image_cube[slice][cut] = 3.0*image_cube[slice][nancut].std()
   cut = ( image_cube[slice] < -3.0*image_cube[slice][nancut].std() ) 
   image_cube[slice][cut] = -3.0*image_cube[slice][nancut].std()
   medianv = median(image_cube[slice][nancut])
   for element in image_cube[slice][cut]:
      print element
   pylab.imshow(image_cube[slice], interpolation='gaussian', extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.savefig(filename1+str(freq)[:3]+'.png')
   pylab.clf()
