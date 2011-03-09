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

# Alternate plotting of frequency for given ra/dec
#for rind,r in enumerate(ra):
#   for dind,d in enumerate(dec):
#       nancut = (image_cube[:,dind,rind]<10e10) & ( image_cube[:,dind,rind] !=NaN)
#       image_cube[:,dind,rind][~nancut] = 0
#       pylab.plot(freqs,image_cube[:,dind,rind])
#       pylab.xlabel('frequency (MHz)')
#       pylab.ylabel('Temperature (K)')
#       pylab.title('Stokes V Temperature for RA = '+str(r)[:3]+' and DEC = '+str(d)[:3])
#       pylab.xlim(670,928)
#       pylab.ylim( -10,10 )
#       pylab.savefig('v_'+filename1+'_'+str(r)[:3]+'_'+str(d)[:3]+'.png')
#       pylab.clf()


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
# Alternate plotting command to set temperature limits
#   pylab.imshow(image_cube[slice], interpolation='gaussian',vmin=-0.5, vmax=0.5, extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.imshow(image_cube[slice], interpolation='gaussian', extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.colorbar()
   pylab.savefig('v_'+filename1+str(freq)[:3]+'.png')
   pylab.clf()
