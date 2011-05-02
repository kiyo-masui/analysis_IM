import pylab
import pyfits
import sys
from numpy import *
import core.algebra as al
import scipy

filename1 = sys.argv[1]

filename2 = filename1.split('.')[0]

array = al.load(filename1)
array = al.make_vect(array)

#   creates a 3D array with indices (freq, ra, dec)

ra = array.get_axis('ra')
dec = array.get_axis('dec')
freqs = array.get_axis('freq')
freqs = freqs/1e6


for slice, freq in enumerate(freqs):
   nancut = (array[slice] < 10e10) & ( array[slice] != NaN )
#   print nancut
#   print image_cube[slice][nancut].std()
#   print median(image_cube[slice][nancut])
   cut = ( array[slice] > 3.0*array[slice][nancut].std() ) 
   array[slice][cut] = 3.0*array[slice][nancut].std()
   cut = ( array[slice] < -3.0*array[slice][nancut].std() ) 
   array[slice][cut] = -3.0*array[slice][nancut].std()

#   Need to rotate array[slice] because axes were flipped
#   print array[0]
   new_array = scipy.transpose(array[slice])
#   print new_array
#   medianv = median(array[slice][nancut])
#   for element in image_cube[slice][cut]:
#      print element
#   Alternate plotting command to set temperature limits
   pylab.imshow(new_array, interpolation='gaussian', vmin=-0.7, vmax=0.7, extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
#   pylab.xlabel('Dec')
#   pylab.ylabel('RA')
#   pylab.imshow(new_array, interpolation='gaussian', extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.colorbar() #For some reason this isn't working, fixed...
   pylab.savefig(filename2+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()
