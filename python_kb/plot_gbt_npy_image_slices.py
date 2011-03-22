import pylab
import pyfits
import sys
from numpy import *
import core.algebra as al

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
#   medianv = median(array[slice][nancut])
#   for element in image_cube[slice][cut]:
#      print element
#   Alternate plotting command to set temperature limits
#   pylab.imshow(array[slice], interpolation='gaussian',vmin=-2, vmax=2, extent=(dec.max(),dec.min(),ra.min(),ra.max()), origin='lower')
   pylab.imshow(array[slice], interpolation='gaussian', extent=(dec.max(),dec.min(),ra.min(),ra.max()), origin='lower')
   pylab.colorbar()
   pylab.savefig(filename2+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()
