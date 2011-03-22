import pylab
import pyfits
import sys
from numpy import *
import core.algebra as al

filename1 = sys.argv[1]

#filename1 = filename1.split('.')[0]

array = al.load(filename1)
array = al.make_vect(array)

#   creates a 3D array with indices (freq, ra, dec)

ra = array.calc_axis('ra')
dec = array.calc_axis('dec')
freqs = array.calc_axis('freq')
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
#   pylab.imshow(image_cube[slice], interpolation='gaussian',vmin=-2, vmax=2, extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.imshow(image_cube[slice], interpolation='gaussian', extent=(ra.max(),ra.min(),dec.min(),dec.max()), origin='lower')
   pylab.colorbar()
   pylab.savefig(filename1+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename1+str(freq)[:3]+'.png')
   pylab.clf()
