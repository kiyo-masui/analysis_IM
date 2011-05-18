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

#print array
#   creates a 3D array with indices (freq, ra, dec)

ras = array.get_axis('ra')
decs = array.get_axis('dec')
freqs = array.get_axis('freq')
freqs = freqs/1e6


for slice, freq in enumerate(freqs):
   nancut = (array[slice] < 10e10) & ( array[slice] != NaN )
   cut = ( array[slice] > 3.0*array[slice][nancut].std() ) 
   array[slice][cut] = 3.0*array[slice][nancut].std()
   cut = ( array[slice] < -3.0*array[slice][nancut].std() ) 
   array[slice][cut] = -3.0*array[slice][nancut].std()

#   Need to rotate array[slice] because axes were flipped
   new_array = scipy.transpose(array[slice])
   medianv = median(array[slice][nancut])

#   Alternate plotting command to set temperature limits
   pylab.imshow(new_array, cmap='bone', vmin=-0.15, vmax=0.15, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.imshow(new_array, interpolation='gaussian', cmap='hot', extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
   pylab.colorbar() #For some reason this isn't working, fixed...
   pylab.savefig(filename2+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()

#for slice, dec in enumerate(decs):
#   print slice
#   print dec
#   nancut = (array[:,:,slice] < 10e10) & ( array[:,:,slice] != NaN )
#   cut = ( array[:,:,slice] > 3.0*array[:,:,slice][nancut].std() )
#   array[:,:,slice][cut] = 3.0*array[:,:,slice][nancut].std()
#   cut = ( array[:,:,slice] < -3.0*array[:,:,slice][nancut].std() )
#   array[:,:,slice][cut] = -3.0*array[:,:,slice][nancut].std()

#   print array[:,:,slice]
#   Need to rotate array[slice] because axes were flipped
#   new_array = scipy.transpose(array[:,:,slice])
#   medianv = median(array[slice][nancut])
 
#   Alternate plotting command to set temperature limits
#   pylab.imshow(array[:,:,slice], aspect = 0.02, interpolation='gaussian', vmin=-0.2, vmax=0.2, extent=(ras.max(),ras.min(),freqs.min(),freqs.max()), origin='lower')
#   pylab.imshow(array[:,:,slice], interpolation='gaussian', extent=(freqs.max(),freqs.min(),ras.min(),ras.max()), origin='lower')
#   pylab.colorbar() #For some reason this isn't working, fixed...
#   pylab.savefig(filename2+str(dec)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
#   pylab.clf()


