#This script generates a set of maps in RA and DEC for each frequency in the map data. Can be used to make plots of maps.

import pylab
import pyfits
import sys
from numpy import *
#from core import algebra as al
#from utils import * 
import core.algebra as al
import scipy

filename1 = sys.argv[1]
filename2 = sys.argv[2]

filename1_trunc = filename1.split('.')[0]
filename2_trunc = filename2.split('.')[0]

array1 = al.load(filename1)
array1 = al.make_vect(array1)
array2 = al.load(filename2)
array2 = al.make_vect(array2)


#print array
#   creates a 3D array with indices (freq, ra, dec)

ras = array1.get_axis('ra')
decs = array1.get_axis('dec')
freqs = array1.get_axis('freq')
freqs = freqs/1e6


for slice, freq in enumerate(freqs):
   nancut = (array1[slice] < 10e10) & ( array1[slice] != NaN )
   cut = ( array1[slice] > 3.0*array1[slice][nancut].std() ) 
   array1[slice][cut] = 3.0*array1[slice][nancut].std()
   cut = ( array1[slice] < -3.0*array1[slice][nancut].std() ) 
   array1[slice][cut] = -3.0*array1[slice][nancut].std()

#   Need to rotate array[slice] because axes were flipped
   new_array1 = scipy.transpose(array1[slice])
   medianv1 = median(array1[slice][nancut])

#   Alternate plotting command to set temperature limits
   pylab.imshow(new_array1, cmap='hot', vmin=-1.0, vmax=1.0, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.imshow(new_array1, interpolation='gaussian', cmap='hot', extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
   pylab.colorbar() #For some reason this isn't working, fixed...
   pylab.title('Frequency: '+str(freq)[:5]+' MHz')
   pylab.savefig(filename1_trunc+ '_' +str(freq*100)[:5]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()

   nancut = (array2[slice] < 10e10) & ( array2[slice] != NaN )
   cut = ( array2[slice] > 3.0*array2[slice][nancut].std() )
   array2[slice][cut] = 3.0*array2[slice][nancut].std()
   cut = ( array2[slice] < -3.0*array2[slice][nancut].std() )
   array2[slice][cut] = -3.0*array2[slice][nancut].std()

#   Need to rotate array[slice] because axes were flipped
   new_array2 = scipy.transpose(array2[slice])
   medianv2 = median(array2[slice][nancut])

#   Alternate plotting command to set temperature limits
   pylab.imshow(new_array2, cmap='hot', vmin=-1.0, vmax=1.0, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.imshow(new_array2, interpolation='gaussian', cmap='hot', extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
   pylab.colorbar() #For some reason this isn't working, fixed...
   pylab.title('Frequency: '+str(freq)[:5]+' MHz')

   pylab.savefig(filename2_trunc+ '_' +str(freq*100)[:5]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()

   pylab.imshow(new_array1-new_array2, interpolation='gaussian', cmap='hot', vmin=-0.5,vmax=0.5, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
   pylab.colorbar() #For some reason this isn't working, fixed...
   pylab.title('Frequency: '+str(freq)[:5]+' MHz')
   pylab.savefig(filename1_trunc+'_diff_' +str(freq*100)[:5]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
   pylab.clf()



#Alternate code if want to plot in term of dec instead of freq.
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


