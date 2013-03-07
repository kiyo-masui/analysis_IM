#This script generates a set of maps in RA and DEC for each frequency in the map data. Can be used to make plots of maps.

import matplotlib.pyplot as plt
import pyfits
import sys
from numpy import *
import core.algebra as al
import scipy
import scipy.optimize as opt

def function(freq,params):
    base = 750./freq
    function = params[0]*power(base,params[1])
    return function

def residuals(params,freq,data):
    err = data-function(freq,params)
    return err

directory = sys.argv[1]
#filename1 = sys.argv[1]

file_base = 'fir_15hr_absorb_clean_map_I_'

file_sets = ('749','768','788','807','827')

full_data = []
full_freqs = []
for i in range(0,len(file_sets)):
    array = al.load(directory+file_base+file_sets[i]+'.npy')
    array = al.make_vect(array)
    #  creates a 3D array with indices (freq, ra, dec)
    ras = array.get_axis('ra')
    decs = array.get_axis('dec')
    freqs = array.get_axis('freq')
    freqs = freqs/1e6
    full_data.append(array)
    full_freqs.append(freqs)

#print shape(full_data)
#print shape(full_freqs)

list_freq = zeros(2000)
list_array = zeros((2000,16,10))
for i in range(0,2000):
    if i<400:
        list_freq[i]=full_freqs[0][399-i]
        list_array[i,:,:] = full_data[0][399-i]
    elif i<800:
        list_freq[i]=full_freqs[1][399-(i-400)]
        list_array[i,:,:] = full_data[1][799-i]
    elif i<1200:
        list_freq[i]=full_freqs[2][399-(i-800)]
        list_array[i,:,:] = full_data[2][1199-i]
    elif i<1600:
        list_freq[i]=full_freqs[3][399-(i-1200)]
        list_array[i,:,:] = full_data[3][1599-i] 
    else:
        list_freq[i]=full_freqs[4][399-(i-1600)]
        list_array[i,:,:] = full_data[4][1999-i]

fits = zeros((2,16,10))
for i in range(0,16):
    for j in range(0,10):
        data = list_array[:,i,j]
        params0 = [-0.5,2.5]
        plsq = opt.leastsq(residuals,params0,args=(list_freq,data),full_output=0,maxfev=5000)
#        print plsq
        fits[0,i,j] = plsq[0][0]
        fits[1,i,j] = plsq[0][1]


plt.subplot(121)
plt.imshow(fits[0],interpolation='nearest')
plt.colorbar()
plt.title('Amplitude')
plt.subplot(122)
plt.imshow(fits[1],interpolation='nearest')
plt.colorbar()
plt.title('Power')
#plt.show()
#print fits
plt.clf()

#plt.plot(list_freq)
#plt.show()
fit_data = zeros((2000,16))
for i in range(0,16):
    fit_data[:,i] = function(list_freq,list_array[:,i,5])

plt.subplot(121)
plt.imshow(list_array[:,:,5],aspect=16.0/2000,vmin=-0.6,vmax=0.6,interpolation='nearest')
plt.colorbar()
plt.title('Before Subtraction')
plt.subplot(122)
plt.imshow(list_array[:,:,5]-fit_data,aspect=16.0/2000,vmin=-0.1,vmax=0.1,interpolation='nearest')
plt.title('After Subtraction')
plt.colorbar()
plt.show()

#for slice, freq in enumerate(freqs):
#   nancut = (array[slice] < 10e10) & ( array[slice] != NaN )
#   cut = ( array[slice] > 3.0*array[slice][nancut].std() ) 
#   array[slice][cut] = 3.0*array[slice][nancut].std()
#   cut = ( array[slice] < -3.0*array[slice][nancut].std() ) 
#   array[slice][cut] = -3.0*array[slice][nancut].std()

#   Need to rotate array[slice] because axes were flipped
#   new_array = scipy.transpose(array[slice])
#   medianv = median(array[slice][nancut])

#   Alternate plotting command to set temperature limits
#  pylab.imshow(new_array, cmap='hot', vmin=-0.6, vmax=0.6, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.imshow(new_array, interpolation='gaussian', cmap='hot', extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.colorbar() #For some reason this isn't working, fixed...
#   pylab.savefig(filename2+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
#   pylab.clf()


