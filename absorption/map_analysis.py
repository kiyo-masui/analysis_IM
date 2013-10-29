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
fitting = True
plotting = False
filtered = True
hit_plot = False
candidate_plot = True

#filename1 = sys.argv[1]

file_base = 'fir_15hr_absorb_clean_map_I_'

file_sets = ('712','732','749','751','768','771','788','790','807','810','827')

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
#    print freqs[0],freqs[-1]

print shape(full_data)
print shape(full_freqs)

list_freq = zeros((len(full_data),len(full_data[0])))
list_array = zeros((len(full_data),len(full_data[0]),16,10))
for i in range(0,len(full_data)):
    for j in range(0,len(full_data[0])):
        list_freq[i,j] = full_freqs[i][399-j]
        list_array[i,j,:,:] = full_data[i][399-j]

if fitting:
    fits = zeros((len(list_freq),2,16,10))
    for i in range(0,16):
        for j in range(0,10):
	    for k in range(0,len(list_freq)):
                data = list_array[k,:,i,j]
                params0 = [-0.5,2.5]
                plsq = opt.leastsq(residuals,params0,args=(list_freq[k],data),full_output=0,maxfev=5000)
#        print plsq
                fits[k,0,i,j] = plsq[0][0]
                fits[k,1,i,j] = plsq[0][1]

    fit_data = zeros((len(list_array),len(list_array[0]),16,10))
    for i in range(0,16):
        for j in range(0,10):
            for k in range(0,len(list_array)):
                fit_data[k,:,i,j] = function(list_freq[k],fits[k,:,i,j])

if plotting:
    for i in range(0,len(list_freq)):
        plt.plot(list_freq[i],list_array[i,:,8,5])
    #plt.show()

#plt.subplot(121)
#plt.imshow(fits[0,0],interpolation='nearest')
#plt.colorbar()
#plt.title('Amplitude')
#plt.subplot(122)
#plt.imshow(fits[0,1],interpolation='nearest')
#plt.colorbar()
#plt.title('Power')
#plt.show()
#print fits
#plt.clf()

    for i in range(0,len(list_freq)):
        plt.plot(list_freq[i],fit_data[i,:,8,5])
    #plt.show()

    plt.clf()
#plt.title('High Resolution Residual Spectra for (15 arcmin)^2 Pixels') 
    for j in range(0,16):
        for k in range(0,10):
            index = j*10+k
#            index = (k-1)*12+(j-2)
#            index = (j-2)*8+(k-1)
            plt.subplot(16,10,index)
            for i in range(0,len(list_freq)):
                plt.plot(list_freq[i],list_array[i,:,j,k]-fit_data[i,:,j,k])
            plt.xlim(700,780)
            plt.ylim(-0.05,0.05)
            plt.grid()
    plt.subplot(16,10,6)
    plt.title('High Resolution Power Law Subtracted (15 arcmin)^2 Spectra')
    plt.savefig('absorber_power_law_sub_maps',dpi=300)
#    plt.show()
    plt.clf()

    for j in range(0,16):
        for k in range(0,10):
            for i in range(0,len(list_freq)):
                plt.plot(list_freq[i],list_array[i,:,j,k]-fit_data[i,:,j,k])
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Single Spatial Pixel Spectrum Residuals (Kelvin)')
            plt.ylim(-0.1,0.1)
            plt.grid()
            plt.title('Spectra for (15 arcmin)^2 Pixel')
            plt.savefig('absorber_power_law_single_map_RA'+str(j)+'_DEC'+str(k),dpi=300)
            plt.clf()


if filtered:
    #Start with a single frequency band and 11 pixel filter
#    for i in range(5,len(list_freq[0])-5):
    full_filter = zeros((len(list_array),len(list_array[0]),16,10))
    for f in range(0,len(list_array)):
        for j in range(0,16):
            for k in range(0,10):
                for i in range(5,len(list_freq[0])-5):
                    filter = list_array[f,i,j,k]-0.1*(list_array[f,i-5,j,k]+list_array[f,i-4,j,k]+list_array[f,i-3,j,k]+list_array[f,i-2,j,k]+list_array[f,i-1,j,k]+list_array[f,i+1,j,k]+list_array[f,i+2,j,k]+list_array[f,i+3,j,k]+list_array[f,i+4,j,k]+list_array[f,i+5,j,k])
                    full_filter[f,i,j,k] = filter
    mean_filter = ma.mean(full_filter, axis=1)
    std_filter = ma.std(full_filter,axis=1)
#    print mean_filter[0]
#    print std_filter[0]
    outliers = zeros((len(list_array),len(list_array[0]),16,10))
    outliers2 = zeros((len(list_array),len(list_array[0]),16,10))
    outliers3 = zeros((len(list_array),len(list_array[0]),16,10))
    for f in range(0,len(list_array)):
        for j in range(0,16):
            for k in range(0,10):
                for i in range(0,len(list_freq[0])):
                    if full_filter[f,i,j,k]<(mean_filter[f,j,k]-std_filter[f,j,k]):
                        outliers[f,i,j,k] = 1.0
                    if full_filter[f,i,j,k]<(mean_filter[f,j,k]-2.*std_filter[f,j,k]):
                        outliers2[f,i,j,k] = 1.0
                    if full_filter[f,i,j,k]<(mean_filter[f,j,k]-3.*std_filter[f,j,k]):
                        outliers3[f,i,j,k] = 1.0

    if hit_plot:
        for i in range(0,len(list_freq[0])):
            for f in range(0,len(list_freq)):
                plt.imshow(outliers[f,i])
                plt.title('One Sigma Outliers for an 11 pixel filter')
                plt.savefig('one_sigma_outliers_'+str(int(list_freq[f,i]*1000))+'kHz',dpi=300)
                plt.clf()
    
    candidate_maps = []
    candidate_freqs = []
    candidate_ra = []
    candidate_dec = []
    for f in range(0,len(list_array)):
        for i in range(0,len(list_freq[0])):
            med_sum2 = ma.sum(outliers2[f,i,:,:],axis=0)
            sum2 = ma.sum(med_sum2)
            med_sum3 = ma.sum(outliers3[f,i,:,:],axis=0)
            sum3 = ma.sum(med_sum3)
            if sum3>=1.0:
                if sum2<=2.0:
                    candidate_maps.append(f)
                    candidate_freqs.append(i)
                    print where(outliers3[f,i]==1.0)
                    candidate_ra.append(where(outliers3[f,i]==1.0)[0][0])
                    candidate_dec.append(where(outliers3[f,i]==1.0)[1][0])
                    if sum3==2.0:
                        candidate_maps.append(f)
                        candidate_freqs.append(i)
                        candidate_ra.append(where(outliers3[f,i]==1.0)[0][1])
                        candidate_dec.append(where(outliers3[f,i]==1.0)[1][1])
                        
    print len(candidate_maps)

    if candidate_plot:
        for g in range(0,len(candidate_maps)):
            min_freq = candidate_freqs[g]-30
            max_freq = candidate_freqs[g]+30
            if candidate_freqs[g]<30:
                min_freq = 0
            elif candidate_freqs[g]>378:
                max_freq = -1
            mean_signal = ma.mean(list_array[candidate_maps[g],min_freq:max_freq,candidate_ra[g],candidate_dec[g]])
            plt.plot(1000*list_freq[candidate_maps[g],min_freq:max_freq],list_array[candidate_maps[g],min_freq:max_freq,candidate_ra[g],candidate_dec[g]]-mean_signal)
            plt.xlabel('Frequency (kHz)')
            plt.axvline(1000*list_freq[candidate_maps[g],candidate_freqs[g]],c='r',ls='--')
            plt.ylabel('Temperature (Kelvin)')
            plt.grid()
            plt.title('Mean Subtracted Absorber Candidate Spectra at RA %0.1f, DEC %0.1f' %(ras[candidate_ra[g]],decs[candidate_dec[g]]))
            plt.savefig('absorber_candidate_at_'+str(int(list_freq[candidate_maps[g],candidate_freqs[g]]*1000))+'kHz',dpi=300)
            plt.clf()

#plt.subplot(121)
#plt.imshow(list_array[:,:,8,5],aspect=len(list_array[0])/len(list_array),vmin=-0.6,vmax=0.6,interpolation='nearest')
#plt.colorbar()
#plt.title('Before Subtraction')
#plt.subplot(122)
#plt.imshow(list_array[:,:,8,5]-fit_data[:,:,8,5],aspect=len(list_array[0])/len(list_array),vmin=-0.1,vmax=0.1,interpolation='nearest')
#plt.title('After Subtraction')
#plt.colorbar()
#plt.show()

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


