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
hit_plot = True
candidate_plot = True

#filename1 = sys.argv[1]

#0ld data
#file_base = 'fir_15hr_absorb_clean_map_I_'
#file_sets = ('712','732','749','751','768','771','788','790','807','810','827')

#new data
file_base = 'fir_15hr_41-80_ptcorr_absorb_clean_map_I_'
#file_sets = ('732','751','771','790','810') #A few frequency pixels are messed up so screwing with the std of the data.  
file_sets = ('732','751','771','790','810','829','849')

full_data = []
full_freqs = []
for i in range(0,len(file_sets)):
    array = al.load(directory+file_base+file_sets[i]+'.npy')
    vect_array = al.make_vect(array)
    #  creates a 3D array with indices (freq, ra, dec)
    ras = vect_array.get_axis('ra')
    decs = vect_array.get_axis('dec')
    freqs = vect_array.get_axis('freq')
    freqs = freqs/1e6
    full_data.append(array)
    full_freqs.append(freqs)
#    print freqs[0],freqs[-1]

print shape(full_data)
print shape(full_freqs)

old_freq = zeros((len(full_data)*len(full_data[0])))
old_array = zeros((len(full_data)*len(full_data[0]),16,10))
for i in range(0,len(full_data)):
    for j in range(0,len(full_data[0])):
        old_freq[400*i+j] = full_freqs[i][399-j]
        old_array[400*i+j,:,:] = full_data[i][399-j]
#        list_freq[i,j] = full_freqs[i][399-j]
#        list_array[i,j,:,:] = full_data[i][399-j]

freqsort = argsort(old_freq)
list_freq = zeros((len(old_freq)))
list_array = zeros((len(old_freq),16,10))
for i in range(0,len(freqsort)):
    list_freq[i] = old_freq[freqsort[i]]
    list_array[i,:,:] = old_array[freqsort[i]]
   

if fitting:
#   fits = zeros((len(list_freq),2,16,10))
    fits = zeros((2,16,10))
    for i in range(0,16):
        for j in range(0,10):
            data = list_array[:,i,j]
            params0 = [-0.5,2.5]
            plsq = opt.leastsq(residuals,params0,args=(list_freq,data),full_output=0,maxfev=5000)
            fits[0,i,j] = plsq[0][0]
            fits[1,i,j] = plsq[0][1]
#	    for k in range(0,len(list_freq)):
#                data = list_array[k,:,i,j]
#                params0 = [-0.5,2.5]
#                plsq = opt.leastsq(residuals,params0,args=(list_freq[k],data),full_output=0,maxfev=5000)
#        print plsq
#                fits[k,0,i,j] = plsq[0][0]
#                fits[k,1,i,j] = plsq[0][1]

#    fit_data = zeros((len(list_array),len(list_array[0]),16,10))
    fit_data = zeros((len(list_array),16,10))
    for i in range(0,16):
        for j in range(0,10):
            fit_data[:,i,j] = function(list_freq,fits[:,i,j])
#            for k in range(0,len(list_array)):
#                fit_data[k,:,i,j] = function(list_freq[k],fits[k,:,i,j])

if plotting:
    plt.plot(list_freq,list_array[:,8,5])
    plt.plot(list_freq,fit_data[:,8,5])
#    for i in range(0,len(list_freq)):
#        plt.plot(list_freq[i],list_array[i,:,8,5])
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

#    for i in range(0,len(list_freq)):
#        plt.plot(list_freq[i],fit_data[i,:,8,5])
    #plt.show()

    plt.clf()
#plt.title('High Resolution Residual Spectra for (15 arcmin)^2 Pixels') 
    for j in range(0,16):
        for k in range(0,10):
            index = j*10+k
#            index = (k-1)*12+(j-2)
#            index = (j-2)*8+(k-1)
            plt.subplot(16,10,index)
            plt.plot(list_freq,list_array[:,j,k]-fit_data[:,j,k])
#            for i in range(0,len(list_freq)):
#                plt.plot(list_freq[i],list_array[i,:,j,k]-fit_data[i,:,j,k])
            plt.xlim(700,780)
            plt.ylim(-0.05,0.05)
            plt.grid()
    plt.subplot(16,10,6)
    plt.title('High Resolution Power Law Subtracted (15 arcmin)^2 Spectra')
    plt.savefig('absorber_power_law_sub_maps',dpi=300)
    plt.show()
    plt.clf()

    for j in range(0,16):
        for k in range(0,10):
#            for i in range(0,len(list_freq)):
#                plt.plot(list_freq[i],list_array[i,:,j,k]-fit_data[i,:,j,k])
            plt.plot(list_freq,list_array[:,j,k]-fit_data[:,j,k])
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
#   full_filter = zeros((len(list_array),len(list_array[0]),16,10))
    full_filter = zeros((len(list_array),16,10))
#    for f in range(0,len(list_array)):
    for j in range(0,16):
        for k in range(0,10):
            for i in range(5,len(list_freq)-5):
                filter = list_array[i,j,k]-0.1*(list_array[i-5,j,k]+list_array[i-4,j,k]+list_array[i-3,j,k]+list_array[i-2,j,k]+list_array[i-1,j,k]+list_array[i+1,j,k]+list_array[i+2,j,k]+list_array[i+3,j,k]+list_array[i+4,j,k]+list_array[i+5,j,k])
                full_filter[i,j,k] = filter
    mean_filter = ma.mean(full_filter, axis=0)
    std_filter = ma.std(full_filter,axis=0)
#    print mean_filter[0]
#    print std_filter[0]
    outliers = zeros((len(list_array),16,10))
    outliers2 = zeros((len(list_array),16,10))
    outliers3 = zeros((len(list_array),16,10))
    posout = zeros((len(list_array),16,10))
    posout2 = zeros((len(list_array),16,10))
    posout3 = zeros((len(list_array),16,10))
#    for f in range(0,len(list_array)):
    for j in range(0,16):
        for k in range(0,10):
            for i in range(0,len(list_freq)):
                if full_filter[i,j,k]<(mean_filter[j,k]-1.5*std_filter[j,k]):
                    outliers[i,j,k] = 1.0
                if full_filter[i,j,k]<(mean_filter[j,k]-2.*std_filter[j,k]):
                    outliers2[i,j,k] = 1.0
                if full_filter[i,j,k]<(mean_filter[j,k]-3.*std_filter[j,k]):
                    outliers3[i,j,k] = 1.0
                if full_filter[i,j,k]>(mean_filter[j,k]+1.5*std_filter[j,k]):
                    posout[i,j,k] = 1.0
                if full_filter[i,j,k]>(mean_filter[j,k]+2.*std_filter[j,k]):
                    posout2[i,j,k] = 1.0
                if full_filter[i,j,k]>(mean_filter[j,k]+3.*std_filter[j,k]):
                    posout3[i,j,k] = 1.0
#            print 'Number of neg 3sig outliers for pixel ', j,k,' is:', ma.sum(outliers3[:,j,k]) 
#            print 'Number of pos 3sig outliers for pixel ', j,k,' is:', ma.sum(posout3[:,j,k])
#            print 'Number of neg 2sig outliers for pixel ', j,k,' is:', ma.sum(outliers2[:,j,k])
#            print 'Number of pos 2sig outliers for pixel ', j,k,' is:', ma.sum(posout2[:,j,k])

    if hit_plot:
#        for i in range(0,len(list_freq[0])):
        for f in range(0,len(list_freq)):
            if ma.sum(ma.sum(outliers3[f,:,:],axis=0))>=1.0:
                print list_freq[f]*1000
                plt.imshow(outliers3[f])
                plt.title('Three Sigma Outliers for an 11 pixel filter')
                plt.savefig('three_sigma_outliers_'+str(int(list_freq[f]*1000))+'kHz',dpi=300)
                plt.clf()
    
    candidate_freqs = []
    candidate_ra = []
    candidate_dec = []
    for f in range(0,len(list_array)):
#        for i in range(0,len(list_freq[0])):
        med_sum1 = ma.sum(outliers[f,:,:],axis=0)
        sum1 = ma.sum(med_sum1)
        med_sum2 = ma.sum(outliers2[f,:,:],axis=0)
        sum2 = ma.sum(med_sum2)
        med_sum3 = ma.sum(outliers3[f,:,:],axis=0)
        sum3 = ma.sum(med_sum3)
        if sum1>=1.0:
            if sum1<=2.0:
#                candidate_maps.append(f)
                candidate_freqs.append(f)
                print where(outliers[f]==1.0)
                candidate_ra.append(where(outliers[f]==1.0)[0][0])
                candidate_dec.append(where(outliers[f]==1.0)[1][0])
#                if sum2==2.0:
#                    candidate_maps.append(f)
#                    candidate_freqs.append(f)
#                    candidate_ra.append(where(outliers2[f]==1.0)[0][1])
#                    candidate_dec.append(where(outliers2[f]==1.0)[1][1])
                print 'Num of Neg outliers at ',int(list_freq[f]*1000),' kHz is:',sum2
                print 'Num of Pos outliers at ',int(list_freq[f]*1000),' kHz is:',ma.sum(posout2[f,:,:])
    
    print len(candidate_freqs)

    if candidate_plot:
        for g in range(0,len(candidate_freqs)):
            min_freq = candidate_freqs[g]-30
            max_freq = candidate_freqs[g]+30
            if candidate_freqs[g]<30:
                min_freq = 0
            elif candidate_freqs[g]>(len(list_freq)-30):
                max_freq = -1
            mean_signal = ma.mean(list_array[min_freq:max_freq,candidate_ra[g],candidate_dec[g]])
            plt.plot(1000*list_freq[min_freq:max_freq],list_array[min_freq:max_freq,candidate_ra[g],candidate_dec[g]]-mean_signal)
            plt.xlabel('Frequency (kHz)')
            plt.axvline(1000*list_freq[candidate_freqs[g]],c='r',ls='--')
            plt.ylabel('Temperature (Kelvin)')
            plt.grid()
            plt.title('Mean Subtracted Absorber Candidate Spectra at RA %0.1f, DEC %0.1f' %(ras[candidate_ra[g]],decs[candidate_dec[g]]))
            plt.savefig('absorber_candidate_at_'+str(int(list_freq[candidate_freqs[g]]*1000))+'kHz',dpi=300)
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


