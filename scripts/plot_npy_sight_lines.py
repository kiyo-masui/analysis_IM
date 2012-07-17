#This script generates a set of maps in RA and DEC for each frequency in the map data. Can be used to make plots of maps.

import pylab
import pyfits
import sys
from numpy import *
#from core import algebra as al
#from utils import * 
import core.algebra as al
import scipy
import scipy.optimize as opt
import scipy.stats as stat
#from pylab import *

def function(freq,params):
   base = 750./freq
   function = params[0]*power(base,params[1])
   return function

def residuals(params,freq,signal):
   err = signal-function(freq,params)
   return err

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

for slice, dec in enumerate(decs):
    for line, ra in enumerate(ras):
        # Now have "slice" being the dec index, "line" being the ra index
#        pylab.plot(freqs,array[:,line,slice])
#        pylab.ylim(-0.5,0.5)
#        pylab.savefig(filename2+'_'+str(ra)+'_'+str(dec)+'_.png')
#        pylab.clf()
        data = array[:,line,slice]
        bad_pts = logical_or(isnan(data),isinf(data))
        good_pts = where(logical_not(bad_pts))
        data = data[good_pts]
#        mean_sig = mean(data)
#        print mean_sig
#        data = data/mean_sig
#        print data
        freqs_good = freqs[good_pts]
        params0=[19.6,0.495]
        plsq = opt.leastsq(residuals,params0,args=(freqs_good,data),full_output=1,maxfev=5000)
#        print ra, dec, plsq[0]
        remainder = data-function(freqs_good,plsq[0])
#        pylab.plot(freqs_good,remainder)
#        pylab.ylim(-0.5,0.5)
#        pylab.savefig(filename2+'_dev_'+str(ra)+'_'+str(dec)+'.png')
#        pylab.clf()
        filter = remainder
        for f in range(5, len(data)-5):
#            filter[f] = filter[f]-filter[f-1]
#            filter[f] = filter[f]-0.5*(filter[f+1]+filter[f-1])
            filter[f] = filter[f]-0.1*(filter[f+1]+filter[f+2]+filter[f+3]+filter[f+4]+filter[f+5]+filter[f-1]+filter[f-2]+filter[f-3]+filter[f-4]+filter[f-5])
        filter[0] = 0
        filter[1] = 0
        filter[2] = 0
        filter[3] = 0
        filter[4] = 0
        filter[-1] = 0
        filter[-2] = 0
        filter[-3] = 0
        filter[-4] = 0
        filter[-5] = 0
#        print len(filter)
#        pylab.hist(filter,200,normed=1)
#        pylab.savefig(filename2+'_hist_'+str(ra)+'_'+str(dec)+'.png')
#        pylab.clf()
#        pylab.hist(filter,bins=200, range=(-1.0,-0.05))
#        pylab.savefig(filename2+'_hist_trunc_'+str(ra)+'_'+str(dec)+'.png')
#        pylab.clf()
#        mean = stat.gmean(filter)
#        var = stat.variation(filter)
        f_hist = stat.histogram(filter,350)
#        print f_hist
        scale = arange(f_hist[1],f_hist[1]+350*f_hist[2],f_hist[2])
#        print len(scale), len(f_hist[0])
        fitfunc = lambda p,x: p[0]*exp(-(p[1]-x)**2/(2*p[2]**2))
        errfunc = lambda p,x,y: fitfunc(p,x)-y
        pylab.plot(scale,f_hist[0])
#        a =1.
#        mean = float(mean)
#        var = float(var)
#        print a
#        p0 = [a,mean,var]
        p0 = [1000.0,0.00001,0.001]
#        print any(isnan(f_hist[0])),any(isinf(f_hist[0])), any(isnan(scale)), any(isinf(scale))
        fit = opt.leastsq(errfunc,p0[:],args=(scale,f_hist[0]),maxfev=500000)
#        print mean,var
#        print fit[0]
#        print '________________'
        pylab.plot(scale,fitfunc(fit[0],scale))
        pylab.savefig(filename2+'_gauss_fit_'+str(ra)+'_'+str(dec)+'.png')
        pylab.clf()
        pylab.plot(scale,abs(fitfunc(fit[0],scale)-f_hist[0]))
        pylab.savefig(filename2+'_gauss_resid_'+str(ra)+'_'+str(dec)+'.png')
        pylab.clf()
        





# Old code for npy plotting maps
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
#   pylab.imshow(new_array, cmap='hot', vmin=-0.1, vmax=0.1, extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.imshow(new_array, interpolation='gaussian', cmap='hot', extent=(ras.max(),ras.min(),decs.min(),decs.max()), origin='lower')
#   pylab.colorbar() #For some reason this isn't working, fixed...
#   pylab.savefig(filename2+str(freq)[:3]+'.png')
#   pylab.savefig('v_'+filename2+str(freq)[:3]+'.png')
#   pylab.clf()

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


