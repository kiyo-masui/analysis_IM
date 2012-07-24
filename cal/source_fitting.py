"""
Procedure to fit a power law curve to source spectroscopy data for sources such as 3C286, 3C048, 3C067, etc."""
import os

import scipy.optimize as opt 
import numpy as np
#import scipy as sp
#import numpy.ma as ma
from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT
import pylab

def function(freq,params):
   base = 750./freq
   function = params[0]*np.power(base,params[1])
   return function

def residuals(params,freq,data,errors):
   err = (data-function(freq,params))/errors
   return err

source = "3C048" # Or whatever source you want to use.
data_filename = source+"_source_data.txt"
data = np.loadtxt(data_filename)
freq = data[:,0]
#print freq
signal = data[:,1]
errors = data[:,2]
params0 = [19.6,0.495]
plsq = opt.leastsq(residuals,params0,args=(freq,signal,errors),full_output=1,maxfev=5000)
print plsq[0]
print plsq[1]
possibles = np.arange(freq[-1],freq[0],10)
#possibles = np.arange(freq[0],freq[-1],10)
m = [0.0,0.0,0.0,0.0,0.0]
merr = [0.0,0.0,0.0,0.0,0.0]
spectral_model = lambda m,x: m[0]*10**(m[1]*np.log10(x/150.0))*10**(m[2]*np.log10(x/150.0)**2)*10**(m[3]*np.log10(x/150.0)**3)*10**(m[4]*np.log10(x/150.0)**4)
if source=="3C048":
    m = [64.768,-0.387,-0.420,0.181,0.0]
    merr = [1.761,0.039,0.031,0.060,0.0]
if source=="3C286":
    m = [27.477,-0.158,0.032,-0.180,0.0]
    merr = [0.746,0.033,0.043,0.052,0.0]
if source=="3C295":
    m = [97.763,-0.582,-0.298,0.583,-0.363]
    merr = [2.787,0.045,0.085,0.116,0.137]

pylab.errorbar(freq,signal,errors,color='b',label='NED data')
pylab.plot(possibles,function(possibles,plsq[0]),color='g',label='Power law fit')
pylab.errorbar(possibles,spectral_model(m,possibles),spectral_model(merr,possibles),color = 'r',label='Spectral_model')
#pylab.plot(possibles,spectral_model(m,possibles),color='r')
#for i in range(0,len(m)):
#    m[i] = m[i]+merr[i]
#pylab.plot(possibles,spectral_model(m,possibles),color='m')
#pylab.plot(possibles,spectral_model(m-merr,possibles),color='m')
#pylab.ylim(15,25)
pylab.ylim(20,30)
pylab.xlim(700,900)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Jy)')
pylab.legend()
pylab.savefig(source+"_flux_fit.png")
pylab.clf()


   
   
