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

def function(freq,params):
   base = 750./freq
   function = params[0]*np.power(base,params[1])
   return function

def residuals(params,freq,data):
   err = data-function(freq,params)
   return err

source = "3C048" # Or whatever source you want to use.
data_filename = source+"_source_data.txt"
data = np.loadtxt(data_filename)
freq = data[:,0]
#print freq
signal = data[:,1]
params0 = [19.6,0.495]
plsq = opt.leastsq(residuals,params0,args=(freq,signal),full_output=1,maxfev=5000)
print plsq[0]
print plsq[1]



   
   
