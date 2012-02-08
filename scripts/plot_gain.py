import pylab
from numpy import *
import scipy as sp

prefix = '41'

gain_params = loadtxt('diff_gain_calc.txt')

size = len(gain_params[:,0])
#print mueller_params[0,:] #is a set of values for highest frequency
freq = gain_params[:,0]
XG = gain_params[:,1]
YG = gain_params[:,2]

#Gain plot
pylab.plot(freq, XG,label='XX Gain')
pylab.plot(freq, YG,label='YY Gain')
pylab.ylim(0,3)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('diff_gain_'+prefix+'.png')
pylab.clf()

