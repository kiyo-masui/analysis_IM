import sys
import pylab
from numpy import *
import scipy as sp
import scipy.interpolate as ip

prefix = '61'

filedir = sys.argv[1]

gain_params = loadtxt(filedir+prefix+'_diff_gain_calc_mod.txt')

size = len(gain_params[:,0])
#print mueller_params[0,:] #is a set of values for highest frequency
freq = gain_params[:,0]
#print len(freq)
freq_new = arange(freq[-1],freq[0]+0.78125,0.78125)
#print freq
#print len(freq)

XG = gain_params[:,1]
YG = gain_params[:,2]
#print XG[100]
#print XG[129]
#print YG[100]
#print YG[129]

for i in range(102,110):
    XG[i] = 0.5*(XG[101]+XG[111])
    YG[i] = 0.5*(YG[101]+YG[111])
for j in range(130,137):
    XG[j] = 0.5*(XG[129]+XG[138])
    YG[j] = 0.5*(YG[129]+YG[138])
#for i in range(1,len(XG)-1):
#    if abs(XG[i]-XG[i-1])>0.5:
#        XG[i] = XG[i-1]
#    if abs(YG[i]-YG[i-1])>0.5:
#        YG[i] = YG[i-1]

XG_mod = ip.UnivariateSpline(freq_new,XG,s=1,k=1)
YG_mod = ip.UnivariateSpline(freq_new,YG,s=1,k=2)
#print XG_mod(800.)
freq_new = arange(freq_new[0],freq_new[-1],.78125)
#print freq[0]
#print freq[-1]
#print freq_new
XG_mod_f = XG_mod(freq_new)
YG_mod_f = YG_mod(freq_new)
XG_new = sp.zeros(len(XG))
YG_new = sp.zeros(len(YG))
for i in range(0,len(XG)):
    XG_new[i] = XG[-i]
    YG_new[i] = YG[-i]
#print XG_mod_f

#Gain plot
pylab.plot(freq, XG_new,label='XX Gain')
pylab.plot(freq, YG_new,label='YY Gain')
pylab.plot(freq_new,XG_mod_f, label='XX_mod')
pylab.plot(freq_new,YG_mod_f, label='YY_mod')
pylab.ylim(0,4)
pylab.xlim(freq[0],freq[-1])
pylab.xlabel('frequency')
pylab.legend()
pylab.grid(b='on')
pylab.savefig('diff_gain_mod_'+prefix+'.png')
pylab.clf()

