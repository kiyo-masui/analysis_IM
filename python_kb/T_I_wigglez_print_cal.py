import pyfits
import pylab
import ephem
import time
from PIL import Image
from numpy import *
import gbt_cal

cal_filename = '/home/scratch/kbandura/GBT10B_036/02_3c48_onoff_161-164.raw.acs.fits'

T_noise,T_sys =  gbt_cal.getcal_3c48(cal_filename)

print T_noise.shape
print T_sys.shape

print "finished loading cal"

pylab.plot(T_noise[3])
pylab.xlabel('freq')
pylab.ylabel('temperature')
pylab.ylim(0,10)
pylab.savefig('wigglez_02_3c48_noise_cal_temp.png', dpi=300)

pylab.clf()
 
pylab.plot(T_sys[3])
pylab.xlabel('freq')
pylab.ylabel('temperature')
pylab.ylim(0,100)
pylab.savefig('wigglez_02_3c48_sys_temp.png', dpi=300)

