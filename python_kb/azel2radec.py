import ephem
from numpy import *
import time
import pyfits

file = pyfits.open('/home/scratch/kbandura/GBT10B_036/00_wigglez1hr_azel_153-160.raw.acs.fits')
data = file[1].data
el = data.field('crval3')
az = data.field('crval2')
times = data.field('DATE-OBS')
dur = data.field('DURATION')

GBT = ephem.Observer()
GBT.long = '-79:50:23.4'
GBT.lat = '38:25:59.23'
GBT.pressure = 0
GBT.temp = 0
el_r = el*pi/180.0
az_r = az*pi/180.0

max_times = times.shape[0]
rag = zeros(max_times)
decg = zeros(max_times)


for i in range(0,max_times):
     t1, t2 = times[i].split(".",1)
     t21 = str(float("0."+t2)+dur[i]/2)
     t23, t22 = t21.split(".",1)
     t3 = time.strptime(t1, "%Y-%m-%dT%H:%M:%S")
     t4 = time.strftime("%Y/%m/%d %H:%M:%S", t3)
     GBT.date = t4 + "." + t22
     rag[i], decg[i] = GBT.radec_of(az_r[i],el_r[i])

rag_d = rag*180/pi
decg_d = decg*180/pi
