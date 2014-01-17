#This code is designed to generate a higher resolution version of a calibration file using a linear fit between points.

import pylab
import sys
from numpy import *
import scipy as sp
import numpy.ma as ma

filedir = sys.argv[1]
fname = '15hr_fdg_calc_avg_new'
fend = '.txt'

old_file = loadtxt(filedir+fname+fend)

old_file = array(old_file)

bin_size = 0.01 #Values in MHz

new_file = sp.zeros((int(200/bin_size),3))

new_file[0,0] = 900. #frequencies are in descending order
for i in range(1,len(new_file)):
    new_file[i,0]=new_file[i-1,0]-bin_size

for i in range(0,len(new_file)):
    if new_file[i,0]>=old_file[0,0]:
        new_file[i,1] = old_file[0,1]
        new_file[i,2] = old_file[0,2]
    elif new_file[i,0]<=old_file[-1,0]:
        new_file[i,1] = old_file[-1,1]
        new_file[i,2] = old_file[-1,2]
    else:
        for j in range(0,len(old_file)):
            if new_file[i,0]==old_file[j,0]:
                new_file[i,1] = old_file[j,1]
                new_file[i,2] = old_file[j,2]
            elif new_file[i,0]<old_file[j,0]:
                large = j

        mx = (old_file[large,1]-old_file[large+1,1])/(old_file[large,0]-old_file[large+1,0])
        bx = old_file[large,1]-mx*old_file[large,0]
        my = (old_file[large,2]-old_file[large+1,2])/(old_file[large,0]-old_file[large+1,0])
        by = old_file[large,2]-my*old_file[large,0]
        new_file[i,1] = mx*new_file[i,0]+bx
        new_file[i,2] = my*new_file[i,0]+by

#print 'x-means', ma.mean(new_file[:,1]), ma.mean(old_file[:,1])
#print 'y-means', ma.mean(new_file[:,2]), ma.mean(old_file[:,2])

pylab.scatter(old_file[:,0],old_file[:,1],c='b',edgecolor='b',s=3,label='X-low')
pylab.scatter(old_file[:,0],old_file[:,2],c='g',edgecolor='g',s=3,label='Y-low')
pylab.scatter(new_file[:,0],new_file[:,1],c='c',edgecolor='c',s=1,label='X-high')
pylab.scatter(new_file[:,0],new_file[:,2],c='r',edgecolor='r',s=1,label='Y-high')
pylab.legend()
pylab.xlim(700.0,900.0)
pylab.xlabel('Frequency (MHz)')
pylab.ylim(0.5,3.0)
pylab.ylabel('Correction Factor')
pylab.savefig('15hr_high_res_convert',dpi=200)

savetxt(filedir+fname+'high_res'+fend, new_file,delimiter=' ')
