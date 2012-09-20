#This script is designed to generate an average calibration correction table (for a flux/differential gain style calibration) using a set of data files such as all files corresponding to a given field.
#If 1hr field: 10B_036 86-89 and 11B_055 01-18 (sesssion 14 is trash so that calibration file is set to match session 13). 
#If 15hr field: 10B_036 41-80

import pylab
import sys
from numpy import *
import scipy as sp

prefix = ''
filedir = sys.argv[1]
suffix = '_diff_gain_calc.txt'
#directory = 'GBT10B_036'
directory2 = 'GBT11B_055'
directory = 'GBT12A_418'
min_sess = 1
max_sess = 23
#min_sess = 41
min_sess2 = 0
#max_sess = 80
max_sess2 = 0
len = max_sess-min_sess
len2 = max_sess2-min_sess2
sessions = []
filenames = []
for i in range(0,len):
    label = min_sess+i
    if label<10:
        label = '0'+str(label)
    else:
        label = str(label)
    filenames.append(filedir+directory+'/'+label+suffix)
    sessions.append(label)
#for i in range(0,len2):
#    label = min_sess2+i
#    if label<10:
#        label = '0'+str(label)
#    else:
#        label = str(label)
#    filenames.append(filedir+directory2+'/'+label+suffix)
#    sessions.append(label)

array = []
test_dir = loadtxt(filenames[0])
scale = shape(test_dir)
tot_len = len+len2
for j in range(0,tot_len):
    source = filenames[j]
    data = loadtxt(source)
    array.append(data)
#print array
freqs = test_dir[:,0]
XG = sp.zeros((tot_len,scale[0]))
YG = sp.zeros((tot_len,scale[0]))
for i in range(0,tot_len):
    for j in range(0,scale[0]):
        XG[i,j] = array[i][j][1]
        YG[i,j] = array[i][j][2]


XG_avg = mean(XG,axis=0)         
#print XG_avg
YG_avg = mean(YG,axis=0)
#print YG_avg

pylab.plot(freqs,XG_avg,'b-',label='Average XX Gain')
pylab.plot(freqs,YG_avg,'g-',label='Average YY Gain')
pylab.xlim(freqs[-1],freqs[0])
pylab.xlabel('frequency (MHz)')
pylab.ylabel('Correction Factor')
pylab.legend()
pylab.savefig('11hr_avg_fdg_correction.png')
pylab.clf()

#extra_file = str(22)+suffix
#extra_data = loadtxt(extra_file)
#print shape(extra_data)
#YG_extra = extra_data[:,2]

for i in range(1,len):
    pylab.plot(freqs,XG[i],'b-',label='Sess %2i XX Gain' %i)
    pylab.plot(freqs,YG[i],'g-',label='Sess %2i YY Gain' %i)

#pylab.plot(freqs,YG_extra,'r-')
pylab.xlim(freqs[-1],freqs[0])
pylab.ylim(0,5)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Correction Factor')
#pylab.legend()
pylab.savefig('11hr_fdg_corrections.png')
pylab.clf()


output = sp.zeros((scale[0],3))
for f in range(0,scale[0]):
    output[f,0] = freqs[f]
    output[f,1] = XG_avg[f]
    output[f,2] = YG_avg[f]

savetxt('11hr_fdg_calc_avg.txt',output,delimiter=' ')
