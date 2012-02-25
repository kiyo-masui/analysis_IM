import pylab
import sys
from numpy import *
import scipy as sp

prefix = ''

#This is $GBT_OUT/diff_gain_params/
filedir = sys.argv[1]
suffix = '_diff_gain_calc.txt'

#Currently set for the range of all 10B_036 with guppi
#May want to rebuild depending what sessions you want to grab.
directory = 'GBT10B_036'
directory2 = 'GBT11B_055'
min_sess = 41
min_sess2 = 01
max_sess = 89
max_sess2 = 18
#Number of sessions used
len = max_sess - min_sess
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
for i in range(0,len2):
   label = min_sess2+i
   if label<10:
      label = '0'+str(label)
   else:
      label = str(label)
   filenames.append(filedir+directory2+'/'+label+suffix)
   sessions.append(label)

#print filenames

#Need a test_dir to determine the number of frequency slices
array =[] 
test_dir = loadtxt(filenames[0])
scale = shape(test_dir)
#print array

#Database of mueller files
tot_len = len+len2
for j in range(0,tot_len):
   source = filenames[j]
   data = loadtxt(source)
   array.append(data)
#print array[0][0][1]
#print array[1][0][1]

#List of frequencies for axis.
freqs = test_dir[:,0]

#print freqs

XG = sp.zeros((tot_len,scale[0]))
YG = sp.zeros((tot_len,scale[0]))

#Building separate tables for each param.
for i in range(0,tot_len):
   for j in range(0,scale[0]):
      XG[i,j] = array[i][j][1]
      YG[i,j] = array[i][j][2]


#Going to make the two gains as subplots
pylab.subplot(211)
#pylab.imshow(m_II/m_II,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.,vmax=1.,extent=(freqs.min(),freqs.max(),min_sess,max_sess))
pylab.imshow(XG,interpolation='gaussian', cmap='hot', origin='lower',vmin=0.75,vmax=3.25,extent=(freqs.max(),freqs.min(),0,tot_len))
pylab.colorbar(orientation='horizontal')
#pylab.savefig('XX_Gain_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf()

pylab.subplot(212)
pylab.imshow(YG,interpolation='gaussian', cmap='hot', origin='lower',vmin=0.75,vmax=3.25,extent=(freqs.max(),freqs.min(),0,tot_len))
pylab.colorbar(orientation='horizontal')
pylab.savefig('XX_and_YY_gain_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
pylab.clf()

