# This script is set to make 3D plots of the XX and YY gain corrections oover a range of sessions.

import pylab
import sys
from numpy import *
import scipy as sp
import scipy.interpolate as ip

prefix = ''

#This is $GBT_OUT/diff_gain_params/
filedir = sys.argv[1]
suffix = '_diff_gain_calc.txt'

#Currently set for the range of all 10B_036 with guppi
#May want to rebuild depending what sessions you want to grab.
directory = 'GBT10B_036'
directory2 = 'GBT11B_055'
directory3 = 'GBT12A_418'
min_sess = 41
min_sess2 = 01
min_sess3 = 01
max_sess = 90
max_sess2 = 19
max_sess3 = 03
#Number of sessions used
len = max_sess - min_sess
len2 = max_sess2-min_sess2
len3 = max_sess3-min_sess3
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
for i in range(0,len3):
   label = min_sess3+i
   if label<10:
      label ='0'+str(label)
   else:
      label = str(label)
   filenames.append(filedir+directory3+'/'+label+suffix)
   sessions.append(label)

#print sessions
#print filenames

#Need a test_dir to determine the number of frequency slices
array =[] 
test_dir = loadtxt(filenames[0])
scale = shape(test_dir)
#print array

#Database of mueller files
tot_len = len+len2+len3
for j in range(0,tot_len):
   source = filenames[j]
   data = loadtxt(source)
   array.append(data)
#print array[0][0][1]
#print array[1][0][1]

#List of frequencies for axis.
freqs = test_dir[:,0]
freqs_new = arange(freqs[-1],freqs[0]+0.78125,0.78125)

#print freqs

XG = sp.zeros((tot_len,scale[0]))
YG = sp.zeros((tot_len,scale[0]))
XG_mod = sp.zeros((tot_len,scale[0]))
YG_mod = sp.zeros((tot_len,scale[0]))
XG_comp = sp.zeros((tot_len,scale[0]))
YG_comp = sp.zeros((tot_len,scale[0]))

#Building separate tables for each param.
for i in range(0,tot_len):
   for j in range(0,scale[0]):
      XG[i,j] = array[i][j][1]
      YG[i,j] = array[i][j][2]
      XG_mod[i,j] = array[i][j][1]
      YG_mod[i,j] = array[i][j][2]
   for k in range(102,110):
      XG_mod[i,k] = 0.5*(XG_mod[i,101]+XG_mod[i,111])
      YG_mod[i,k] = 0.5*(YG_mod[i,101]+YG_mod[i,111])
   for k in range(130,137):
      XG_mod[i,k] = 0.5*(XG_mod[i,129]+XG_mod[i,138])
      YG_mod[i,k] = 0.5*(YG_mod[i,129]+YG_mod[i,138])
   for k in range(1,scale[0]-1):
      if abs(XG_mod[i,k]-XG_mod[i,k-1])>0.5:
         XG_mod[i,k] = XG_mod[i,k-1]
      if abs(YG_mod[i,k]-YG_mod[i,k-1])>0.5:
         YG_mod[i,k] = YG_mod[i,k-1]
   Xmod = ip.UnivariateSpline(freqs_new,XG_mod[i,:],s=1,k=1)
   Ymod = ip.UnivariateSpline(freqs_new,YG_mod[i,:],s=1,k=2)
   XG_mod[i] = Xmod(freqs_new)
   YG_mod[i] = Ymod(freqs_new)
   
for i in range(0,tot_len):
   for j in range(0,scale[0]):
      XG_comp[i,j] = abs(XG[i,j]-XG[67,j])/XG[67,j]
      YG_comp[i,j] = abs(YG[i,j]-YG[67,j])/YG[67,j]

#Going to make the two gains as subplots
pylab.suptitle('XX and YY Gain Correction Factors (Scale is 3.25-0.75)')
pylab.subplot(211)
#pylab.imshow(m_II/m_II,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.,vmax=1.,extent=(freqs.min(),freqs.max(),min_sess,max_sess))
pylab.imshow(XG, cmap='hot', origin='lower',vmin=0.75,vmax=3.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal')
#pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Session')
#pylab.savefig('XX_Gain_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf()

pylab.subplot(212)
pylab.imshow(YG, cmap='hot', origin='lower',vmin=0.75,vmax=3.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Session')
pylab.savefig('XX_and_YY_gain_'+str(min_sess)+'-'+str(max_sess-1)+'_and_'+str(min_sess2)+'-'+str(max_sess2-1)+'_and_'+str(min_sess3)+ '-'+str(max_sess3-1)+'_compare.png')
pylab.clf()

pylab.subplot(211)
pylab.imshow(XG_mod,cmap='hot',origin='lower',vmin=0.75,vmax=3.25,extent=(freqs_new.max(),freqs_new.min(),0,tot_len))
pylab.colorbar(orientation='horizontal')

pylab.subplot(212)
pylab.imshow(YG_mod,cmap='hot',origin='lower',vmin=0.75,vmax=3.25,extent=(freqs_new.max(),freqs_new.min(),0,tot_len))
pylab.colorbar(orientation='horizontal')
pylab.savefig('XX_and_YY_gain_smoothed_'+str(min_sess)+'-'+str(max_sess2-1)+'_compare.png')
pylab.clf()

# Plotting deviation from a reference correction factor 
pylab.suptitle('XX and YY Gain Correction Factor Deviation (deltaG/G, max = 15%) ')
pylab.subplot(211)
#pylab.imshow(m_II/m_II,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.,vmax=1.,extent=(freqs.min(),freqs.max(),min_sess,max_sess))
pylab.imshow(XG_comp, cmap='hot', origin='lower',vmin=0,vmax=0.15,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal')
#pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Session')
#pylab.savefig('XX_Gain_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf()

pylab.subplot(212)
pylab.imshow(YG_comp, cmap='hot', origin='lower',vmin=0,vmax=0.15,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Session')
pylab.savefig('XX_and_YY_gain_'+str(min_sess)+'-'+str(max_sess-1)+'_and_'+str(min_sess2)+'-'+str(max_sess2-1)+'_and_'+str(min_sess3)+ '-'+str(max_sess3-1)+'_correction_variation_from_2(418).png')
pylab.clf()

