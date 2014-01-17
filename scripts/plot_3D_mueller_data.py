# This script is set to plot the generated mueller matrix elements for a set of sessions. Currently set so that it plots all matrix elements other than mII as the element divided by mII so that flux is separated from polarization.
import pylab
import sys
from numpy import *
import scipy as sp

prefix = ''

#This is $GBT_OUT/mueller_params/
filedir = sys.argv[1]
suffix = '_mueller_matrix_from_inverted_params.txt'

#Currently set for the range of all 10B_036 with guppi
#May want to rebuild depending what sessions you want to grab.
directory = 'GBT10B_036'
directory2 = 'GBT11B_055'
#directory3 = 'GBT12A_418'
min_sess = 41
min_sess2 = 01
#min_sess3 = 01
max_sess = 89
max_sess2 = 18
#max_sess3 = 02
#Number of sessions used
len = max_sess - min_sess
len2 = max_sess2-min_sess2
#len3 = max_sess3-min_sess3
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
#for i in range(0,len3):
#   label=min_sess3+i
#   if label<10:
#      label ='0'+str(label)
#   else:
#      label = str(label)
#   filenames.append(filedir+directory3+'/'+label+suffix)
#   sessions.append(label)

#print filenames

#Need a test_dir to determine the number of frequency slices
array =[] 
test_dir = loadtxt(filenames[0])
scale = shape(test_dir)
#print array

#Database of mueller files
tot_len = len+len2
#tot_len = len+len2+len3
for j in range(0,tot_len):
   source = filenames[j]
   data = loadtxt(source)
   array.append(data)
#print array[0][0][1]
#print array[1][0][1]

#List of frequencies for axis.
freqs = test_dir[:,0]

#print freqs

m_II = sp.zeros((tot_len,scale[0]))
m_IQ = sp.zeros((tot_len,scale[0]))
m_IU = sp.zeros((tot_len,scale[0]))
m_IV = sp.zeros((tot_len,scale[0]))
m_QI = sp.zeros((tot_len,scale[0]))
m_QQ = sp.zeros((tot_len,scale[0]))
m_QU = sp.zeros((tot_len,scale[0]))
m_QV = sp.zeros((tot_len,scale[0]))
m_UI = sp.zeros((tot_len,scale[0]))
m_UQ = sp.zeros((tot_len,scale[0]))
m_UU = sp.zeros((tot_len,scale[0]))
m_UV = sp.zeros((tot_len,scale[0]))
m_VI = sp.zeros((tot_len,scale[0]))
m_VQ = sp.zeros((tot_len,scale[0]))
m_VU = sp.zeros((tot_len,scale[0]))
m_VV = sp.zeros((tot_len,scale[0]))

#Building separate tables for each matrix element.
for i in range(0,tot_len):
   for j in range(0,scale[0]):
      m_II[i,j] = array[i][j][1]
      #All terms other than m_II normalized w/ respect to m_II
      m_IQ[i,j] = array[i][j][2]/m_II[i,j]
      m_IU[i,j] = array[i][j][3]/m_II[i,j]
      m_IV[i,j] = array[i][j][4]/m_II[i,j]
      m_QI[i,j] = array[i][j][5]/m_II[i,j]
      m_QQ[i,j] = array[i][j][6]/m_II[i,j]
      m_QU[i,j] = array[i][j][7]/m_II[i,j]
      m_QV[i,j] = array[i][j][8]/m_II[i,j]
      m_UI[i,j] = array[i][j][9]/m_II[i,j]
      m_UQ[i,j] = array[i][j][10]/m_II[i,j]
      m_UU[i,j] = array[i][j][11]/m_II[i,j]
      m_UV[i,j] = array[i][j][12]/m_II[i,j]
      m_VI[i,j] = array[i][j][13]/m_II[i,j]
      m_VQ[i,j] = array[i][j][14]/m_II[i,j]
      m_VU[i,j] = array[i][j][15]/m_II[i,j]
      m_VV[i,j] = array[i][j][16]/m_II[i,j]

#print m_II

pylab.suptitle('Gain Correction Factor (From flux and pol cal generation)')
#M_II
#pylab.subplot(441)
#pylab.imshow(m_II/m_II,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.,vmax=1.,extent=(freqs.min(),freqs.max(),min_sess,max_sess))
pylab.imshow(m_II,interpolation='gaussian', cmap='hot', origin='lower',vmin=0.75,vmax=3.25,extent=(freqs.max(),freqs.min(),0,tot_len))
pylab.colorbar(orientation='horizontal')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Session')
pylab.savefig('Flux_'+str(min_sess)+'-'+str(max_sess)+'_and_'+str(min_sess2)+'-'+str(max_sess2)+'_compare.png')
pylab.clf()

#M_IQ
#pylab.subplot(442)
#pylab.imshow(m_IQ,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.5,vmax=0.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal')
#pylab.savefig('m_IQ_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf()

#M_IU 
#pylab.subplot(443)
#pylab.imshow(m_IU,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.06,vmax=0.06,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_IU_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_IV 
#pylab.subplot(444)
#pylab.imshow(m_IV,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.1,vmax=0.1,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_IV_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_QI 
#pylab.subplot(445)
#pylab.imshow(m_QI,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.5,vmax=0.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_QI_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_QQ
#pylab.subplot(446)
#pylab.imshow(m_QQ,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.5,vmax=1.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_QQ_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_QU
#pylab.subplot(447)
#pylab.imshow(m_QU,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.5,vmax=1.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_QU_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_QV
#pylab.subplot(448)
#pylab.imshow(m_QV,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.5,vmax=1.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_QV_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_UI 
#pylab.subplot(449)
#pylab.imshow(m_UI,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.5,vmax=0.5,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_UI_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_UQ
#pylab.subplot(4,4,10)
#pylab.imshow(m_UQ,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.25,vmax=1.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_UQ_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_UU 
#pylab.subplot(4,4,11)
#pylab.imshow(m_UU,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.25,vmax=1.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_UU_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_UV
#pylab.subplot(4,4,12)
#pylab.imshow(m_UV,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.25,vmax=1.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_UV_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_VI
#pylab.subplot(4,4,13)
#pylab.imshow(m_VI,interpolation='gaussian', cmap='hot', origin='lower',vmin=-0.25,vmax=0.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_VI_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_VQ 
#pylab.subplot(4,4,14)
#pylab.imshow(m_VQ,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1,vmax=1,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_VQ_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_VU 
#pylab.subplot(4,4,15)
#pylab.imshow(m_VU,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.25,vmax=1.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_VU_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#M_VV
#pylab.subplot(4,4,16)
#pylab.imshow(m_VV,interpolation='gaussian', cmap='hot', origin='lower',vmin=-1.25,vmax=1.25,extent=(freqs.max(),freqs.min(),0,tot_len))
#pylab.colorbar(orientation='horizontal') 
#pylab.savefig('m_VV_'+str(min_sess)+'-'+str(max_sess2)+'_compare.png')
#pylab.clf() 

#pylab.savefig('mueller_'+directory+'_'+str(min_sess)+'-'+str(max_sess)+'_compare.png')
#pylab.clf()
