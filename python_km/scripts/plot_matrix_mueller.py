import pylab
from numpy import *
import scipy as sp
import sys

prefix = '41-42'

mm_1_row = loadtxt('mueller_matrix_first_row.txt')
mm_2_row = loadtxt('mueller_matrix_second_row.txt')
mm_3_row = loadtxt('mueller_matrix_third_row.txt')
mm_4_row = loadtxt('mueller_matrix_fourth_row.txt')
source = loadtxt('mueller_params_calc.txt')
jones = loadtxt('jones_bin1.txt')
filename = sys.argv[1]
mm_dir = loadtxt(filename)

freq = source[:,0]
freq_j = sp.zeros(len(jones[:,0]))
mj_00 = sp.zeros(len(jones[:,0]))
mj_01 = sp.zeros(len(jones[:,0]))
mj_02 = sp.zeros(len(jones[:,0]))
mj_03 = sp.zeros(len(jones[:,0]))
mj_10 = sp.zeros(len(jones[:,0]))
mj_11 = sp.zeros(len(jones[:,0]))
mj_12 = sp.zeros(len(jones[:,0]))
mj_13 = sp.zeros(len(jones[:,0]))
mj_20 = sp.zeros(len(jones[:,0]))
mj_21 = sp.zeros(len(jones[:,0]))
mj_22 = sp.zeros(len(jones[:,0]))
mj_23 = sp.zeros(len(jones[:,0]))
mj_30 = sp.zeros(len(jones[:,0]))
mj_31 = sp.zeros(len(jones[:,0]))
mj_32 = sp.zeros(len(jones[:,0]))
mj_33 = sp.zeros(len(jones[:,0]))

for i in range(0,len(jones[:,0])):
    freq_j[i] = 900.0-jones[i,0]*200/jones[-1,0]
#    print freq_j[i]
    mj_00[i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]+jones[i,3]*jones[i,3]+jones[i,4]*jones[i,4]+jones[i,5]*jones[i,5]+jones[i,6]*jones[i,6]+jones[i,7]*jones[i,7]+jones[i,8]*jones[i,8])
    mj_01[i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]-jones[i,3]*jones[i,3]-jones[i,4]*jones[i,4]+jones[i,5]*jones[i,5]+jones[i,6]*jones[i,6]-jones[i,7]*jones[i,7]-jones[i,8]*jones[i,8])   
    mj_02[i] = 0.5*(jones[i,1]*jones[i,3]+jones[i,2]*jones[i,4]+jones[i,5]*jones[i,7]+jones[i,6]*jones[i,8])
    mj_03[i] = 0.5*(-jones[i,2]*jones[i,3]+jones[i,1]*jones[i,4]+jones[i,8]*jones[i,5]-jones[i,7]*jones[i,6])
    mj_10[i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]+jones[i,3]*jones[i,3]+jones[i,4]*jones[i,4]-jones[i,5]*jones[i,5]-jones[i,6]*jones[i,6]-jones[i,7]*jones[i,7]-jones[i,8]*jones[i,8])
    mj_11[i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]-jones[i,3]*jones[i,3]-jones[i,4]*jones[i,4]-jones[i,5]*jones[i,5]-jones[i,6]*jones[i,6]+jones[i,7]*jones[i,7]+jones[i,8]*jones[i,8])
    mj_12[i] = 0.5*(jones[i,1]*jones[i,3]+jones[i,2]*jones[i,4]-jones[i,5]*jones[i,7]-jones[i,6]*jones[i,8])
    mj_13[i] = 0.5*(-jones[i,2]*jones[i,3]+jones[i,1]*jones[i,4]-jones[i,8]*jones[i,5]+jones[i,7]*jones[i,6])
    mj_20[i] = 0.5*(jones[i,1]*jones[i,5]+jones[i,2]*jones[i,6]+jones[i,3]*jones[i,7]+jones[i,4]*jones[i,8])
    mj_21[i] = 0.5*(jones[i,1]*jones[i,5]+jones[i,2]*jones[i,6]-jones[i,3]*jones[i,7]-jones[i,4]*jones[i,8])
    mj_22[i] = 0.5*(jones[i,1]*jones[i,7]+jones[i,2]*jones[i,8]+jones[i,3]*jones[i,5]+jones[i,4]*jones[i,6])
    mj_23[i] = 0.5*(-jones[i,2]*jones[i,7]+jones[i,1]*jones[i,8]+jones[i,4]*jones[i,5]-jones[i,3]*jones[i,6])
    mj_30[i] = 0.5*(jones[i,2]*jones[i,5]-jones[i,1]*jones[i,6]+jones[i,4]*jones[i,7]-jones[i,3]*jones[i,8])
    mj_31[i] = 0.5*(jones[i,2]*jones[i,5]-jones[i,1]*jones[i,6]-jones[i,4]*jones[i,7]+jones[i,3]*jones[i,8]) 
    mj_32[i] = 0.5*(jones[i,2]*jones[i,7]-jones[i,1]*jones[i,8]+jones[i,4]*jones[i,5]-jones[i,3]*jones[i,6])
    mj_33[i] = 0.5*(jones[i,1]*jones[i,7]+jones[i,2]*jones[i,8]-jones[i,3]*jones[i,5]-jones[i,4]*jones[i,6]) 

mm_00 = mm_1_row[0]
mm_01 = mm_1_row[1]
mm_02 = mm_1_row[2]
mm_03 = mm_1_row[3]
mm_10 = mm_2_row[0]
mm_11 = mm_2_row[1]
mm_12 = mm_2_row[2]
mm_13 = mm_2_row[3]
mm_20 = mm_3_row[0]
mm_21 = mm_3_row[1]
mm_22 = mm_3_row[2]
mm_23 = mm_3_row[3]
mm_30 = mm_4_row[0]
mm_31 = mm_4_row[1]
mm_32 = mm_4_row[2]
mm_33 = mm_4_row[3]

#M00
pylab.plot(freq_j, mj_00, label='Jones')
pylab.plot(freq, mm_00, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,1], label='Direct')
pylab.ylim(0,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_II_compare_'+prefix+'.png')
pylab.clf()

#M01
pylab.plot(freq_j, mj_01, label='Jones')
pylab.plot(freq, mm_01, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,2], label='Direct')
pylab.ylim(-0.25,0.75)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_IQ_compare_'+prefix+'.png')
pylab.clf()

#M02
pylab.plot(freq_j, mj_02, label='Jones')
pylab.plot(freq, mm_02, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,3], label='Direct')
pylab.ylim(-0.35,0.35)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_IU_compare_'+prefix+'.png')
pylab.clf() 
 
#M03
pylab.plot(freq_j, mj_03, label='Jones')
pylab.plot(freq, mm_03, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,4], label='Direct')
pylab.ylim(-0.25,0.25)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_IV_compare_'+prefix+'.png')
pylab.clf() 

#M10
pylab.plot(freq_j, mj_10, label='Jones')
pylab.plot(freq, mm_10, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,5], label='Direct')
pylab.ylim(-0.35,0.35)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_QI_compare_'+prefix+'.png')
pylab.clf()
  
#M11
pylab.plot(freq_j, mj_11, label='Jones')
pylab.plot(freq, mm_11, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,6], label='Direct')
pylab.ylim(-1.5,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_QQ_compare_'+prefix+'.png')
pylab.clf()
  
#M12
pylab.plot(freq_j, mj_12, label='Jones')
pylab.plot(freq, mm_12, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,7], label='Direct')
pylab.ylim(-0.35,0.35)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_QU_compare_'+prefix+'.png')
pylab.clf()
  
#M13
pylab.plot(freq_j, mj_13, label='Jones')
pylab.plot(freq, mm_13, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,8], label='Direct')
pylab.ylim(-0.35,0.35)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_QV_compare_'+prefix+'.png')
pylab.clf()
  
#M20
pylab.plot(freq_j, mj_20, label='Jones')
pylab.plot(freq, mm_20, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,9], label='Direct')
pylab.ylim(-0.35,0.35)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_UI_compare_'+prefix+'.png')
pylab.clf()
  
#M21
pylab.plot(freq_j, mj_21, label='Jones')
pylab.plot(freq, mm_21, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,10], label='Direct')
pylab.ylim(-0.25,0.25)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_UQ_compare_'+prefix+'.png')
pylab.clf()
  
#M22
pylab.plot(freq_j, mj_22, label='Jones')
pylab.plot(freq, mm_22, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,11], label='Direct')
pylab.ylim(-1.5,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_UU_compare_'+prefix+'.png')
pylab.clf()
  
#M23
pylab.plot(freq_j, mj_23, label='Jones')
pylab.plot(freq, mm_23, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,12], label='Direct')
pylab.ylim(-1.5,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_UV_compare_'+prefix+'.png')
pylab.clf()
  
#M30
pylab.plot(freq_j, mj_30, label='Jones')
pylab.plot(freq, mm_30, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,13], label='Direct')
pylab.ylim(-0.25,0.25)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_VI_compare_'+prefix+'.png')
pylab.clf()
  
#M31
pylab.plot(freq_j, mj_31, label='Jones')
pylab.plot(freq, mm_31, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,14], label='Direct')
pylab.ylim(-0.25,0.25)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_VQ_compare_'+prefix+'.png')
pylab.clf()
  
#M32
pylab.plot(freq_j, mj_32, label='Jones')
pylab.plot(freq, mm_32, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,15], label='Direct')
pylab.ylim(-1.5,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_VU_compare_'+prefix+'.png')
pylab.clf()
  
#M33
pylab.plot(freq_j, mj_33, label='Jones')
pylab.plot(freq, mm_33, label='Mueller')
#pylab.plot(mm_dir[:,0],mm_dir[:,16], label='Direct')
pylab.ylim(-1.5,1.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_VV_compare_'+prefix+'.png')
pylab.clf()
