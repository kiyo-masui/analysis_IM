"""Script for generating mueller matrix elements from jones data"""

import scipy as sp
import numpy.ma as ma
import pylab as pl
import numpy as np

jones = np.loadtxt('Jones1019_2_fluxcal.txt')

# Generates Inverse Mueller Matrix for use
freq_limit = len(jones[:,0])
#     print freq_limit
m_total = sp.zeros((4,4,freq_limit),float)
m_tot = sp.zeros((freq_limit,17),float)
for i in range(0,freq_limit):
    m_total[0,0,i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]+jones[i,3]*jones[i,3]+jones[i,4]*jones[i,4]+jones[i,5]*jones[i,5]+jones[i,6]*jones[i,6]+jones[i,7]*jones[i,7]+jones[i,8]*jones[i,8])
    m_total[0,1,i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]-jones[i,3]*jones[i,3]-jones[i,4]*jones[i,4]+jones[i,5]*jones[i,5]+jones[i,6]*jones[i,6]-jones[i,7]*jones[i,7]-jones[i,8]*jones[i,8])
    m_total[0,2,i] = 0.5*(jones[i,1]*jones[i,3]+jones[i,2]*jones[i,4]+jones[i,5]*jones[i,7]+jones[i,6]*jones[i,8])
    m_total[0,3,i] = 0.5*(-jones[i,2]*jones[i,3]+jones[i,1]*jones[i,4]+jones[i,8]*jones[i,5]-jones[i,7]*jones[i,6])
    m_total[1,0,i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]+jones[i,3]*jones[i,3]+jones[i,4]*jones[i,4]-jones[i,5]*jones[i,5]-jones[i,6]*jones[i,6]-jones[i,7]*jones[i,7]-jones[i,8]*jones[i,8])
    m_total[1,1,i] = 0.25*(jones[i,1]*jones[i,1]+jones[i,2]*jones[i,2]-jones[i,3]*jones[i,3]-jones[i,4]*jones[i,4]-jones[i,5]*jones[i,5]-jones[i,6]*jones[i,6]+jones[i,7]*jones[i,7]+jones[i,8]*jones[i,8])
    m_total[1,2,i] = 0.5*(jones[i,1]*jones[i,3]+jones[i,2]*jones[i,4]-jones[i,5]*jones[i,7]-jones[i,6]*jones[i,8])
    m_total[1,3,i] = 0.5*(-jones[i,2]*jones[i,3]+jones[i,1]*jones[i,4]-jones[i,8]*jones[i,5]+jones[i,7]*jones[i,6])
    m_total[2,0,i] = 0.5*(jones[i,1]*jones[i,5]+jones[i,2]*jones[i,6]+jones[i,3]*jones[i,7]+jones[i,4]*jones[i,8])
    m_total[2,1,i] = 0.5*(jones[i,1]*jones[i,5]+jones[i,2]*jones[i,6]-jones[i,3]*jones[i,7]-jones[i,4]*jones[i,8])
    m_total[2,2,i] = 0.5*(jones[i,1]*jones[i,7]+jones[i,2]*jones[i,8]+jones[i,3]*jones[i,5]+jones[i,4]*jones[i,6])
    m_total[2,3,i] = 0.5*(-jones[i,2]*jones[i,7]+jones[i,1]*jones[i,8]+jones[i,4]*jones[i,5]-jones[i,3]*jones[i,6])
    m_total[3,0,i] = 0.5*(jones[i,2]*jones[i,5]-jones[i,1]*jones[i,6]+jones[i,4]*jones[i,7]-jones[i,3]*jones[i,8])
    m_total[3,1,i] = 0.5*(jones[i,2]*jones[i,5]-jones[i,1]*jones[i,6]-jones[i,4]*jones[i,7]+jones[i,3]*jones[i,8])
    m_total[3,2,i] = 0.5*(jones[i,2]*jones[i,7]-jones[i,1]*jones[i,8]+jones[i,4]*jones[i,5]-jones[i,3]*jones[i,6])
    m_total[3,3,i] = 0.5*(jones[i,1]*jones[i,7]+jones[i,2]*jones[i,8]-jones[i,3]*jones[i,5]-jones[i,4]*jones[i,6])

    M_total = sp.mat(m_total[:,:,i])
#   print M_total
    M_total = M_total.I
#   M_astron = sp.mat(m_astron[:,:,i])
#   M_total = M_astron*M_total
    m_tot[i,0] = jones[i,0]
    m_tot[i,1] = M_total[0,0]
    m_tot[i,2] = M_total[0,1]
    m_tot[i,3] = M_total[0,2]
    m_tot[i,4] = M_total[0,3]
    m_tot[i,5] = M_total[1,0]
    m_tot[i,6] = M_total[1,1]
    m_tot[i,7] = M_total[1,2]
    m_tot[i,8] = M_total[1,3]
    m_tot[i,9] = M_total[2,0]
    m_tot[i,10] = M_total[2,1]
    m_tot[i,11] = M_total[2,2]
    m_tot[i,12] = M_total[2,3]
    m_tot[i,13] = M_total[3,0]
    m_tot[i,14] = M_total[3,1]
    m_tot[i,15] = M_total[3,2]
    m_tot[i,16] = M_total[3,3]
        
prefix = '16'
#path = '$GBT10B_OUT/mueller_params/'
suffix = '_mueller_matrix_from_jones.txt'
filename = prefix+suffix
np.savetxt(filename, m_tot[:,:], delimiter=' ')

