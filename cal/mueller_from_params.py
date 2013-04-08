"""Script for generating mueller matrix elements from mueller params (generated using any script), expected params are (in order): frequency, deltaG, alpha, psi, phi, epsilon, chi, and flux.
Run from analysis_IM: python cal/mueller_from_params.py (note that it uses the file mueller_params_calc.txt that is currently in that folder). 
"""

import scipy as sp
import numpy.ma as ma
import pylab as pl
import numpy as np

# m_total is the final Mueller matrix that needs to multiply by the stokes parameters.

# Note: I have it set up so that the third index represents the frequency 
# bin where 0 is the lowest frequency bin and 7 is the highest, and the 
# first and second indices represent the mueller matrix for each frquency bin. 


mp = np.loadtxt('mueller_params_calc.txt')
#     print mp
#This is a file with first index being freq, second index being parameter where the values are:
# 0 = Freq, 1 = deltaG, 2 = alpha, 3 = psi, 4 = phi, 5 = epsilon, 6 = Q_src, 7 = U_src, 8 = chi

#Generates Inverse Mueller Matrix for use
freq_limit = len(mp[:,0])
#     print freq_limit
m_total = sp.zeros((4,4,freq_limit),float)
m_tot = sp.zeros((freq_limit,17),float)
for i in range(0,freq_limit):
#Note that this version is based on the combined matrix I derived corrected for linear feed:
    dG = mp[i,1]
    al = mp[i,2]*sp.pi/180
    ps = mp[i,3]*sp.pi/180
    ph = mp[i,4]*sp.pi/180
    ep = mp[i,5]
    ch = mp[i,6]*sp.pi/180
    flux = mp[i,7]     
    
    m_total[0,0,i] = flux
    m_total[0,1,i] = flux*(0.5*dG*sp.cos(2*al)-2*ep*sp.sin(2*al)*sp.cos(ph-ch))
    m_total[0,2,i] = flux*(0.5*dG*sp.sin(2*al)*sp.cos(ch)+2*ep*(sp.cos(al)*sp.cos(al)*sp.sin(ph)-sp.sin(al)*sp.sin(al)*sp.cos(ph-2*ch)))
    m_total[0,3,i] = flux*(0.5*dG*sp.sin(2*al)*sp.sin(ch)+2*ep*(sp.cos(al)*sp.cos(al)*sp.sin(ph)+sp.sin(al)*sp.sin(al)*sp.sin(ph-2*ch)))
    m_total[1,0,i] = flux*0.5*dG
    m_total[1,1,i] = flux*sp.cos(2*al)
    m_total[1,2,i] = flux*sp.sin(2*al)*sp.cos(ch)
    m_total[1,3,i] = flux*sp.sin(2*al)*sp.sin(ch)
    m_total[2,0,i] = flux*2*ep*sp.cos(ph+ps)
    m_total[2,1,i] = -flux*sp.sin(2*al)*sp.cos(ps+ch)
    m_total[2,2,i] = flux*(sp.cos(al)*sp.cos(al)*sp.cos(ps)-sp.sin(al)*sp.sin(al)*sp.cos(ps+2*ch))
    m_total[2,3,i] = flux*(-sp.cos(al)*sp.cos(al)*sp.sin(ps)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2*ch))
    m_total[3,0,i] = flux*2*ep*sp.sin(ps+ph)
    m_total[3,1,i] = -flux*sp.sin(2*al)*sp.sin(ps+ch)
    m_total[3,2,i] = flux*(sp.cos(al)*sp.cos(al)*sp.sin(ps)-sp.sin(al)*sp.sin(al)*sp.sin(ps+2*ch))
    m_total[3,3,i] = flux*(sp.cos(al)*sp.cos(al)*sp.cos(ps)+sp.sin(al)*sp.sin(al)*sp.cos(ps+2*ch))
    M_total = sp.mat(m_total[:,:,i])
#   print M_total
#    M_total = M_total.I
#   M_astron = sp.mat(m_astron[:,:,i])
#   M_total = M_astron*M_total
    m_tot[i,0] = mp[i,0]
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
        
prefix = '18'
#path = '$GBT10B_OUT/mueller_params/'
suffix = '_mueller_matrix_from_inverted_params.txt'
filename = prefix+suffix
np.savetxt(filename, m_tot[:,:], delimiter=' ')

